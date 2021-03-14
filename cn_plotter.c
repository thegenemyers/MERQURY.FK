/********************************************************************************************
 *
 *  Routine to produce CN-spectra plots
 *
 *  Author:  Gene Myers
 *  Date  :  March, 2021
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#define DEBUG

#include "libfastk.h"

int64 cnplot(char  *OUT, char  *ASM, char  *READS,
             double XDIM, double YDIM,
             double XREL, double YREL,
             int    XMAX, int64  YMAX,
             int    PDF, int ZGRAM, int LINE, int FILL, int STACK,
             char  *troot, int NTHREADS)
        
{ char      command[1000];  //  buffers for commands & parts thereof
  char      what[1000];
  char      extra[1000];
  Histogram *H[7];
  int        i, k;
  int64      errors;
  char      *Label[] = { "read-only", "1", "2", "3", "4", ">4" };
  FILE      *f;

  ZGRAM = (ZGRAM != 0);

  //  Call Logex to make desired histograms

  sprintf(what,"'%s.0=B-A' '%s.1=B&.A[1]' '%s.2=B&.A[2]' '%s.3=B&.A[3]' '%s.4=B&.A[4]'",
               troot,troot,troot,troot,troot);
  if (ZGRAM)
    sprintf(extra," '%s.5=B&.A[5-]' '%s.6=A-B'",troot,troot);
  else
    sprintf(extra," '%s.5=B&.A[5-]'",troot);
  sprintf(command,"Logex -T%d -H1000 %s%s %s %s",NTHREADS,what,extra,ASM,READS);
#ifdef DEBUG
  printf("%s\n",command);
  fflush(stdout);
#endif
  system(command);

  //  Open all histograms

  for (i = 0; i <= 5+ZGRAM; i++)
    { sprintf(command,"%s.%d",troot,i);
      H[i] = Load_Histogram(command);
      Modify_Histogram(H[i],H[i]->low,H[i]->high,0);
    }

  //  If relative x- or y-max then must find peak x,y

  if (XMAX == 0 || YMAX == 0)
    { int64 sum, last;
      int64 ym, ymax;
      int   xm, xmax;   

      ymax = xmax = -1;
      if (STACK)
        { int  low  = H[0]->low;
          int  high = H[0]->high;

          last = 0;
          for (i = 0; i <= 5; i++)
            last += H[i]->hist[low];
          for (k = low+1; k < high; k++) 
            { sum = 0;
              for (i = 0; i <= 5; i++)
                sum += H[i]->hist[k];
              if (sum > last)
                break;
              last = sum;
            }
          ymax = sum;
          xmax = k;
          for ( ; k < high; k++)
            { sum = 0;
              for (i = 0; i <= 5; i++)
                sum += H[i]->hist[k];
              if (sum >= ymax)
                { ymax = sum;
                  xmax = k;
                }
            }
        }
      if (FILL+LINE > 0)
        for (i = 0; i <= 5; i++)
          { int    low  = H[i]->low;
            int    high = H[i]->high;
            int64 *hist = H[i]->hist;

            for (k = low+1; hist[k] < hist[k-1]; k++) 
              if (k >= high)
                break;
            ym = hist[k];
            xm = k;
            for ( ; k < high; k++)
              if (hist[k] >= ym)
                { ym = hist[k];
                  xm = k;
                }
            if (ym > ymax)
              { ymax = ym;
                xmax = xm;
              }
          }
       if (XMAX == 0)
         XMAX = xmax*XREL;
       if (YMAX == 0)
         YMAX = ymax*YREL;
#ifdef DEBUG
       printf("x,y-peak = %d, %lld\n",xmax,ymax);
#endif
    }
#ifdef DEBUG
   printf("x,y-max  = %d, %lld\n",XMAX,YMAX);
#endif

  //  Merge histograms into 1 for R plotter

  f = fopen(Catenate(troot,".cni","",""),"w");
#ifdef DEBUG
  if (f == NULL)
    printf("Could not open %s\n",Catenate(troot,".cni","",""));
  else
    printf("Writing %s\n",Catenate(troot,".cni","",""));
  fflush(stdout);
#endif
  fprintf(f,"Copies\tkmer_multiplicity\tCount\n");
  for (i = 0; i <= 5; i++)
    { int    low  = H[i]->low;
      int    high = H[i]->high;
      int64 *hist = H[i]->hist;

      for (k = low; k < high; k++)
        if (hist[k] > 0)
          fprintf(f,"%s\t%d\t%lld\n",Label[i],k,hist[k]);
    }
  fclose(f);

  //  If -z option then use A-B histogram to produce needed data file

  if (ZGRAM)
    { int    high = H[6]->high;
      int64 *hist = H[6]->hist;

      f = fopen(Catenate(troot,".cnz","",""),"w");
#ifdef DEBUG
      if (f == NULL)
        printf("Could not open %s\n",Catenate(troot,".cnz","",""));
      else
        printf("Writing %s\n",Catenate(troot,".cnz","",""));
      fflush(stdout);
#endif
      fprintf(f,"1\t%lld\n",hist[1]);
      errors = 0;
      for (k = 2; k <= high; k++)
        errors += hist[k];
      fprintf(f,"2\t%lld\n",errors);
      fclose(f);

      errors += hist[1];
    }
  else
    errors = 0;

  //  Remove all the FastK histograms except A-B

  sprintf(command,"Fastrm %s.*.hist",troot);
  system(command);

  //  Call the R plotter with arguments

  sprintf(what,"plot_spectra_cn.R -f %s.cni -o %s%s -x %g -y %g -m%d -n %lld",
               troot,OUT,PDF?" -p":" ",XDIM,YDIM,XMAX,YMAX);
  sprintf(extra," -z %s.cnz",troot);
  if (LINE+FILL+STACK == 3)
    { sprintf(command,"%s%s 2>/tmp/NULL",what,ZGRAM?extra:"");
#ifdef DEBUG
      printf("%s\n",command);
      fflush(stdout);
#endif
      system(command);
    }
  else
    { if (LINE)
        { sprintf(command,"%s -t line%s 2>/tmp/NULL",what,ZGRAM?extra:"");
#ifdef DEBUG
          printf("%s\n",command);
          fflush(stdout);
#endif
          system(command);
        }
      if (FILL)
        { sprintf(command,"%s -t fill%s 2>/tmp/NULL",what,ZGRAM?extra:"");
#ifdef DEBUG
          printf("%s\n",command);
          fflush(stdout);
#endif
          system(command);
        }
      if (STACK)
        { sprintf(command,"%s -t stack%s 2>/tmp/NULL",what,ZGRAM?extra:"");
#ifdef DEBUG
          printf("%s\n",command);
          fflush(stdout);
#endif
          system(command);
        }
    }

  sprintf(command,"rm -f %s.cni %s.cnz",troot,troot);
  system(command);

  return (errors);
}
