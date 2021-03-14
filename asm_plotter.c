/********************************************************************************************
 *
 *  Command line utility to produce Assembly-spectra plots
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

void asmplot(char  *OUT, char  *ASM1, char *ASM2, char  *READS,
              double XDIM, double YDIM,
              double XREL, double YREL,
              int    XMAX, int64  YMAX,
              int    PDF, int ZGRAM, int LINE, int FILL, int STACK,
              char  *troot, int NTHREADS)

{ char      command[1000];
  char      what[1000];
  char      extra[1000];
  Histogram *H[6];
  int        i, k, nhist;
  char      *Label[] = { "read-only", "", "", "shared" };
  FILE      *f;

  //  Call Logex to make desired histograms

  ZGRAM = (ZGRAM != 0);

  if (ASM2 == NULL)
    { sprintf(what,"'%s.0=B-A' '%s.1=B&.A'",troot,troot);
      if (ZGRAM)
        sprintf(extra," '%s.2=A-B'",troot);
      else
        sprintf(extra,"");
      sprintf(command,"Logex -T%d -H32767 %s%s %s %s",NTHREADS,extra,what,ASM1,READS);
      nhist = 2;
      Label[1] = Root(ASM1,".ktab");
      Label[2] = Root(ASM1,".ktab");
    }
  else
    { sprintf(what,"'%s.0=C-#(A|B)' '%s.1=C&.(A-B)' '%s.2=C&.(B-A)' '%s.3=C&.#(A&B)'",
                   troot,troot,troot,troot);
      if (ZGRAM)
        sprintf(extra," '%s.4=(A|+B)-C'",troot);
      else
        sprintf(extra,"");
      sprintf(command,"Logex -T%d -H32767 %s%s %s %s %s",NTHREADS,what,extra,ASM1,ASM2,READS);
      Label[1] = Root(ASM1,".ktab");
      Label[2] = Root(ASM2,".ktab");
      nhist = 4;
    }
#ifdef DEBUG
  printf("%s\n",command);
  fflush(stdout);
#endif
  system(command);

  //  Open all histograms

  for (i = 0; i < nhist+ZGRAM; i++)
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
          for (i = 0; i < nhist; i++)
            last += H[i]->hist[low];
          for (k = low+1; k < high; k++) 
            { sum = 0;
              for (i = 0; i < nhist; i++)
                sum += H[i]->hist[k];
              if (sum > last)
                break;
              last = sum;
            }
          ymax = sum;
          xmax = k;
          for ( ; k < high; k++)
            { sum = 0;
              for (i = 0; i < nhist; i++)
                sum += H[i]->hist[k];
              if (sum >= ymax)
                { ymax = sum;
                  xmax = k;
                }
            }
        }
      if (FILL+LINE > 0)
        for (i = 0; i < nhist; i++)
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
fflush(stdout);
#endif
    }
#ifdef DEBUG
   printf("x,y-max  = %d, %lld\n",XMAX,YMAX);
#endif

  //  Merge histograms into 1 for R plotter

  f = fopen(Catenate(troot,".asmi","",""),"w");
#ifdef DEBUG
  if (f == NULL)
    printf("Could not open %s\n",Catenate(troot,".asmi","",""));
  else
    printf("Writing %s\n",Catenate(troot,".asmi","",""));
  fflush(stdout);
#endif
  fprintf(f,"Copies\tkmer_multiplicity\tCount\n");
  for (i = 0; i < nhist; i++)
    { int    low  = H[i]->low;
      int    high = H[i]->high;
      int64 *hist = H[i]->hist;

      if (nhist == 4 && (i == 1 || i == 2))
        { for (k = low; k < high; k++)
            if (hist[k] > 0)
              fprintf(f,"%s-only\t%d\t%lld\n",Label[i],k,hist[k]);
        }
      else
        { for (k = low; k < high; k++)
            if (hist[k] > 0)
              fprintf(f,"%s\t%d\t%lld\n",Label[i],k,hist[k]);
        }
    }
  free(Label[1]);
  free(Label[2]);
  fclose(f);

  //  If -z option then use A-B histogram to produce needed data file

  if (ZGRAM)
    { int    high = H[nhist]->high;
      int64 *hist = H[nhist]->hist;
      int64  sum;

      f = fopen(Catenate(troot,".asmz","",""),"w");
#ifdef DEBUG
      if (f == NULL)
        printf("Could not open %s\n",Catenate(troot,".asmz","",""));
      else
        printf("Writing %s\n",Catenate(troot,".asmz","",""));
      fflush(stdout);
#endif
      fprintf(f,"1\t%lld\n",hist[1]);
      sum = 0;
      for (k = 2; k <= high; k++)
        sum += hist[k];
      fprintf(f,"2\t%lld\n",sum);
      fclose(f);
    }

  sprintf(command,"Fastrm %s.*.hist",troot);
  system(command);

  //  Call the plotter with arguments

  sprintf(what,"plot_spectra_cn.R -f %s.asmi -o %s%s -x %g -y %g -m%d -n %lld",
               troot,OUT,PDF?" -p":" ",XDIM,YDIM,XMAX,YMAX);
  sprintf(extra," -z %s.asmz",troot);
  if (LINE+FILL+STACK == 3)
    { sprintf(command,"%s%s 2>/tmp/NULL",what,extra);
#ifdef DEBUG
      printf("%s\n",command);
      fflush(stdout);
#endif
      system(command);
    }
  else
    { if (LINE)
        { sprintf(command,"%s -t line%s 2>/tmp/NULL",what,extra);
#ifdef DEBUG
          printf("%s\n",command);
          fflush(stdout);
#endif
          system(command);
        }
      if (FILL)
        { sprintf(command,"%s -t fill%s 2>/tmp/NULL",what,extra);
#ifdef DEBUG
          printf("%s\n",command);
          fflush(stdout);
#endif
          system(command);
        }
      if (STACK)
        { sprintf(command,"%s -t stack%s 2>/tmp/NULL",what,extra);
#ifdef DEBUG
          printf("%s\n",command);
          fflush(stdout);
#endif
          system(command);
        }
    }

  sprintf(command,"rm -f %s.asmi %s.asmz",troot,troot);
  system(command);
}
