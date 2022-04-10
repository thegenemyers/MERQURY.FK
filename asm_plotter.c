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

#undef DEBUG

#include "libfastk.h"

#include "cn_plot.R.h"

static int64 GetCount(Histogram *H)
{ int64 sum;
  int   k;

  sum = 0;
  for (k = H->low; k <= H->high; k++)
    sum += H->hist[k];
  return (sum);
}

void asm_plot(char  *OUT, char  *ASM1, char *ASM2, char  *READS,
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
      sprintf(command,"Logex -T%d -H1000 %s%s %s %s",NTHREADS,what,extra,ASM1,READS);
      nhist = 2;
      Label[1] = Root(ASM1,".ktab");
      Label[2] = Root(ASM1,".ktab");
    }
  else
    { sprintf(what,"'%s.0=C-#(A|B)' '%s.1=C&.(A-B)' '%s.2=C&.(B-A)' '%s.3=C&.#(A&B)'",
                   troot,troot,troot,troot);
      if (ZGRAM)
        sprintf(extra," '%s.4=(A&+B)-C' '%s.5=(A-B)-C' '%s.6=(B-A)-C'",troot,troot,troot);
      else
        sprintf(extra,"");
      sprintf(command,"Logex -T%d -H1000 %s%s %s %s %s",NTHREADS,what,extra,ASM1,ASM2,READS);
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

  for (i = 0; i < nhist*(1+ZGRAM)-ZGRAM; i++)
    { sprintf(command,"%s.%d",troot,i);
      H[i] = Load_Histogram(command);
    }

  //  If relative x- or y-max then must find peak x,y

  if (XMAX == 0 || YMAX == 0)
    { int64 sum, last;
      int64 ym, ymax;
      int   xm, xmax;   
      int       xsec;

      ymax = xmax = -1;
      if (STACK)
        { int  low  = H[0]->low;
          int  high = H[0]->high;

          last = 0;                     //  Find highest y away from origin
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

          xsec = xmax;                 //  Find any seconday peak 10% of y-max
          last = 0;
          for (i = 0; i < nhist; i++)
            last += H[i]->hist[xmax];
          k = xmax+1;
          sum = 0;
          for (i = 0; i < nhist; i++)
            sum += H[i]->hist[k];
          while (k < high)
            { while (sum <= last)
                { if (k >= high)
                    break;
                  k += 1;
                  last = sum;
                  sum = 0;
                  for (i = 0; i < nhist; i++)
                    sum += H[i]->hist[k];
                }
              while (sum >= last)
                { if (k >= high)
                    break;
                  k += 1;
                  last = sum;
                  sum = 0;
                  for (i = 0; i < nhist; i++)
                    sum += H[i]->hist[k];
                }
              if (last >= .1*ymax)
                xsec = k-1;
            }
        }

      else //  FILL || LINE
        { for (i = 0; i < nhist; i++)
            { int    low  = H[i]->low;
              int    high = H[i]->high;
              int64 *hist = H[i]->hist;

              for (k = low+1; hist[k] < hist[k-1]; k++)    //  Find largest y away from 0
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

          xsec = xmax;                     //  Find highest y away from origin
          for (i = 0; i < nhist; i++)
            { int    high = H[i]->high;
              int64 *hist = H[i]->hist;

              k = xmax+1;
              while (k < high)
                { while (hist[k] <= hist[k-1])
                    { if (k >= high)
                        break;
                      k += 1;
                    }
                  while (hist[k] >= hist[k-1])
                    { if (k >= high)
                        break;
                      k += 1;
                    }
                  if (hist[k-1] >= .1*ymax)
                    xsec = k-1;
                }
            }
        }

      if (XMAX == 0)
        XMAX = ((xmax+xsec)*XREL)/2;
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
  fclose(f);

  //  If -z option then use A-B histogram to produce needed data file

  if (ZGRAM)
    { f = fopen(Catenate(troot,".asmz","",""),"w");
#ifdef DEBUG
      if (f == NULL)
        printf("Could not open %s\n",Catenate(troot,".asmz","",""));
      else
        printf("Writing %s\n",Catenate(troot,".asmz","",""));
      fflush(stdout);
#endif

      if (nhist == 4)
        { fprintf(f,"%s-only\t0\t%lld\n",ASM1,GetCount(H[nhist]));
          fprintf(f,"%s-only\t0\t%lld\n",ASM2,GetCount(H[nhist+1]));
          fprintf(f,"shared\t0\t%lld\n",GetCount(H[nhist+2]));
        }
      else
        fprintf(f,"%s\t0\t%lld\n",ASM1,GetCount(H[nhist]));
      fclose(f);
    }

  sprintf(command,"Fastrm %s.*.hist",troot);
  system(command);

  free(Label[1]);
  free(Label[2]);

  //  Generate the R plot script

  f = fopen(Catenate(troot,".R","",""),"w");
#ifdef DEBUG
  if (f == NULL)
    printf("Could not open %s\n",Catenate(troot,".R","",""));
  else
    printf("Generating %s\n",Catenate(troot,".R","",""));
  fflush(stdout);
#endif
  fwrite(cn_plot_script,strlen(cn_plot_script),1,f);
  fclose(f);

  //  Call the plotter with arguments

  sprintf(what,"Rscript %s.R -f %s.asmi -o %s%s -x %g -y %g -m%d -n %lld",
               troot,troot,OUT,PDF?" -p":"",XDIM,YDIM,XMAX,YMAX);
  if (ZGRAM)
    sprintf(extra," -z %s.asmz",troot);
  else
    sprintf(extra,"");
  if (LINE)
    { sprintf(command,"%s -t line%s 2>/dev/null",what,extra);
#ifdef DEBUG
      printf("%s\n",command);
      fflush(stdout);
#endif
      system(command);
    }
  if (FILL)
    { sprintf(command,"%s -t fill%s 2>/dev/null",what,extra);
#ifdef DEBUG
      printf("%s\n",command);
      fflush(stdout);
#endif
      system(command);
    }
  if (STACK)
    { sprintf(command,"%s -t stack%s 2>/dev/null",what,extra);
#ifdef DEBUG
      printf("%s\n",command);
      fflush(stdout);
#endif
      system(command);
    }

  sprintf(command,"rm -f %s.asmi %s.asmz %s.R",troot,troot,troot);
  system(command);
}
