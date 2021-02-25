/********************************************************************************************
 *
 *  Example code for opening and fetching compressed profiles produced by FastK
 *
 *  Author:  Gene Myers
 *  Date  :  October, 2020
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

#include "cn.R.h"

static char *Usage[] = { " [-w<double(6.0)>] [-h<double(4.5)>] [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]",
                         " [-lfs] [-P] [-z] [-o<output>] <asm>[.ktab] <reads>[.ktab]"
                       };

static char template[15] = "._CN.XXXX";

/****************************************************************************************
 *
 *  Main Routine
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ int    LINE, FILL, STACK;
  int    PDF;
  int    ZGRAM;
  double XDIM, YDIM;
  double XREL, YREL;
  int    XMAX;
  int64  YMAX;
  char  *OUT;
  
  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("CNpLot");

    XDIM = 6.0;
    YDIM = 4.5;
    XREL = 2.1;
    YREL = 1.1;
    XMAX = 0;
    YMAX = 0;
    OUT  = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("lfsPz")
            break;
          case 'o':
            OUT = argv[i]+2;
            break;
          case 'X':
            ARG_POSITIVE(XMAX,"x max");
            break;
          case 'x':
            ARG_REAL(XREL);
            if (XREL <= 0.)
              { fprintf(stderr,"%s: max x scaling factor must be > 0\n",Prog_Name);
                exit (1);
              }
            break;
          case 'Y':
            { int ymax;

               ARG_POSITIVE(ymax,"y max");
               YMAX = ymax;
               break;
            }
          case 'y':
            ARG_REAL(YREL);
            if (YREL <= 0.)
              { fprintf(stderr,"%s: max y scaling factor must be > 0\n",Prog_Name);
                exit (1);
              }
            break;
          case 'w':
            ARG_REAL(XDIM);
            break;
          case 'h':
            ARG_REAL(YDIM);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    PDF   = flags['P'];
    LINE  = flags['l'];
    FILL  = flags['f'];
    STACK = flags['s'];
    ZGRAM = flags['z'];

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -w: width in inches of plots\n");
        fprintf(stderr,"      -h: height in inches of plots\n");
        fprintf(stderr,"      -x: max x as a real-valued multiple of x* with max\n");
        fprintf(stderr,"              count 'peak' away from the origin\n");
        fprintf(stderr,"      -X: max x as an int value in absolute terms\n");
        fprintf(stderr,"      -y: max y as a real-valued multiple of max count\n");
        fprintf(stderr,"              'peak' away from the origin\n");
        fprintf(stderr,"      -Y: max y as an int value in absolute terms\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: draw line plot\n");
        fprintf(stderr,"      -f: draw fill plot\n");
        fprintf(stderr,"      -s: draw stack plot\n");
        fprintf(stderr,"          any combo allowed, none => draw all\n");
        fprintf(stderr,"\n");
	fprintf(stderr,"      -P: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
	fprintf(stderr,"      -o: root name for output plots\n");
	fprintf(stderr,"          default is root path of <asm> argument\n");
        exit (1);
      }

    if (LINE+FILL+STACK == 0)
      LINE = FILL = STACK = 1;
    if (OUT == NULL)
      OUT = Root(argv[1],".ktab");
  }

  //  Open assembly file ASM, assembly profile AP, and all reference profiles RP[i]

  { char      command[1000];
    char      what[1000];
    char     *troot;
    Histogram *H[6];
    int        i, k;
    char      *Label[] = { "read-only", "1", "2", "3", "4", ">4" };
    FILE      *f;

    troot = mktemp(template);

    //  Call Logex to make desired histograms

    sprintf(what,"'%s.0=B-A' '%s.1=B&.A[1]' '%s.2=B&.A[2]' '%s.3=B&.A[3]' '%s.4=B&.A[4]' '%s.5=B&.A[5-]'",
                 troot,troot,troot,troot,troot,troot);
    sprintf(command,"Logex -T6 -H32767 %s %s %s",what,argv[1],argv[2]);
#ifdef DEBUG
    printf("%s\n",command);
    fflush(stdout);
#endif
    system(command);

    //  Open all histograms

    for (i = 0; i <= 5; i++)
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

    f = fopen(Catenate(troot,".histo","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".histo","",""));
    else
      printf("Writing %s\n",Catenate(troot,".histo","",""));
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

    sprintf(command,"Fastrm %s.*.hist",troot);
    system(command);

    //  Generate the plot script

    f = fopen(Catenate(troot,".R","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".R","",""));
    else
      printf("Generating %s\n",Catenate(troot,".R","",""));
    fflush(stdout);
#endif
    fwrite(cn_plot,strlen(cn_plot),1,f);
    fclose(f);

    //  Call the plotter with arguments

    if (LINE+FILL+STACK == 3)
      { sprintf(command,"RScript %s.R -f %s.histo -o %s -x %g -y %g -m %d -n %lld %s 2>/tmp/NULL",
                        troot,troot,OUT,XDIM,YDIM,XMAX,YMAX,PDF?"-p":"");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        system(command);
      }
    else
      { if (LINE)
          { sprintf(command,"RScript %s.R -f %s.histo -o %s -t line -x %g -y %g -m %d -n %lld%s",
                            troot,troot,OUT,XDIM,YDIM,XMAX,YMAX,PDF?" -p":"");
#ifdef DEBUG
            printf("%s\n",command);
            fflush(stdout);
#endif
            system(command);
          }
        if (FILL)
          { sprintf(command,"RScript %s.R -f %s.histo -o %s -t fill -x %g -y %g -m %d -n %lld%s",
                            troot,troot,OUT,XDIM,YDIM,XMAX,YMAX,PDF?" -p":"");
#ifdef DEBUG
            printf("%s\n",command);
            fflush(stdout);
#endif
            system(command);
          }
        if (STACK)
          { sprintf(command,"RScript %s.R -f %s.histo -o %s -t stack -x %g -y %g -m %d -n %lld%s",
                            troot,troot,OUT,XDIM,YDIM,XMAX,YMAX,PDF?" -p":"");
#ifdef DEBUG
            printf("%s\n",command);
            fflush(stdout);
#endif
            system(command);
          }
      }

    sprintf(command,"rm %s.histo %s.R",troot,troot);
    system(command);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
