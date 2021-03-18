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
#include "asm_plotter.h"

static char *Usage[] = { " [-w<double(6.0)>] [-h<double(4.5)>]",
                         " [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]",
                         " [-lfs] [-pdf] [-z] [-T<int(4)>]",
                         " [-o<output>] <asm1>[.ktab] [<asm2>[.ktab]] <reads>[.ktab]"
                       };

static char template[15] = "._ASM.XXXX";

static void check_table(char *name)
{ static int KMER = 0;
  int   kmer;
  FILE *f;
  
  if (strcmp(name+(strlen(name)-5),".ktab") == 0)
    name = Catenate(name,".ktab","","");

  f = fopen(name,"r");
  if (f == NULL)
    { fprintf(stderr,"%s: Cannot find FastK table %s\n",Prog_Name,name);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (KMER == 0)
        KMER = kmer;
      else if (kmer != KMER)
        { fprintf(stderr,"%s: Kmer (%d) of table %s != %d\n",Prog_Name,kmer,name,KMER);
          exit (1);
        }
      fclose(f);
    }
}

int main(int argc, char *argv[])
{ int    LINE, FILL, STACK;
  int    PDF;
  int    ZGRAM;
  double XDIM, YDIM;
  double XREL, YREL;
  int    XMAX;
  int64  YMAX;
  char  *OUT;
  char  *ASM1, *ASM2;
  char  *READS;
  int    NTHREADS;
  
  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("ASMpLot");

    XDIM = 6.0;
    YDIM = 4.5;
    XREL = 2.1;
    YREL = 1.1;
    XMAX = 0;
    YMAX = 0;
    PDF  = 0;
    OUT  = NULL;
    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("lfsz")
            break;
          case 'h':
            ARG_REAL(YDIM);
            break;
          case 'o':
            OUT = argv[i]+2;
            break;
          case 'p':
            if (strcmp("df",argv[i]+2) == 0)
              PDF = 1;
            else
              { fprintf(stderr,"%s: don't recognize option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
            break;
          case 'w':
            ARG_REAL(XDIM);
            break;
          case 'x':
            ARG_REAL(XREL);
            if (XREL <= 0.)
              { fprintf(stderr,"%s: max x scaling factor must be > 0\n",Prog_Name);
                exit (1);
              }
            break;
          case 'y':
            ARG_REAL(YREL);
            if (YREL <= 0.)
              { fprintf(stderr,"%s: max y scaling factor must be > 0\n",Prog_Name);
                exit (1);
              }
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
          case 'X':
            ARG_POSITIVE(XMAX,"x max");
            break;
          case 'Y':
            { int ymax;

               ARG_POSITIVE(ymax,"y max");
               YMAX = ymax;
               break;
            }
        }
      else
        argv[j++] = argv[i];
    argc = j;

    LINE  = flags['l'];
    FILL  = flags['f'];
    STACK = flags['s'];
    ZGRAM = flags['z'];

    if (argc != 3 && argc != 4)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[3]);
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
	fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
	fprintf(stderr,"      -o: root name for output plots\n");
	fprintf(stderr,"          default is root path of <asm> argument\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        exit (1);
      }

    if (LINE+FILL+STACK == 0)
      LINE = FILL = STACK = 1;
    if (OUT == NULL)
      OUT = Root(argv[1],".ktab");
    ASM1 = argv[1];
    check_table(ASM1);
    if (argc == 4)
      { ASM2  = argv[2];
        check_table(ASM2);
        READS = argv[3];
      }
    else
      { ASM2  = NULL;
        READS = argv[2];
      }
    check_table(READS);
  }

  { char *troot;

    troot = mktemp(template);

    asmplot(OUT,ASM1,ASM2,READS,XDIM,YDIM,XREL,YREL,XMAX,YMAX,
            PDF,ZGRAM,LINE,FILL,STACK,troot,NTHREADS);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
