/********************************************************************************************
 *
 *  Command line utility to produce CN-spectra plots
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

#include "libfastk.h"
#include "cn_plotter.h"

static char *Usage[] = { " [-w<double(6.0)>] [-h<double(4.5)>]",
                         " [-[xX]<number(x2.1)>] [-[yY]<number(y1.1)>]",
                         " [-v] [-lfs] [-pdf] [-z] [-T<int(4)>] [-P<dir(/tmp)>]",
                         " <reads>[.ktab] <asm>:.dna> <out>"
                       };

static char template[15] = "._CN.XXXX";

static int check_table(char *name, int lmer)
{ int   kmer;
  FILE *f;

  f = fopen(name,"r");
  if (f == NULL)
    { fprintf(stderr,"%s: Cannot find FastK table %s\n",Prog_Name,name);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (lmer != 0 && kmer != lmer)
        { fprintf(stderr,"%s: Kmer (%d) of table %s != %d\n",Prog_Name,kmer,name,lmer);
          exit (1);
        }
      fclose(f);
      return (kmer);
    }
}

int main(int argc, char *argv[])
{ int    KMER;
  int    VERBOSE;
  int    LINE, FILL, STACK;
  int    PDF;
  int    ZGRAM;
  double XDIM, YDIM;
  double XREL, YREL;
  int    XMAX;
  int64  YMAX;
  char  *OUT;
  char  *ASM;
  char  *READS;
  int    NTHREADS;
  char  *SORT_PATH;
  
  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("CNpLot");

    XDIM = 6.0;
    YDIM = 4.5;
    XREL = 2.1;
    YREL = 1.1;
    XMAX = 0;
    YMAX = 0;
    PDF  = 0;
    NTHREADS = 4;
    SORT_PATH = "/tmp";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vlfsz")
            break;
          case 'h':
            ARG_REAL(YDIM);
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
          case 'P':
            SORT_PATH = argv[i]+2;
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

    VERBOSE = flags['v'];
    LINE    = flags['l'];
    FILL    = flags['f'];
    STACK   = flags['s'];
    ZGRAM   = flags['z'];

    if (argc != 4)
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
        fprintf(stderr,"      -z: plot counts of k-mers unique to assembly\n");
        fprintf(stderr,"\n");
	fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose output to stderr\n");
	fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
      }

    if (LINE+FILL+STACK == 0)
      LINE = FILL = STACK = 1;

    READS = argv[1];
    ASM   = argv[2];
    OUT   = argv[3];
  }

  { char *suffix[9] = { ".gz", ".fa", ".fq", ".fasta", ".fastq", ".db", ".sam", ".bam", ".cram" };
    int   j, len;

    READS = Root(READS,".ktab");

    KMER = check_table(Catenate(READS,".ktab","",""),0);

    for (j = 0; j < 9; j++)
      { len = strlen(ASM) - strlen(suffix[j]);
        if (strcmp(ASM+len,suffix[j]) == 0)
          ASM[len] = '\0';
      }
  }

  { char *troot;
    char  command[5000];

    troot = mktemp(template);

    if (VERBOSE)
      fprintf(stderr,"\n Making k-mer table for assembly %s\n",ASM);

    sprintf(command,"FastK -k%d -T%d -P%s -t1 %s",KMER,NTHREADS,SORT_PATH,ASM);
    system(command);

    if (VERBOSE)
      fprintf(stderr,"\n Making spectra histograms and plotting\n");

    cn_plot(OUT,ASM,READS,XDIM,YDIM,XREL,YREL,XMAX,YMAX,PDF,ZGRAM,LINE,FILL,STACK,troot,NTHREADS);

    sprintf(command,"Fastrm %s",ASM);
    system(command);

  }

  free(READS);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
