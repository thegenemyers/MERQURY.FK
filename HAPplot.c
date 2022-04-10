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
#include "hap_plotter.h"

static char *Usage[] =
    { " [-v] [-w<double(6.0)>] [-h<double(4.5)>] [-pdf] [-T<int(4)>] [-P<dir(/tmp)>]",
      " <mat>[.hap[.ktab]] <pat>[.hap[.ktab]] <asm1:dna> [<asm2:dna>] <out>"
    };

static char template[16] = "._HAP.XXXX";

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
  int    PDF;
  int    VERBOSE;
  double XDIM, YDIM;
  char  *OUT;
  char  *MAT;
  char  *PAT;
  char  *ASM[2];
  int    NTHREADS;
  char  *SORT_PATH;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("HAPpLot");

    XDIM = 6.0;
    YDIM = 4.5;
    PDF  = 0;
    NTHREADS = 4;
    SORT_PATH = "/tmp";

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
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
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 5 || argc > 6)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -w: width in inches of plots\n");
        fprintf(stderr,"      -h: height in inches of plots\n");
        fprintf(stderr,"\n");
	fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose output to stderr\n");
	fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        exit (1);
      }

    switch (argc)
      { case 5:
        if (VERBOSE)
          fprintf(stderr,"\n Single diploid assembly\n");
        ASM[0] = argv[3];
        ASM[1] = NULL;
        break;
      case 6:
        if (VERBOSE)
          fprintf(stderr,"\n Two haploid assemblies\n");
        ASM[0] = argv[3];
        ASM[1] = argv[4];
        break;
      default:
        ;
    }
    MAT = argv[1];
    PAT = argv[2];
    OUT = argv[argc-1];
  }

  { char *troot, *x;
    int   i, j, len;
    char *suffix[9] = { ".gz", ".fa", ".fq", ".fasta", ".fastq", ".db", ".sam", ".bam", ".cram" };
    char  command[5000];

    troot = mktemp(template);

    x = Root(MAT,".ktab");
    MAT = Root(x,".hap");
    free(x);
    x = Root(PAT,".ktab");
    PAT = Root(x,".hap");
    free(x);

    KMER = check_table(Catenate(MAT,".hap",".ktab",""),0);
    KMER = check_table(Catenate(PAT,".hap",".ktab",""),KMER);

    for (i = 0; i < 2; i++)
      { if (ASM[i] == NULL)
          continue;

        for (j = 0; j < 9; j++)
          { len = strlen(ASM[i]) - strlen(suffix[j]);
            if (strcmp(ASM[i]+len,suffix[j]) == 0)
              ASM[i][len] = '\0';
          }

        if (VERBOSE)
          fprintf(stderr,"\n Computing k-table and profiles for assembly %s\n",ASM[i]);

        sprintf(command,"FastK -k%d -T%d -P%s -p %s",KMER,NTHREADS,SORT_PATH,ASM[i]);
        system(command);
        sprintf(command,"FastK -k%d -T%d -P%s -p:%s.hap %s -N%s.%s",
                        KMER,NTHREADS,SORT_PATH,MAT,ASM[i],ASM[i],MAT);
        system(command);
        sprintf(command,"FastK -k%d -T%d -P%s -p:%s.hap %s -N%s.%s",
                        KMER,NTHREADS,SORT_PATH,PAT,ASM[i],ASM[i],PAT);
        system(command);
      }

    if (VERBOSE)
      fprintf(stderr,"\n Creating blob table and plotting\n");

    hap_plot(OUT,MAT,PAT,ASM,XDIM,YDIM,PDF,troot);

    for (i = 0; i < 2; i++)
      { if (ASM[i] == NULL)
          continue;
        sprintf(command,"Fastrm %s %s.%s %s.%s",ASM[i],ASM[i],MAT,ASM[i],PAT);
        system(command);
      }

    free(PAT);
    free(MAT);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
