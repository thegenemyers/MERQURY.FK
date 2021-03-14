/********************************************************************************************
 *
 *  Refactoring of Merqury CN-spectra script as a command line tool using FastK
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
#include "cn_plotter.h"
#include "asm_plotter.h"

static char *Usage = " [-v] [-T<int(4)>] [-pdf] [-lfs] <read> <asm1> [<asm2>] <out>";

//  Expected inputs from FastK ...
//    READS.hist     FastK -t1 -kKMER [...] <read_data> -NREAD
//    READS.ktab
//    ASM[i].ktab    FastK -t1 -p -kKMER [...] <assembly_i> -NASM[i]
//    ASM[i].prof
//    ASM[i].READS.prof   FastK -p:READS -kKMER [...] <assembly_i> -NASM[i].READS
//
//  Outputs:
//    OUT.ASM[i].spectra-cn.*
//    OUT.spectra-asm.*
//    OUT.spectra-cn.* (if 2 haploids)
//    OUT.qv
//    OUT.ASM[i].qv
//    OUT.completeness-stat
//    OUT.ASM_only.bed

static char template[15] = "._CN.XXXX";

static int   VERBOSE;
static int   KMER;
static int   SOLID_THRESH;
static int64 SOLID_COUNT;

static int64 scan_asm(char *asmb, char *reads, char *out)
{ Profile_Index *AP, *RP;
  FILE   *qvs, *bed;
  uint16 *aprof, *rprof;
  int     pmax, plen;
  int     i, x;
  int64   miss, tots;
  int64   TOTS;
  double  err, qv;

  //  In a scan of an assembly's profile and relative read profile:
  //    Compute per-scaffold qv, count the number of *non*-solid k-mers,
  //    and output a bed file of 0-read count intervals

  if (VERBOSE)
    fprintf(stderr,"\n Making .qv and .bed files for assembly %s\n",asmb);
  
  AP = Open_Profiles(asmb);
  RP = Open_Profiles(reads);

  pmax  = 20000;
  aprof = Malloc(2*pmax*sizeof(uint16),"Profile array");
  rprof = aprof + pmax*sizeof(uint16);

  bed = fopen(Catenate(asmb,"_only.bed","",""),"w");
  qvs = fopen(Catenate(out,".",asmb,".qv"),"w");

  fprintf(qvs,"Assembly Only\tTotal\tError %%\tQV\n");

  TOTS = 0;
  for (i = 0; i < AP->nreads; i++)
    { plen = Fetch_Profile(AP,i,pmax,aprof);
      if (plen > pmax)
        { pmax  = 1.2*plen + 1000;
          aprof = Realloc(aprof,2*pmax*sizeof(uint16),"Profile array");
          rprof = aprof + pmax*sizeof(uint16);
          Fetch_Profile(AP,i,pmax,aprof);
        }
      plen = Fetch_Profile(RP,0,pmax,rprof);

      miss = tots = 0;
      for (x = 0; x < plen; x++)
        { if (aprof[x] != 0)
            { if (rprof[x] == 0)
                { miss += 1;
                  fprintf(bed,"%d\t%d\t%d\n",i,x-KMER,x);
                }
              tots += 1;
            }
	}
      err = 1. - pow(1.-(1.*miss)/tots,1./KMER);
      qv  = -10.*log10(err); 
      fprintf(qvs,"%lld\t%lld\t%.2f\t%.1f",miss,tots,err,qv);
      TOTS += tots;
    }
  fclose(bed);
  fclose(qvs);

  Free_Profiles(RP);
  Free_Profiles(AP);

  return (TOTS);
}

/****************************************************************************************
 *
 *  Main Routine
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ char  *READS, *ASM[2], *OUT;
  int    LINE, FILL, STACK;
  int    PDF;
  int    NTHREADS;
  
  //  Command line processing

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) eptr;

    ARG_INIT("CNSpectra");

    PDF      = 0;
    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vlfs")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
          case 'p':
            if (strcmp("df",argv[i]+2) == 0)
              PDF = 1;
            else
              { fprintf(stderr,"%s: don't recognize option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    LINE    = flags['l'];
    FILL    = flags['f'];
    STACK   = flags['s'];

    if (argc < 4 || argc > 5)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
    else if (VERBOSE)
      { switch (argc)
        { case 4:
            fprintf(stderr,"\n Single diploid assembly\n");
            break;
          case 5:
            fprintf(stderr,"\n Two haploid assemblies\n");
            break;
          default:
            ;
        }
      }

    READS = argv[1];
    if (argc%2 == 1)
      { ASM[0] = argv[argc-3];
        ASM[1] = argv[argc-2];
      }
    else
      { ASM[0] = argv[argc-2];
        ASM[1] = NULL;
      }
    OUT = argv[argc-1];
  }

  //  Action ..

  { char       command[1000];
    char       Out[1000];
    char       ArelR[1000];
    char       A1uA2[1000];
    char       Hname[1000];
    char      *troot;
    Histogram *Rhist;
    int        i, k;

    troot = mktemp(template);

    //  Load read histogram and from it get KMER size and infer SOLID threshold and count

    Rhist = Load_Histogram(READS);
    Modify_Histogram(Rhist,Rhist->low,Rhist->high,1);

    KMER = Rhist->kmer;
    if (VERBOSE)
      fprintf(stderr,"\n Kmer size is %d\n",KMER);

    { int    low  = Rhist->low;
      int    high = Rhist->high;
      int64 *hist = Rhist->hist;

      for (k = low+1; hist[k] < hist[k-1]; k++) 
        if (k >= high)
          break;
      SOLID_THRESH = k;
      SOLID_COUNT  = 0;
      while (k <= high)
        SOLID_COUNT += hist[k++];

      if (VERBOSE)
        fprintf(stderr,"\n Solid k-mer cutoff is %d\n",SOLID_THRESH);
    }

    Free_Histogram(Rhist);

    //  For each assembly

    for (i = 0; i < 2; i++)
      { char  *A = ASM[i];
        int64  miss, total;

        if (A == NULL)
          continue;

        //  Make a CN-spectra plot

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plots for %s\n",A);

        sprintf(Out,"%s.%s.spectra-cn",OUT,A);
        miss = cnplot(Out,A,READS,6.0,4.5,2.1,1.1,0,0,PDF,1,LINE,FILL,STACK,troot,NTHREADS);

        //  Compute scaffold QV's and make bed file

        sprintf(ArelR,"%s.%s",A,READS);
        total = scan_asm(A,ArelR,OUT);

        //  Output global qv

        { FILE *qvs;

          double err, qv;
          qvs = fopen(Catenate(OUT,"","",".qv"),"a");
          err = 1. - pow(1.-(1.*miss)/total,1./KMER);
          qv  = -10.*log10(err); 
          if (i == 0)
            fprintf(qvs,"Assembly\tNo Supprt\tTotal\tError %%\tQV\n");
          fprintf(qvs,"%s\t%lld\t%lld\t%.2f\t%.1f",A,miss,total,100.*err,qv);
          fclose(qvs);
        }
      }

    //  1 diploid assembly ...

    if (ASM[1] == NULL)
      {
        //  Compute & output completeness stat

        { Histogram *H;
          FILE *cps;

          sprintf(command,"Logex -H1 -T%d '%s.0 = A-B[%d-]' %s %s",
                          NTHREADS,troot,SOLID_THRESH,ASM[0],READS);
      
          cps = fopen(Catenate(OUT,"","",".completeness_stat"),"w");
          fprintf(cps,"Assembly\t%% Covered\n");

          sprintf(Hname,"%s.0",troot);
          H = Load_Histogram(Hname);
          fprintf(cps,"%s\t%.2f\n",ASM[0],(SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
          Free_Histogram(H);

          fclose(cps);

          sprintf(command,"Fastrm %s.*.hist",troot);
          system(command);
        }

        //   Produce assembly-spectra plots

        if (VERBOSE)
          fprintf(stderr,"\n Making Assembly-spectra plot for %s\n",*ASM);

        asmplot(Out,ASM[0],NULL,READS,6.0,4.5,2.1,1.1,0,0,PDF,1,LINE,FILL,STACK,troot,NTHREADS);
      }

    //  2 haploid assemblies ...

    else
      { int64  miss, total;

        //  Form the "k-mer union" of the two assemblies

        sprintf(command,"Logex -T%d -h1 '%s.U = A|+B' %s %s",NTHREADS,troot,ASM[0],ASM[1]);
        system(command);

        //  Make a CN spectra plot for the union

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plot for %s U %s\n",ASM[0],ASM[1]);

        sprintf(Out,"%s.spectra-cn",OUT);
        sprintf(A1uA2,"%s.U",troot);
        miss = cnplot(Out,A1uA2,READS,6.0,4.5,2.1,1.1,0,0,PDF,1,LINE,FILL,STACK,troot,NTHREADS);

        total = 0;  // the histogram entry of troot.U

        //  Output global qv

        { FILE *qvs;

          double err, qv;
          qvs = fopen(Catenate(OUT,"","",".qv"),"a");
          err = 1. - pow(1.-(1.*miss)/total,1./KMER);
          qv  = -10.*log10(err); 
          fprintf(qvs,"Both\t%lld\t%lld\t%.2f\t%.1f",miss,total,100.*err,qv);
          fclose(qvs);
        }

        //  Compute & output completeness stats

        { Histogram *H;
          FILE      *cps;

          sprintf(command,
                  "Logex -H1 -T%d '%s.0 = A-D[%d-]' '%s.1=B-D[%d-]' '%s.2=C-D[%d-]' %s %s %s %s",
                  NTHREADS,troot,SOLID_THRESH,troot,SOLID_THRESH,troot,SOLID_THRESH,
                  ASM[0],ASM[1],A1uA2,READS);
      
          cps = fopen(Catenate(OUT,"","",".completeness_stat"),"w");
          fprintf(cps,"Assembly\t%% Covered\n");

          sprintf(Hname,"%s.0",troot);
          H = Load_Histogram(Hname);
          fprintf(cps,"%s\t%.2f\n",ASM[0],(SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
          Free_Histogram(H);

          sprintf(Hname,"%s.1",troot);
          H = Load_Histogram(Hname);
          fprintf(cps,"%s\t%.2f\n",ASM[1],(SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
          Free_Histogram(H);

          sprintf(Hname,"%s.2",troot);
          H = Load_Histogram(Hname);
          fprintf(cps,"Both\t%.2f\n",(SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
          Free_Histogram(H);

          fclose(cps);

          sprintf(command,"Fastrm %s.*.hist",troot);
          system(command);
        }

        sprintf(command,"Fastrm %s.U",troot);
        system(command);

        //   Produce assembly-spectra plots

        if (VERBOSE)
          fprintf(stderr,"\n Making Assembly-spectra plot for %s aand %s\n",ASM[0],ASM[1]);

        asmplot(Out,ASM[0],ASM[1],READS,6.0,4.5,2.1,1.1,0,0,PDF,1,LINE,FILL,STACK,troot,NTHREADS);
      }
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
