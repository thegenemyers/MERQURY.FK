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

  //  Plotting constants

#define PLOTW  6.0    //  Plot width & height (in inches)
#define PLOTH  4.5

#define PLOTx  2.1    //  Peak relative x,y scale
#define PLOTy  1.1

#define PLOTX    0    //  0 => use relative x,y above
#define PLOTY    0

#define ZOPTION  1    //  0 if do not want -z

  //  Usage

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
//    OUT.completeness-stats
//    OUT.ASM_only.bed

static char template[16] = "._CNS.XXXX";

static int   VERBOSE;
static int   KMER;
static int   SOLID_THRESH;
static int64 SOLID_COUNT;

static int64 scan_asm(char *asmb, char *reads, char *out)
{ Profile_Index *AP, *RP;
  FILE   *qvs, *bed;
  uint16 *aprof, *rprof;
  int     pmax, plen;
  int     i, x, last;
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
  if (AP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK profile %s.prof\n",Prog_Name,asmb);
      exit (1);
    }
  if (RP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK relative profile %s.prof\n",Prog_Name,reads);
      exit (1);
    }
  if (AP->kmer != KMER)
    { fprintf(stderr,"\n%s: Profile %s.prof has wrong k-mer size %d != %d\n",
                     Prog_Name,asmb,AP->kmer,KMER);
      exit (1);
    }
  if (RP->kmer != KMER)
    { fprintf(stderr,"\n%s: Relative profile %s.prof has wrong k-mer size %d != %d\n",
                     Prog_Name,reads,RP->kmer,KMER);
      exit (1);
    }

  pmax  = 20000;
  aprof = Malloc(2*pmax*sizeof(uint16),"Profile array");
  rprof = aprof + pmax;

  bed = fopen(Catenate(out,".",asmb,"_only.bed"),"w");
  qvs = fopen(Catenate(out,".",asmb,".qv"),"w");

  fprintf(qvs,"Assembly Only\tTotal\tError %%\tQV\n");

  TOTS = 0;
  for (i = 0; i < AP->nreads; i++)
    { plen = Fetch_Profile(AP,i,pmax,aprof);
      if (plen > pmax)
        { pmax  = 1.2*plen + 1000;
          aprof = Realloc(aprof,2*pmax*sizeof(uint16),"Profile array");
          if (aprof == NULL)
            exit (1);
          rprof = aprof + pmax;
          Fetch_Profile(AP,i,pmax,aprof);
        }
      Fetch_Profile(RP,i,pmax,rprof);

      last = -1;
      miss = tots = 0;
      for (x = 0; x < plen; x++)
        { if (aprof[x] != 0)
            { if (rprof[x] == 0)
                { miss += 1;
                  if (x > last)
                    { if (last > 0)
                        fprintf(bed,"\t%d\n",last);
                      fprintf(bed,"%d\t%d",i,x);
                    }
                  last = x+KMER;
                }
              tots += 1;
            }
	}
      if (last > 0)
        fprintf(bed,"\t%d\n",last);

      err = 1. - pow(1.-(1.*miss)/tots,1./KMER);
      qv  = -10.*log10(err); 
      fprintf(qvs,"%lld\t%lld\t%.4f\t%.1f\n",miss,tots,err,qv);
      TOTS += tots;
    }
  fclose(bed);
  fclose(qvs);

  Free_Profiles(RP);
  Free_Profiles(AP);

  return (TOTS);
}

static void check_table(char *root)
{ FILE *f;
  int   kmer;

  f = fopen(Catenate(root,".ktab","",""),"r");
  if (f == NULL)
    { fprintf(stderr,"\n%s: Cannot find FastK table %s.ktab\n",Prog_Name,root);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (kmer != KMER)
        { fprintf(stderr,"\n%s: Kmer (%d) of table %s.ktab != %d\n",Prog_Name,kmer,root,KMER);
          exit (1);
        }
      fclose(f);
    }
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
          case 'p':
            if (strcmp("df",argv[i]+2) == 0)
              PDF = 1;
            else
              { fprintf(stderr,"%s: don't recognize option %s\n",Prog_Name,argv[i]);
                exit (1);
              }
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
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

    if (LINE+FILL+STACK == 0)
      LINE = FILL = STACK = 1;
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
    if (Rhist == NULL)
      { fprintf(stderr,"\n%s: Cannot find FastK histograam %s.hist\n",Prog_Name,READS);
        exit (1);
      }

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

    check_table(READS);

    for (i = 0; i < 2; i++)
      { char  *A = ASM[i];
        int64  miss, total;

        if (A == NULL)
          continue;

        check_table(A);

        //  Make a CN-spectra plot

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plots for %s\n",A);

        sprintf(Out,"%s.%s.spectra-cn",OUT,A);
        miss = cnplot(Out,A,READS,
                      PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                      LINE,FILL,STACK,troot,NTHREADS);

        //  Compute scaffold QV's and make bed file

        sprintf(ArelR,"%s.%s",A,READS);
        total = scan_asm(A,ArelR,OUT);

        //  Output global qv

        { FILE *qvs;

          double err, qv;
          if (i == 0)
            qvs = fopen(Catenate(OUT,"","",".qv"),"w");
          else
            qvs = fopen(Catenate(OUT,"","",".qv"),"a");
          err = 1. - pow(1.-(1.*miss)/total,1./KMER);
          qv  = -10.*log10(err); 
          if (i == 0)
            fprintf(qvs,"Assembly\tNo Supprt\tTotal\tError %%\tQV\n");
          fprintf(qvs,"%s\t%lld\t%lld\t%.4f\t%.1f\n",A,miss,total,100.*err,qv);
          fclose(qvs);
        }
      }

    //  1 diploid assembly ...

    if (ASM[1] == NULL)
      {
        //  Compute & output completeness stat

        if (VERBOSE)
          fprintf(stderr,"\n Computing completeness stats for %s\n",ASM[0]);

        { Histogram *H;
          FILE *cps;

          sprintf(command,"Logex -H1 -T%d '%s.0 = A-B[%d-]' %s %s",
                          NTHREADS,troot,SOLID_THRESH,ASM[0],READS);
          system(command);
      
          cps = fopen(Catenate(OUT,"","",".completeness_stats"),"a");

          sprintf(Hname,"%s.0",troot);
          H = Load_Histogram(Hname);
          fprintf(cps,"%s\tall\t%lld\t%lld\t%.2f\n",ASM[0],(SOLID_COUNT-H->hist[1]),SOLID_COUNT,
                                                    (SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
          Free_Histogram(H);

          fclose(cps);

          sprintf(command,"Fastrm %s.*.hist",troot);
          system(command);
        }

        //   Produce assembly-spectra plots

        if (VERBOSE)
          fprintf(stderr,"\n Making Assembly-spectra plot for %s\n",*ASM);

        sprintf(Out,"%s.spectra-asm",OUT);

        asmplot(Out,ASM[0],NULL,READS,
                PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                LINE,FILL,STACK,troot,NTHREADS);
      }

    //  2 haploid assemblies ...

    else
      { int64      miss, total;
        Histogram *H;
        FILE      *cps, *qvs;
        double     err, qv;

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plot for %s U %s\n",ASM[0],ASM[1]);

        //  Form the "k-mer union" of the two assemblies

        sprintf(command,"Logex -T%d -h1 '%s.U = A|+B' %s %s",NTHREADS,troot,ASM[0],ASM[1]);
        system(command);

        //  Make a CN spectra plot for the union

        sprintf(Out,"%s.spectra-cn",OUT);
        sprintf(A1uA2,"%s.U",troot);

        H = Load_Histogram(A1uA2);
        total = H->hist[1];
        Free_Histogram(H);

        miss = cnplot(Out,A1uA2,READS,
                      PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                      LINE,FILL,STACK,troot,NTHREADS);

        //  Output global qv

        qvs = fopen(Catenate(OUT,"","",".qv"),"a");
        err = 1. - pow(1.-(1.*miss)/total,1./KMER);
        qv  = -10.*log10(err); 
        fprintf(qvs,"Both\t%lld\t%lld\t%.4f\t%.1f\n",miss,total,100.*err,qv);
        fclose(qvs);

        //  Compute & output completeness stats

        if (VERBOSE)
          fprintf(stderr,"\n Computing completeness stats for %s and %s\n",ASM[0],ASM[1]);

        sprintf(command,
                "Logex -H1 -T%d '%s.0 = A-D[%d-]' '%s.1=B-D[%d-]' '%s.2=C-D[%d-]' %s %s %s %s",
                NTHREADS,troot,SOLID_THRESH,troot,SOLID_THRESH,troot,SOLID_THRESH,
                ASM[0],ASM[1],A1uA2,READS);
        system(command);
      
        cps = fopen(Catenate(OUT,"","",".completeness_stats"),"a");
        fprintf(cps,"Assembly\t%% Covered\n");

        sprintf(Hname,"%s.0",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"%s\tall\t%lld\t%lld\t%.2f\n",ASM[0],(SOLID_COUNT-H->hist[1]),SOLID_COUNT,
                                                  (SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
        Free_Histogram(H);

        sprintf(Hname,"%s.1",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"%s\tall\t%lld\t%lld\t%.2f\n",ASM[1],(SOLID_COUNT-H->hist[1]),SOLID_COUNT,
                                                  (SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
        Free_Histogram(H);

        sprintf(Hname,"%s.2",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"both\tall\t%lld\t%lld\t%.2f\n",(SOLID_COUNT-H->hist[1]),SOLID_COUNT,
                                                    (SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
        Free_Histogram(H);

        fclose(cps);

        //  Clean up

        sprintf(command,"Fastrm %s.*.hist %s.U",troot,troot);
        system(command);

        //   Produce assembly-spectra plots

        if (VERBOSE)
          fprintf(stderr,"\n Making Assembly-spectra plot for %s and %s\n",ASM[0],ASM[1]);

        sprintf(Out,"%s.spectra-asm",OUT);

        asmplot(Out,ASM[0],ASM[1],READS,
                PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                LINE,FILL,STACK,troot,NTHREADS);
      }
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
