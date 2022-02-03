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
#include "hap_plotter.h"
#include "blk_plot.R.h"

  //  Plotting constants

#define PLOTW  6.0    //  Plot width & height (in inches)
#define PLOTH  4.5

#define PLOTx  2.1    //  Peak relative x,y scale
#define PLOTy  1.1

#define PLOTX    0    //  0 => use relative x,y above
#define PLOTY    0

#define ZOPTION  1    //  0 if do not want -z

  //  Phase Block parameters

static int ANCHOR_MARK    = 5;
static int ANCHOR_LENGTH  = 20000;

  //  Usage

static char *Usage[3] = { " [-v] [-lfs] [-pdf] [-T<int(4)>] [-P<dir(tmp)>]",
                          " <read>[.ktab] [ <mat>[.hap[.ktab]] <pat>[.hap[.ktab]] ]",
                          " <asm1:dna> [<asm2:dna>] <out>"
                        };

//  Expected inputs from FastK ...
//    READS.hist     FastK -t1 -kKMER [...] <read_data> -NREADS
//    READS.ktab
//    MAT.hap.ktab   HAPmaker MAT PAT READS
//    PAT.hap.ktab
//
//  Outputs:
//    OUT.ASM[i].spectra-cn.*
//    OUT.spectra-asm.*
//    OUT.spectra-cn.* (if 2 haploids)
//    OUT.qv
//    OUT.ASM[i].qv
//    OUT.ASM[i]_only.bed
//    OUT.completeness.stats
//  TRIO Outputs:
//    OUT.ASM[i].HAP.spectra-cn.*
//    OUT.hapmers.blob.*
//    OUT.ASM[i].phased_block.bed
//    OUT.ASM[i].phased_block.stats
//    OUT.ASM[i].phased_block.blob.*
//    OUT.ASM[i].block.N.[pdf|png]
//    OUT.ASM[i].continuity.N.[pdf|png]


static char template[16] = "._MQY.XXXX";

static int   VERBOSE;
static int   KMER;
static int   SOLID_THRESH;
static int64 SOLID_COUNT;

  //  In a scan of an assembly's profile and relative read profile:
  //    Compute per-scaffold qv in <out>.<asmb>.qv, and output a bed
  //    file of 0-read count intervals in <out>.<asmb>_only.bed
  //    Return the total # of bp's and the total # of 0-read counts.

static int64 *scan_asm(char *asmb, char *reads, char *out)
{ static int64 summary[2];

  Profile_Index *AP, *RP;
  FILE   *qvs, *bed;
  uint16 *aprof, *rprof;
  int     pmax, plen;
  int     i, x, last;
  int64   miss, tots;
  int64   TOTS, MISS;
  double  err, qv;

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
  MISS = 0;
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
      MISS += miss;
    }
  fclose(bed);
  fclose(qvs);

  Free_Profiles(RP);
  Free_Profiles(AP);

  summary[0] = MISS;
  summary[1] = TOTS;
  return (summary);
}


  //  Phase block codes: mark record

typedef struct
  { int beg;     //  interval [beg,end]
    int end;
    int mrk;     //  sum of majority marks  (+ = hap 1, - = hap 2)
    int opp;     //  sum of minority marks  (- = hap 1, + = hap 2)
  } Mark;

  //  mark[0..mtop) is an array of *pure* (opp=0) blocks, partition and merge blocks
  //    so that no block is unreliable save for the 1st and last, doing so in a way
  //    that minimized the impurity of the blocks.

static int merge_blocks(int mtop, Mark *mark)
{ int i, j, k;
  int last, first;

  //  First merge adjacent reliable blocks of the same polarity.

  last = -1;
  j = 0;
  for (i = 0; i < mtop; i++, j++)
    { mark[j] = mark[i];
      if (mark[i].end-mark[i].beg >= ANCHOR_LENGTH || abs(mark[i].mrk) >= ANCHOR_MARK)
        { if (last >= 0 && (mark[i].mrk < 0) == (mark[last].mrk < 0))
            { for (k = last+1; k < j; k += 2)
                { mark[last].opp += mark[k].mrk;
                  mark[last].mrk += mark[k+1].mrk;
                }
              mark[last].end  = mark[j].end;
              j = last;
            }
          else
            { if (last < 0)
                first = j;
              last = j;
            }
        }
    }
  mtop = j;

  //  Now decide how to divide any unreliable blocks between opposite polarity reliable
  //    blocks to arrive at a final partitioning

  if ((mtop-last) % 2 == 0)
    { mark[mtop].beg = mark[mtop].end = mark[mtop-1].end;
      mark[mtop].mrk = mark[mtop].opp = 0;
      mtop += 1;
    }
  mark[mtop].mrk = mark[mtop].opp = 0;
  mark[mtop].end = mark[mtop-1].end;
  mark[mtop].beg = mark[mtop].end - ANCHOR_LENGTH;

  if (first % 2 != 0)
    { mark[-1].beg = mark[-1].end = mark[0].beg;
      mark[-1].mrk = mark[-1].opp = 0;
      last = -2;
    }
  else
    last = -1;

  j = 0;
  for (i = 0; i <= mtop; i++, j++)
    { mark[j] = mark[i];
      if (mark[i].end-mark[i].beg >= ANCHOR_LENGTH || abs(mark[i].mrk) >= ANCHOR_MARK)
        { int score, best, bidx;

          score = 0;
          for (k = j-1; k > last; k -= 2)
            score += mark[k].mrk;

          best = abs(score);
          bidx = k;
          for (k += 2 ; k < j; k += 2)
            { score -= (mark[k-1].mrk + mark[k].mrk);
              if (abs(score) < best)
                { best = abs(score);
                  bidx = k;
                }
              else if (abs(score) == best && last >= 0 && abs(mark[last].mrk) > abs(mark[j].mrk))
                { best = abs(score);
                  bidx = k;
                }
            }

          if (last >= 0)
            { for (k = last+2; k <= bidx; k += 2) 
                { mark[last].mrk += mark[k].mrk;
                  mark[last].opp += mark[k-1].mrk;
                }
              mark[last++].end = mark[bidx].end;
            }
          else if (bidx >= 0)
            { int tmrk, topp;

              tmrk = topp = 0;
              for (k = last+2; k <= bidx; k += 2) 
                { tmrk += mark[k].mrk;
                  topp += mark[k-1].mrk;
                }
              mark[0].mrk = tmrk;
              mark[0].opp = topp;
              mark[0].end = mark[bidx].end;
              last = 1;
            }
          else
            last = 0;
          
          for (k = bidx+2; k < j; k += 2)
            { mark[j].mrk += mark[k-1].mrk;
              mark[j].opp += mark[k].mrk;
            }
          mark[j].beg = mark[bidx+1].beg;

          mark[last] = mark[j];
          j = last;
        }
    }

  if (mark[j-1].mrk == 0)
    j -= 1;

  return (j);
}

static int64 BTOT;   //  block length sum, #, min & max
static int   NBLK;   //  set by phase_blocks, consumed by block_stats
static int   BMIN;
static int   BMAX;

  //  In a scan of an assembly's profile and relative hapmer profiles:
  //    Compute a phased-block partitioning of the contigs of the assembly
  //    and output them to <out>.<asmb>.phased_block.bed.
  //    Output temporary files containing the sizes of blocks, contigs,
  //    and scaffolds in <troot>.[blk+ctg+scf].un (these are later sorted
  //    and used to produce N-curves.
  //    Output a listing of each block with its hapmer composition for the
  //    production of a blob plot of the blocks.

static int phase_blocks(char *asmb, char *mat, char *pat, char *out, char *troot)
{ Profile_Index *AP, *MP, *PP;
  FILE   *bed, *plot, *fscf, *fctg, *fblk;
  uint16 *aprof, *mprof, *pprof;
  int     pmax, plen;
  Mark   *mark;
  int     nctg, nscf;
  int64   stot, mtot;
  int     p, i, x;

  AP = Open_Profiles(Catenate(asmb,"","",""));
  MP = Open_Profiles(Catenate(asmb,".",mat,""));
  PP = Open_Profiles(Catenate(asmb,".",pat,""));
  if (AP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK profile %s.prof\n",Prog_Name,asmb);
      exit (1);
    }
  if (MP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK relative profile %s.%s.prof\n",
                     Prog_Name,asmb,mat);
      exit (1);
    }
  if (PP == NULL)
    { fprintf(stderr,"\n%s: Cannot open/find FastK relative profile %s.%s.prof\n",
                     Prog_Name,asmb,pat);
      exit (1);
    }

  pmax  = 20000;
  aprof = Malloc((3*pmax+2)*sizeof(uint16) + (pmax/KMER+4)*sizeof(Mark),"Profile array");
  if (aprof == NULL)
    exit (1);
  mprof = aprof + (pmax+2);
  pprof = mprof + pmax;
  mark  = ((Mark *) (pprof + pmax)) + 1;

  fscf = fopen(Catenate(troot,".scf.un","",""),"w");
  fctg = fopen(Catenate(troot,".ctg.un","",""),"w");
  fblk = fopen(Catenate(troot,".blk.un","",""),"w");

  bed  = fopen(Catenate(out,".",asmb,".phased_block.bed"),"w");
  fprintf(bed,"Scaffold\tStart\tEnd\tPhase\tPurity\tSwitches\tMarkers\n");

  plot = fopen(Catenate(troot,".hpi","",""),"w");
  fprintf(plot,"Block\tRange\t%s\t%s\tSize\n",mat,pat);

  stot = 0;
  mtot = 0;
  BTOT = 0;
  NBLK = 0;
  BMIN = 0x7fffffff;
  BMAX = 0;
  nctg = 0;
  nscf = AP->nreads;
  for (p = 0; p < nscf; p++)
    { int beg, end, frst;
      int mrk, sum;
      int eoc, mtop;

      plen = Fetch_Profile(AP,p,pmax,aprof);
      if (plen > pmax)
        { pmax  = 1.2*plen + 1000;
          aprof = Realloc(aprof,(3*pmax+2)*sizeof(uint16) + (pmax/KMER+4)*sizeof(Mark),
                                "Profile array");
          if (aprof == NULL)
            exit (1);
          mprof = aprof + (pmax+2);
          pprof = mprof + pmax;
          mark  = ((Mark *) (pprof + pmax)) + 1;
          Fetch_Profile(AP,p,pmax,aprof);
        }
      Fetch_Profile(MP,p,pmax,mprof);
      Fetch_Profile(PP,p,pmax,pprof);
      aprof[plen] = 0;
      aprof[plen+1] = 1;
      fprintf(fscf,"scaffold\t%s\t%d\n",asmb,plen);

      //  For each contig, build stack of hap sites (= overlapping hap-mers), merging adjacent
      //    sites of the same polarity and then call process_blocks to finish the partitioning
      //    into putative phased blocks

      eoc  = 0;
      end  = 0;
      frst = 0;
      mtop = 0;
      for (x = 0; x <= plen; x++)
        { int d;

          if (aprof[x] == 0)
            eoc = 1;
          else
            { if (mprof[x] > 0)
                d = 1;
              else if (pprof[x] > 0)
                d = -1;
              else
                continue;
              if (x < end)
                { mrk += d;
                  end  = x+KMER;
                  sum += 1;
                  continue;
                }
            }
          if (end > 0 && abs(mrk) >= (sum+3)/4)
            { if (mtop > 0 && (mrk < 0) == (mark[mtop-1].mrk < 0))
                { mark[mtop-1].end = end;
                  if (mrk < 0)
                    mark[mtop-1].mrk += (sum+mrk)/2 - (end-beg);
                  else
                    mark[mtop-1].mrk += (end-beg) - (sum-mrk)/2;
                }
              else
                { mark[mtop].beg = beg;
                  mark[mtop].end = end;
                  mark[mtop].opp = 0;
                  if (mrk < 0)
                    mark[mtop].mrk = (sum+mrk)/2 - (end-beg);
                  else
                    mark[mtop].mrk = (end-beg) - (sum-mrk)/2;
                  mtop += 1;
                }
            }
          if (eoc)
            { int ns, to, len;

              if (mtop > 0)
                { nctg += 1;
                  fprintf(fctg,"contig\t%s\t%d\n",asmb,(x-frst)+(KMER-1));

                  mtop = merge_blocks(mtop,mark);

                  mark[0].beg      = frst;
                  mark[mtop-1].end = x + (KMER-1);
                  for (i = 0; i < mtop; i++)
                    { len = mark[i].end - mark[i].beg;
                      if (mark[i].mrk > 0)
                        { ns = -mark[i].opp;
                          to = ns + mark[i].mrk;
                          fprintf(bed,"%d\t%d\t%d\t%s\t%.3f\t%d\t%d\n",
                                      nctg,mark[i].beg,mark[i].end,mat,(100.*ns)/to,ns,to);
                          fprintf(plot,"%s\t%d\t%d\t%d\t%d\n",mat,NBLK+i,to-ns,ns,len);
                          fprintf(fblk,"block\t%s\t%d\n",mat,len);
                        }
                      else
                        { ns = mark[i].opp;
                          to = ns - mark[i].mrk;
                          fprintf(bed,"%d\t%d\t%d\t%s\t%.3f\t%d\t%d\n",
                                      nctg,mark[i].beg,mark[i].end,pat,(100.*ns)/to,ns,to);
                          fprintf(plot,"%s\t%d\t%d\t%d\t%d\n",pat,NBLK+i,ns,to-ns,len);
                          fprintf(fblk,"block\t%s\t%d\n",pat,len);
                        }
                      stot += ns;
                      mtot += to;
                      BTOT += len;
                      if (len > BMAX)
                        BMAX = len;
                      else if (len < BMIN)
                        BMIN = len;
                    }
                  NBLK += mtop;
                }
               
              while (aprof[x+1] == 0)
                x += 1;
              eoc  = 0;
              end  = 0;
              mtop = 0;
              frst = x+1;
            }
          else
            { beg = x;
              end = x+KMER;
              mrk = d;
              sum = 1;
            }
        }
    }
  fprintf(bed,"total\t-t-\t-\t%.3f\t%lld\t%lld\n",
              (100.*stot)/mtot,stot,mtot);

  fclose(plot);
  fclose(bed);
  fclose(fblk);
  fclose(fctg);
  fclose(fscf);

  Free_Profiles(MP);
  Free_Profiles(PP);
  Free_Profiles(AP);

  return (nctg != nscf);
}

static void block_stats(char *asmb, char *out, char *troot)
{ FILE *sfile, *nfile;
  int   bn50;
  int64 sum, thr;

  nfile = fopen(Catenate(troot,".block.sizes","",""),"r");
  sum = 0;
  thr = BTOT/2;
  while (fscanf(nfile,"block\t%*s\t%d\n",&bn50) == 1)
    { sum += bn50;
      if (sum >= thr)
        break;
    }
  fclose(nfile);
  
  sfile = fopen(Catenate(out,".",asmb,".phased_block.stats"),"w");
  fprintf(sfile,"# Blocks\tSum\tMin.\tAvg.\tN50\tMax.\n");
  fprintf(sfile,"%d\t%lld\t%d\t%lld\t%d\t%d\n",NBLK,BTOT,BMIN,BTOT/NBLK,bn50,BMAX);
  fclose(sfile);
}


/****************************************************************************************
 *
 *  Main Routine
 *
 *****************************************************************************************/

static int check_table(char *name, int lmer)
{ int   kmer;
  FILE *f;

  f = fopen(name,"r");
  if (f == NULL)
    { fprintf(stderr,"\n%s: Cannot find FastK table %s\n",Prog_Name,name);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (lmer != 0 && kmer != lmer)
        { fprintf(stderr,"\n%s: Kmer (%d) of table %s != %d\n",Prog_Name,kmer,name,lmer);
          exit (1);
        }
      fclose(f);
      return (kmer);
    }
}

int main(int argc, char *argv[])
{ char  *READS, *ASM[2], *MAT, *PAT, *OUT;
  int    LINE, FILL, STACK;
  int    PDF;
  int    NTHREADS;
  char  *SORT_PATH;

  char   command[1000];
  char  *troot;
  
  //  Command line processing

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) eptr;

    ARG_INIT("MerquryFK");

    PDF      = 0;
    NTHREADS = 4;
    SORT_PATH = "/tmp";

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
    LINE    = flags['l'];
    FILL    = flags['f'];
    STACK   = flags['s'];

    if (argc < 4 || argc > 7)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose output to stderr\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        fprintf(stderr,"      -P: Place all temporary files in directory -P.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: draw line plot\n");
        fprintf(stderr,"      -f: draw fill plot\n");
        fprintf(stderr,"      -s: draw stack plot\n");
        fprintf(stderr,"          any combo allowed, none => draw all\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        exit (1);
      }

    if (LINE+FILL+STACK == 0)
      LINE = FILL = STACK = 1;

    READS = argv[1];
    switch (argc)
    { case 4:
        if (VERBOSE)
          fprintf(stderr,"\n Single diploid assembly, no trio data\n");
        MAT    = NULL;
        ASM[0] = argv[2];
        ASM[1] = NULL;
        break;
      case 5:
        if (VERBOSE)
          fprintf(stderr,"\n Two haploid assemblies, no trio data\n");
        MAT    = NULL;
        ASM[0] = argv[2];
        ASM[1] = argv[3];
        break;
      case 6:
        if (VERBOSE)
          fprintf(stderr,"\n Single diploid assembly with trio data\n");
        MAT    = argv[2];
        PAT    = argv[3];
        ASM[0] = argv[4];
        ASM[1] = NULL;
        break;
      case 7:
        if (VERBOSE)
          fprintf(stderr,"\n Two haploid assemblies with trio data\n");
        MAT    = argv[2];
        PAT    = argv[3];
        ASM[0] = argv[4];
        ASM[1] = argv[5];
        break;
      default:
        ;
    }
    OUT = argv[argc-1];

    //  Remove any suffixes from argument names

    { char *suffix[9] = { ".gz", ".fa", ".fq", ".fasta", ".fastq", ".db", ".sam", ".bam", ".cram" };

      troot = mktemp(template);

      for (i = 0; i < 2; i++)
        { char *A = ASM[i];
          int   len;

          if (A == NULL)
            continue;

          for (j = 0; j < 9; j++)
            { len = strlen(A) - strlen(suffix[j]);
              if (strcmp(A+len,suffix[j]) == 0)
                A[len] = '\0';
            }
        }

      READS = Root(READS,".ktab");
      if (MAT != NULL)
        { char *x;

          x   = Root(MAT,".ktab");
          MAT = Root(x,".hap");
          free(x);
          x   = Root(PAT,".ktab");
          PAT = Root(x,".hap");
          free(x);
        }

      KMER = check_table(Catenate(READS,".ktab","",""),0);
      if (MAT != NULL)
        { KMER = check_table(Catenate(MAT,".hap",".ktab",""),KMER);
          KMER = check_table(Catenate(PAT,".hap",".ktab",""),KMER);
        }

      if (VERBOSE)
        fprintf(stderr,"\n Kmer size is %d\n",KMER);

      ANCHOR_MARK *= KMER;
    }
  }

  //  Non-Trio action ...

  { char       Out[1000];
    char       ArelR[1000];
    char       A1uA2[1000];
    char       Hname[1000];
    Histogram *Rhist;
    int        i, k;

    //  Load read histogram and infer SOLID threshold and count

    Rhist = Load_Histogram(READS);
    if (Rhist == NULL)
      { fprintf(stderr,"\n%s: Cannot find FastK histograam %s.hist\n",Prog_Name,READS);
        exit (1);
      }
    if (KMER != Rhist->kmer)
      { fprintf(stderr,"\n%s: Kmer size of %s.hist doesn't match %s.ktab\n",Prog_Name,READS,READS);
        exit (1);
      }

    Modify_Histogram(Rhist,Rhist->low,Rhist->high,1);

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

    { FILE  *qvs;
      int64 *summary, totl, miss;

      qvs = fopen(Catenate(OUT,"","",".qv"),"w");
      fprintf(qvs,"Assembly\tNo Support\tTotal\tError %%\tQV\n");

      for (i = 0; i < 2; i++)
        { char  *A = ASM[i];
          double err, qv;

          if (A == NULL)
            continue;

          //  Create assembly tables and profiles

          sprintf(command,"FastK -k%d -T%d -P%s -t1 -p %s",KMER,NTHREADS,SORT_PATH,A); 
          system(command);

          sprintf(command,"FastK -k%d -T%d -P%s -p:%s %s -N%s.%s",
                          KMER,NTHREADS,SORT_PATH,READS,A,A,READS); 
          system(command);

          //  Make a CN-spectra plot

          if (VERBOSE)
            fprintf(stderr,"\n Making CN-spectra plots for %s\n",A);

          sprintf(Out,"%s.%s.spectra-cn",OUT,A);
          cn_plot(Out,A,READS,
                  PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                  LINE,FILL,STACK,troot,NTHREADS);

          //  Compute scaffold QV's and make bed file

          sprintf(ArelR,"%s.%s",A,READS);
          summary = scan_asm(A,ArelR,OUT);

          //  Output global qv

          err = 1. - pow(1.-(1.*summary[0])/summary[1],1./KMER);
          qv  = -10.*log10(err); 
          fprintf(qvs,"%s\t%lld\t%lld\t%.4f\t%.1f\n",A,summary[0],summary[1],100.*err,qv);

          if (i == 0)
            { miss = summary[0];
              totl = summary[1];
            }
          else
            { miss += summary[0];
              totl += summary[1];
              err = 1. - pow(1.-(1.*miss)/totl,1./KMER);
              qv  = -10.*log10(err); 
              fprintf(qvs,"both\t%lld\t%lld\t%.4f\t%.1f\n",miss,totl,100.*err,qv);
            }
        }

      fclose(qvs);
    }

    //  1 diploid assembly ...

    if (ASM[1] == NULL)
      { Histogram *H;
        FILE      *cps;
       
        //  Compute & output completeness stat

        if (VERBOSE)
          fprintf(stderr,"\n Computing completeness stats for %s\n",ASM[0]);

        sprintf(command,"Logex -H1 -T%d '%s.0 = A-B[%d-]' %s %s",
                        NTHREADS,troot,SOLID_THRESH,ASM[0],READS);
        system(command);
      
        cps = fopen(Catenate(OUT,"","",".completeness.stats"),"a");
        fprintf(cps,"Assembly\tRegion\tFound\tTotal\t%% Covered\n");

        sprintf(Hname,"%s.0",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"%s\tall\t%lld\t%lld\t%.2f\n",ASM[0],(SOLID_COUNT-H->hist[1]),SOLID_COUNT,
                                                  (SOLID_COUNT-H->hist[1])/(.01*SOLID_COUNT));
        Free_Histogram(H);

        fclose(cps);

        sprintf(command,"Fastrm %s.*.hist",troot);
        system(command);

        //   Produce assembly-spectra plots

        if (VERBOSE)
          fprintf(stderr,"\n Making Assembly-spectra plot for %s\n",*ASM);

        sprintf(Out,"%s.spectra-asm",OUT);

        asm_plot(Out,ASM[0],NULL,READS,
                 PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                 LINE,FILL,STACK,troot,NTHREADS);
      }

    //  2 haploid assemblies ...

    else
      { Histogram *H;
        FILE      *cps;

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plot for %s U %s\n",ASM[0],ASM[1]);

        //  Form the "k-mer union" of the two assemblies

        sprintf(command,"Logex -T%d '%s.U = A|+B' %s %s",NTHREADS,troot,ASM[0],ASM[1]);
        system(command);

        //  Make a CN spectra plot for the union

        sprintf(Out,"%s.spectra-cn",OUT);
        sprintf(A1uA2,"%s.U",troot);

        cn_plot(Out,A1uA2,READS,
                PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                LINE,FILL,STACK,troot,NTHREADS);

        //  Compute & output completeness stats

        if (VERBOSE)
          fprintf(stderr,"\n Computing completeness stats for %s and %s\n",ASM[0],ASM[1]);

        sprintf(command,
                "Logex -H1 -T%d '%s.0 = A&.D[%d-]' '%s.1=B&.D[%d-]' '%s.2=C&.D[%d-]' %s %s %s %s",
                NTHREADS,troot,SOLID_THRESH,troot,SOLID_THRESH,troot,SOLID_THRESH,
                ASM[0],ASM[1],A1uA2,READS);
        system(command);
      
        cps = fopen(Catenate(OUT,"","",".completeness.stats"),"a");
        fprintf(cps,"Assembly\tRegion\tFound\tTotal\t%% Covered\n");

        sprintf(Hname,"%s.0",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"%s\tall\t%lld\t%lld\t%.2f\n",ASM[0],H->hist[1],SOLID_COUNT,
                                                  H->hist[1]/(.01*SOLID_COUNT));
        Free_Histogram(H);

        sprintf(Hname,"%s.1",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"%s\tall\t%lld\t%lld\t%.2f\n",ASM[1],H->hist[1],SOLID_COUNT,
                                                  H->hist[1]/(.01*SOLID_COUNT));
        Free_Histogram(H);

        sprintf(Hname,"%s.2",troot);
        H = Load_Histogram(Hname);
        fprintf(cps,"both\tall\t%lld\t%lld\t%.2f\n",H->hist[1],SOLID_COUNT,
                                                    H->hist[1]/(.01*SOLID_COUNT));
        Free_Histogram(H);

        fclose(cps);

        //  Clean up

        sprintf(command,"Fastrm %s.*.hist %s.U",troot,troot);
        system(command);

        //   Produce assembly-spectra plots

        if (VERBOSE)
          fprintf(stderr,"\n Making Assembly-spectra plot for %s and %s\n",ASM[0],ASM[1]);

        sprintf(Out,"%s.spectra-asm",OUT);

        asm_plot(Out,ASM[0],ASM[1],READS,
                 PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                 LINE,FILL,STACK,troot,NTHREADS);
      } 
  }

  if (MAT == NULL)
    goto clean_up;

  //  Trio actions ...

  { char       Out[1000];
    char       Hname[1000];
    FILE      *cps;
    int        i;

    //  CN spectra of the hapmers versus the assembly(ies)

    cps = fopen(Catenate(OUT,"","",".completeness.stats"),"a");

    for (i = 0; i < 2; i++)
      { char      *A = ASM[i];
        Histogram *H, *G;

        if (A == NULL)
          continue;

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plots for %s versus assembly %s\n",MAT,A);

        sprintf(Out,"%s.%s.%s.spectra-cn",OUT,A,MAT);
        cn_plot(Out,A,MAT,
                PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                LINE,FILL,STACK,troot,NTHREADS);

        if (VERBOSE)
          fprintf(stderr,"\n Making CN-spectra plots for %s versus assembly %s\n",PAT,A);

        sprintf(Out,"%s.%s.%s.spectra-cn",OUT,A,PAT);
        cn_plot(Out,A,PAT,
                PLOTW,PLOTH,PLOTx,PLOTy,PLOTX,PLOTY,PDF,ZOPTION,
                LINE,FILL,STACK,troot,NTHREADS);

        if (VERBOSE)
          fprintf(stderr,"\n Determining hap-mer completeness of assembly %s\n",A);

        sprintf(command,
                "Logex -H1 -T%d '%s.M = A&.B' '%s.P=A&.C' '%s.0=B' '%s.1=C' %s %s %s",
                NTHREADS,troot,troot,troot,troot,A,MAT,PAT);
        system(command);

        sprintf(Hname,"%s.M",troot);   //  unique or instances?
        H = Load_Histogram(Hname);
        sprintf(Hname,"%s.0",troot);
        G = Load_Histogram(Hname);

        fprintf(cps,"%s\t%s\t%lld\t%lld\t%.2f\n",A,MAT,H->hist[1],G->hist[1],
                                                       H->hist[1]/(.01*G->hist[1]));

        Free_Histogram(G);
        Free_Histogram(H);

        sprintf(Hname,"%s.P",troot);   //  unique or instances?
        H = Load_Histogram(Hname);
        sprintf(Hname,"%s.1",troot);
        G = Load_Histogram(Hname);

        fprintf(cps,"%s\t%s\t%lld\t%lld\t%.2f\n",A,PAT,H->hist[1],G->hist[1],
                                                       H->hist[1]/(.01*G->hist[1]));

        Free_Histogram(G);
        Free_Histogram(H);

        //  Clean up

        sprintf(command,"Fastrm %s.*.hist",troot);
        system(command);
      }

    fclose(cps);

    //  Generate the R plot script for block.N and continuity plots

    cps = fopen(Catenate(troot,".N.R","",""),"w");
    fwrite(blk_plot_script,strlen(blk_plot_script),1,cps);
    fclose(cps);

    //  For each assembly, determine the block phasing and blob plot data file

    for (i = 0; i < 2; i++)
      { char *A = ASM[i];
        int   do_scaffs;

        if (A == NULL)
          continue;

        if (VERBOSE)
          fprintf(stderr,"\n Producing relative profiles for phasing block calculation on %s\n",A);

        sprintf(command,"FastK -T%d -P%s -k%d -p:%s.hap %s -N%s.%s",
                        NTHREADS,SORT_PATH,KMER,MAT,A,A,MAT);
        system(command);

        sprintf(command,"FastK -T%d -P%s -k%d -p:%s.hap %s -N%s.%s",
                        NTHREADS,SORT_PATH,KMER,PAT,A,A,PAT);
        system(command);

        if (VERBOSE)
          fprintf(stderr,"\n Computing phasing blocks for assemlby %s\n",A);

        do_scaffs = phase_blocks(A, MAT, PAT, OUT, troot);

        if (VERBOSE)
          fprintf(stderr,"\n Producing phased blob plot of the imputed blocks of assembly %s\n",A);

        sprintf(Out,"%s.%s.phased_block.blob",OUT,A);
        hap_plot(Out,MAT,PAT,ASM,PLOTW,PLOTH,PDF,troot);

        if (VERBOSE)
          fprintf(stderr,"\n Sorting phase block and assembly sizes for %s\n",A);

        sprintf(command,"sort -nr -k3 %s.scf.un >%s.scaff.sizes",troot,troot);
        system(command);

        sprintf(command,"sort -nr -k3 %s.ctg.un >%s.contig.sizes",troot,troot);
        system(command);

        sprintf(command,"sort -nr -k3 %s.blk.un >%s.block.sizes",troot,troot);
        system(command);

        //  Output block partition stats now that have block.sizes needed to compute N50

        block_stats(A,OUT,troot);

        //  Plot the block.N and continuity N-curves

        sprintf(command,"Rscript %s.N.R -b %s.block.sizes -c %s.contig.sizes",troot,troot,troot);
        if (do_scaffs)
          sprintf(command+strlen(command)," -s %s.scaff.sizes",troot);
        sprintf(command+strlen(command)," -o %s.%s%s -x %g -y %g 2>/dev/null",
                                        OUT,A,PDF?" -p":" ",PLOTW,PLOTH);
        system(command);

        sprintf(command,"rm %s.scf.un %s.scaff.sizes",troot,troot);
        system(command);

        sprintf(command,"rm  %s.ctg.un %s.contig.sizes",troot,troot);
        system(command);

        sprintf(command,"rm %s.blk.un %s.block.sizes",troot,troot);
        system(command);
      }

    sprintf(command,"rm -f %s.N.R",troot);
    system(command);

    //  Plot phased contig blobs

    if (VERBOSE)
      { if (ASM[1] == NULL)
          fprintf(stderr,"\n Producing phased blob plots of the contigs of the assembly\n");
        else
          fprintf(stderr,"\n Producing phased plots of the contigs of the assemblies\n");
      }

    sprintf(Out,"%s.hapmers.blob",OUT);
    hap_plot(Out,MAT,PAT,ASM,PLOTW,PLOTH,PDF,troot);
  }

clean_up:
  { int  i;

    for (i = 0; i < 2; i++)
      { char *A = ASM[i];
 
        if (A == NULL)
          continue;

        if (MAT != NULL)
          sprintf(command,"Fastrm %s %s.%s %s.%s %s.%s",A,A,READS,A,MAT,A,PAT);
        else
          sprintf(command,"Fastrm %s %s.%s",A,A,READS);
        system(command);
      }

    if (MAT != NULL)
      { free(MAT);
        free(PAT);
      }
    free(READS);
  }

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
