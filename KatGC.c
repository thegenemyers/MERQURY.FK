/*******************************************************************************************
 *
 *  KatGC
 *
 *  Author:  Gene Myers
 *  Date  :  May, 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <pthread.h>

#undef   DEBUG
#undef   DEBUG_THREADS

#include "libfastk.h"

#include "kgc_plot.R.h"

static char *Usage[] = { " [-w<double(6.0)>] [-h<double(4.5)>]",
                         " [-[xX]<number(x2.1)>] [-lfs] [-pdf] [-T<int(4)>]",
                         " <source>[.ktab] <out>"
                       };

static int NTHREADS;
static int KMER;
static int HMAX;

/****************************************************************************************
 *
 *  Streaming eval
 *
 *****************************************************************************************/

typedef struct
  { int          tid;
    Kmer_Stream *T;
    int64        beg;
    int64        end;
    int64      **plot;
  } TP;

static int    GC[256];
static int    GCR[256];

static void gc_setup(int kmer)
{ static int isgc[4] = { 0, 1, 1, 0 };
  uint32 x;

  for (x = 0; x < 256; x++)
    { GC[x] = isgc[(x>>6)&0x3] + isgc[(x>>4)&0x3] + isgc[(x>>2)&0x3] + isgc[x&0x3];
      switch( kmer % 4)
      { case 0:
          GCR[x] = GC[x];
          break;
        case 1:
          GCR[x] = isgc[(x>>6)&0x3];
          break;
        case 2:
          GCR[x] = isgc[(x>>6)&0x3] + isgc[(x>>4)&0x3];
          break;
        case 3:
          GCR[x] = isgc[(x>>6)&0x3] + isgc[(x>>4)&0x3] + isgc[(x>>2)&0x3];
          break;
      }
    }
}

static inline int gcontent(uint8 *a, int kbyte)
{ int i, cnt;

  cnt = 0;
  for (i = 1; i < kbyte; i++)
    cnt += GC[*a++];
  return (cnt+GCR[*a]);
}

static void *merge_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  Kmer_Stream  *T     = parm->T;
  int64         beg   = parm->beg;
  int64         end   = parm->end;

  int kbyte = T->kbyte;

  int64 **plot = NULL;
  uint8  *entry;
  int     i;

#ifdef DEBUG_TRACE
  char *buffer;
#endif

  plot    = Malloc(sizeof(int64 *)*(KMER+1),"Allocating thread working memory");
  plot[0] = Malloc(sizeof(int64)*(HMAX+1)*(KMER+1),"Allocating plot");
  for (i = 1; i <= KMER; i++)
    plot[i] = plot[i-1] + (HMAX+1);
  bzero(plot[0],sizeof(int64)*(HMAX+1)*(KMER+1));

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  printf(" [%lld-%lld]",beg,end);
  printf("\n");
#endif

  if (tid != 0)
    T = Clone_Kmer_Stream(T);

  entry = Current_Entry(T,NULL);

  GoTo_Kmer_Index(T,beg);

  while (T->cidx < end)
    { int tn, un;

      tn = Current_Count(T);
      un = gcontent(Current_Entry(T,entry),kbyte);
      if (tn <= HMAX)
        plot[un][tn] += 1;
      Next_Kmer_Entry(T);
    }

  if (tid != 0)
    Free_Kmer_Stream(T);

  parm->plot = plot;
  return (NULL);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

static char template[15] = "._KGC.XXXX";

static int check_table(char *name, int lmer)
{ int   kmer;
  FILE *f;

  f = fopen(name,"r");
  if (f == NULL)
    { fprintf(stderr,"%s: Cannot find FastK file %s\n",Prog_Name,name);
      exit (1);
    }
  else
    { fread(&kmer,sizeof(int),1,f);
      if (lmer != 0 && kmer != lmer)
        { fprintf(stderr,"%s: Kmer (%d) of %s != %d\n",Prog_Name,kmer,name,lmer);
          exit (1);
        }
      fclose(f);
    }
  return (kmer);
}

int main(int argc, char *argv[])
{ int     LINE, FILL, BOTH;
  int     PDF;
  double  XDIM, YDIM;
  double  XREL;
  int     XMAX;
  char   *OUT;
  char   *SOURCE;

  int64 **PLOT;
  int64   ZMAX;

  //  Process command line
 
  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("KatGC");

    XDIM = 6.0;
    YDIM = 4.5;
    XREL = 2.1;
    XMAX = 0;
    PDF  = 0;
    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("lfs")
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
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
          case 'X':
            ARG_POSITIVE(XMAX,"x max");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    LINE = flags['l'];
    FILL = flags['f'];
    BOTH = flags['s'];

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -w: width in inches of plots\n");
        fprintf(stderr,"      -h: height in inches of plots\n");
        fprintf(stderr,"      -x: max x as a real-valued multiple of x* with max\n");
        fprintf(stderr,"              count 'peak' away from the origin\n");
        fprintf(stderr,"      -X: max x as an int value in absolute terms\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: draw a contour map\n");
        fprintf(stderr,"      -f: draw a heat map\n");
        fprintf(stderr,"      -s: draw a heat map with contour overlay\n");
        fprintf(stderr,"          any combo allowed, none => draw all\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        exit (1);
      }

    if (LINE+FILL+BOTH == 0)
      LINE = FILL = BOTH = 1;
    SOURCE = Root(argv[1],".ktab");
    OUT    = argv[2];

    KMER = check_table(Catenate(SOURCE,".ktab","",""),0);

    if (XMAX == 0)
      HMAX = 1000;
    else
      HMAX = XMAX;
  }

  //  Compute the GC vs K-mer matrix

  { Kmer_Stream  *T;
    int64     range[NTHREADS+1];
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif
    TP        parm[NTHREADS];
    char     *seq;
    uint8    *ent;
    int       t, a, i;
    int64     p;

    T = Open_Kmer_Stream(SOURCE);

    gc_setup(KMER);

    range[0] = 0;
    range[NTHREADS] = T->nels;

    seq = Current_Kmer(T,NULL);
    ent = Current_Entry(T,NULL);
    for (t = 1; t < NTHREADS; t++)
      { p = (T->nels*t)/NTHREADS; 
        GoTo_Kmer_Index(T,p);
#ifdef DEBUG
        printf("\n%d: %0*x\n",t,2*T->ibyte,T->cpre);
        printf(" %lld: %s\n",p,Current_Kmer(T,seq));
#endif
        ent = Current_Entry(T,ent);                //  Break at prefix boundaries
        for (i = T->ibyte; i < T->kbyte; i++)
          ent[i] = 0;
        GoTo_Kmer_Entry(T,ent);
#ifdef DEBUG
        printf(" %lld: %s\n",T->cidx,Current_Kmer(T,seq));
#endif
        range[t] = T->cidx;
      }
    free(seq);

    for (t = 0; t < NTHREADS; t++)
      { parm[t].tid   = t;
        parm[t].T     = T;
        parm[t].beg   = range[t];
        parm[t].end   = range[t+1];
      }

#ifdef DEBUG_THREADS
    for (t = 0; t < NTHREADS; t++)
      merge_thread(parm+t);
#else
    for (t = 1; t < NTHREADS; t++)
      pthread_create(threads+t,NULL,merge_thread,parm+t);
    merge_thread(parm);
    for (t = 1; t < NTHREADS; t++)
      pthread_join(threads[t],NULL);
#endif

    { int64 *plot0, *plott;

      for (t = 1; t < NTHREADS; t++)
        for (i = 0; i <= KMER; i++)
          { plot0 = parm[0].plot[i];
            plott = parm[t].plot[i];
            for (a = 0; a <= HMAX; a++)
              plot0[a] += plott[a];
          }

      for (t = 1; t < NTHREADS; t++)
        { free(parm[t].plot[0]);
          free(parm[t].plot);
        }
    }

    PLOT = parm[0].plot;

    Free_Kmer_Stream(T);
  }

  //  Find peak away from 0 in x-dimension

  { int   i, k, xm, xmax;
    int64 ym, *row;
 
    ZMAX = 0;
    xmax = 0;
    for (i = 0; i <= KMER; i++)
      { row = PLOT[i];
        for (k = 2; row[k] < row[k-1]; k++)
          if (k >= HMAX)
            break;
        ym = row[k];
        xm = k;
        for ( ; k <= HMAX; k++)
          if (row[k] >= ym)
            { ym = row[k];
              xm = k;
            }
        if (ym > ZMAX)
          { ZMAX = ym;
            xmax = xm;
          }
      }

    if (xmax == 0 || xmax >= HMAX)
      { fprintf(stderr,"%s: No maximal peak away from 0 in histogrom interval [1,%d] of %s\n",
                       Prog_Name,HMAX,SOURCE);
        exit (1);
      }

    if (XMAX == 0)
      { XMAX = XREL*xmax;
        if (XMAX > HMAX)
          XMAX = HMAX;
      }
  }

  //  Plot matrix

  { char  *troot;
    int    i, a;
    int64 *row0, *row1, val;
    FILE  *f;
    char   command[1000], *capend;

    //  Output matrix to temp file

    troot = mktemp(template);

    f = fopen(Catenate(troot,".kgc","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".kgc","",""));
    else
      printf("Writing %s\n",Catenate(troot,".kgc","",""));
    fflush(stdout);
#endif
    fprintf(f,"GCP\tKF\tCount\n");
    for (i = 0; i < KMER; i++)
      { row0 = PLOT[i];
        row1 = PLOT[i+1];
        for (a = 0; a < XMAX; a++)
          { val = (row0[a]+row0[a+1]+row1[a]+row1[a+1])/4;
            if (val > ZMAX)
              fprintf(f,"%d.5\t%d.5\t%lld\n",i,a,ZMAX);
            else
              fprintf(f,"%d.5\t%d.5\t%lld\n",i,a,val);
          }
      }
    fclose(f);

    //  Generate the R plot script in another temp file

    f = fopen(Catenate(troot,".R","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".R","",""));
    else
      printf("Generating %s\n",Catenate(troot,".R","",""));
    fflush(stdout);
#endif
    fwrite(kgc_plot,strlen(kgc_plot),1,f);
    fclose(f);

    //  Call the R plotter with arguments

    sprintf(command,"Rscript %s.R -f %s.kgc -o %s%s -x %g -y %g -s %s",
                    troot,troot,OUT,PDF?" -p":" ",XDIM,YDIM,SOURCE);
    capend = command+strlen(command);
    if (LINE)
      { sprintf(capend," -t contour 2>/dev/null");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        SystemX(command);
      }
    if (FILL)
      { sprintf(capend," -t heat 2>/dev/null");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        SystemX(command);
      }
    if (BOTH)
      { sprintf(capend," -t combo 2>/dev/null");
#ifdef DEBUG
        printf("%s\n",command);
        fflush(stdout);
#endif
        SystemX(command);
      }

    //  Remove the temp files

    sprintf(command,"rm -f %s.kgc %s.R",troot,troot);
    SystemX(command);
  }

  free(SOURCE);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
