/*******************************************************************************************
 *
 *  KatComp
 *
 *  Author:  Gene Myers
 *  Date  :  April, 2021
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <pthread.h>

#undef   DEBUG
#undef   DEBUG_THREADS

#include "libfastk.h"

#include "kx_plot.R.h"

static char *Usage[] = { " [-w<double(6.0)>] [-h<double(4.5)>]",
                         " [-[xX]<number(x2.1)>] [-[yY]<number(y2.1)>]",
                         " [-lfs] [-pdf] [-T<int(4)>]",
                         " <source1>[.ktab] <source2>[.ktab] <out>"
                       };

static int NTHREADS;
static int KMER;
static int HMAX, JMAX;


/****************************************************************************************
 *
 *  Streaming eval
 *
 *****************************************************************************************/

typedef struct
  { int          tid;
    Kmer_Stream *T;
    Kmer_Stream *U;
    int64       *begs;
    int64       *ends;
    int64      **plot;
  } TP;

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n-- > 0)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static void *merge_thread(void *args)
{ TP *parm = (TP *) args;
  int           tid   = parm->tid;
  Kmer_Stream  *T     = parm->T;
  Kmer_Stream  *U     = parm->U;
  int64        *begs  = parm->begs;
  int64        *ends  = parm->ends;

  int hbyte = T->hbyte;

  int64 **plot = NULL;
  int     i;

#ifdef DEBUG_TRACE
  char *buffer;
#endif

  plot    = Malloc(sizeof(int64 *)*(JMAX+1),"Allocating thread working memory");
  plot[0] = Malloc(sizeof(int64)*(JMAX+1)*(HMAX+1),"Allocating plot");
  for (i = 1; i <= JMAX; i++)
    plot[i] = plot[i-1] + (HMAX+1);
  bzero(plot[0],sizeof(int64)*(HMAX+1)*(JMAX+1));

#ifdef DEBUG_THREADS
  printf("Doing %d:",tid);
  for (c = 0; c < 2; c++)
    printf(" [%lld-%lld]",begs[c],ends[c]);
  printf("\n");
#endif

  if (tid != 0)
    { T = Clone_Kmer_Stream(T);
      U = Clone_Kmer_Stream(U);
    }

#ifdef DEBUG_TRACE
  buffer = Current_Kmer(T,NULL);
#endif

  GoTo_Kmer_Index(T,begs[0]);
  GoTo_Kmer_Index(U,begs[1]);

  while (1)
    { int x, tn, un;

      if (T->cidx < ends[0])
        if (U->cidx < ends[1])
          { if (T->cpre < U->cpre)
              x = -1;
            else if (T->cpre > U->cpre)
              x = 1;
            else
              x = mycmp(T->csuf,U->csuf,hbyte);
          }
        else
          x = -1;
      else
        if (U->cidx < ends[1])
          x = 1;
        else
          break;

      if (x < 0)
        { tn = Current_Count(T);
          if (tn <= HMAX)
            plot[0][tn] += 1;
          Next_Kmer_Entry(T);
        }
      else if (x > 0)
        { un = Current_Count(U);
          if (un <= JMAX)
            plot[un][0] += 1;
          Next_Kmer_Entry(U);
        }
      else
        { tn = Current_Count(T);
          un = Current_Count(U);
          if (tn <= HMAX && un <= JMAX)
            plot[un][tn] += 1;
          Next_Kmer_Entry(T);
          Next_Kmer_Entry(U);
        }
    }

  if (tid != 0)
    { Free_Kmer_Stream(T);
      Free_Kmer_Stream(U);
    }

  parm->plot = plot;
  return (NULL);
}


/****************************************************************************************
 *
 *  Main
 *
 *****************************************************************************************/

static char template[15] = "._KX.XXXX";

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
{ int    LINE, FILL, BOTH;
  int    PDF;
  double XDIM, YDIM;
  double XREL, YREL;
  int    XMAX, YMAX;
  char  *OUT;
  char  *SOURCE1;
  char  *SOURCE2;

  int64 **PLOT;
  int64   ZMAX;

  //  Process the comand line options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("KatComp");

    XDIM = 6.0;
    YDIM = 4.5;
    XREL = 2.1;
    YREL = 2.1;
    XMAX = 0;
    YMAX = 0;
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
            ARG_POSITIVE(YMAX,"y max");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    LINE  = flags['l'];
    FILL  = flags['f'];
    BOTH  = flags['s'];

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
        fprintf(stderr,"      -y: max y as a real-valued multiple of y* with max\n");
        fprintf(stderr,"              count 'peak' away from the origin\n");
        fprintf(stderr,"      -Y: max y as an int value in absolute terms\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: draw line plot\n");
        fprintf(stderr,"      -f: draw fill plot\n");
        fprintf(stderr,"      -s: draw stack plot\n");
        fprintf(stderr,"          any combo allowed, none => draw all\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"    -pdf: output .pdf (default is .png)\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -T: number of threads to use\n");
        exit (1);
      }

    if (LINE+FILL+BOTH == 0)
      LINE = FILL = BOTH = 1;
    SOURCE1 = Root(argv[1],".ktab");
    SOURCE2 = Root(argv[2],".ktab");
    OUT     = argv[3];

    KMER = check_table(Catenate(SOURCE1,".ktab","",""),0);
    KMER = check_table(Catenate(SOURCE2,".ktab","",""),KMER);

    if (XMAX == 0)
      HMAX = 1000;
    else
      HMAX = XMAX;
    if (YMAX == 0)
      JMAX = 1000;
    else
      JMAX = YMAX;
  }

  //  Compute the K-mer cross product matrix (JMAX x HMAX)

  { Kmer_Stream  *T, *U;
    int64     range[NTHREADS+1][2];
#ifndef DEBUG_THREADS
    pthread_t threads[NTHREADS];
#endif
    TP        parm[NTHREADS];
    char     *seq;
    uint8    *ent;
    int       t, a, i;
    int64     p;
  
    T = Open_Kmer_Stream(SOURCE1);
    U = Open_Kmer_Stream(SOURCE2);

    range[0][0] = 0;
    range[NTHREADS][0] = T->nels;
    range[0][1] = 0;
    range[NTHREADS][1] = U->nels;

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
        range[t][0] = T->cidx;
        GoTo_Kmer_Entry(U,ent);
#ifdef DEBUG
        printf(" %lld: %s\n",U->cidx,Current_Kmer(U,seq));
#endif
        range[t][1] = U->cidx;
      }
    free(seq);

    for (t = 0; t < NTHREADS; t++)
      { parm[t].tid   = t;
        parm[t].T     = T;
        parm[t].U     = U;
        parm[t].begs  = range[t];
        parm[t].ends  = range[t+1];
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
        for (i = 0; i <= JMAX; i++)
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
    Free_Kmer_Stream(U);
  }
        
  //  Find peaks away from 0 in both the x- and y-dimensions
    
  { int   i, k, xm, xmax, ymax;
    int64 ym, *row, zmax; 
      
    zmax = 0;
    xmax = 0;
    for (i = 0; i <= JMAX; i++)
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
        if (ym > zmax)
          { zmax = ym;
            xmax = xm;
          }
      }
        
    if (xmax == 0 || xmax >= HMAX)
      { fprintf(stderr,"%s: No maximal peak away from 0 in histogrom interval [1,%d] of %s\n",
                       Prog_Name,HMAX,SOURCE1);
        exit (1);
      }   

    if (XMAX == 0)
      XMAX = XREL*xmax;
    ZMAX = zmax;

    zmax = 0;
    ymax = 0;
    for (i = 0; i <= HMAX; i++)
      { for (k = 2; PLOT[k][i] < PLOT[k-1][i]; k++)
          if (k >= JMAX)
            break;
        ym = PLOT[k][i];
        xm = k;
        for ( ; k <= JMAX; k++)
          if (PLOT[k][i] >= ym)
            { ym = PLOT[k][i];
              xm = k;
            }
        if (ym > zmax)
          { zmax = ym;
            ymax = xm;
          }
      }
        
    if (ymax == 0 || ymax >= JMAX)
      { fprintf(stderr,"%s: No maximal peak away from 0 in histogrom interval [1,%d] of %s\n",
                       Prog_Name,JMAX,SOURCE2);
        exit (1);
      }   

    if (YMAX == 0)
      YMAX = YREL*ymax;
    if (ZMAX < zmax)
      ZMAX = zmax;
  }       

  //  Plot matrix

  { char  *troot;
    int    i, a;
    int64 *row0, *row1, val;
    FILE  *f;
    char   command[1000], *capend;

    //  Output matrix to temp file

    troot = mktemp(template);

    f = fopen(Catenate(troot,".kx","",""),"w");
#ifdef DEBUG
    if (f == NULL)
      printf("Could not open %s\n",Catenate(troot,".kx","",""));
    else
      printf("Writing %s\n",Catenate(troot,".kx","",""));
    fflush(stdout);
#endif
    fprintf(f,"KF1\tKF2\tCount\n");
    for (i = 0; i < YMAX; i++)
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
    fwrite(kx_plot,strlen(kx_plot),1,f);
    fclose(f);

    //  Call the R plotter with arguments

    sprintf(command,"Rscript %s.R -f %s.kx -o %s%s -x %g -y %g -s1 %s -s2 %s",
                    troot,troot,OUT,PDF?" -p":" ",XDIM,YDIM,SOURCE1,SOURCE2);
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

    sprintf(command,"rm -f %s.kx %s.R",troot,troot);
    SystemX(command);
  }

  free(SOURCE2);
  free(SOURCE1);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
