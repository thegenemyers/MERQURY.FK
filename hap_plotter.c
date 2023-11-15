/********************************************************************************************
 *
 *  Routine to produce haplotype blob plots
 *
 *  Author:  Gene Myers
 *  Date  :  September, 2021
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#undef KAMIL   //  if set then instead of plots output R-script tables as OUT.hpi

#undef DEBUG

#include "libfastk.h"

#include "hap_plot.R.h"

void hap_plot(char  *OUT, char  *MAT, char  *PAT, char **ASM,
              double XDIM, double YDIM, int PDF, char  *troot)
        
{ char           command[1000];  //  buffers for commands & parts thereof
  char           name[1000];
  Profile_Index *AI, *MI, *PI;
  uint16        *aprof, *mprof, *pprof;
  int64          pmax, plen;
  int            nct, mct, pct, cln;
  int            i, c, x;
  FILE          *f;

  MAT = Root(MAT,"");
  PAT = Root(PAT,"");

  f = fopen(Catenate(troot,".hpi","",""),"r");
  if (f != NULL)
    goto plot;

  //  Scan profiles to produce blob table

  f = fopen(Catenate(troot,".hpi","",""),"w");
#ifdef DEBUG
  if (f == NULL)
    printf("Could not open %s\n",Catenate(troot,".hpi","",""));
  else
    printf("Writing %s\n",Catenate(troot,".hpi","",""));
  fflush(stdout);
#endif
  fprintf(f,"Assembly\tContig\t%s\t%s\tSize\n",MAT,PAT);

  pmax  = 20000;
  aprof = Malloc(3*pmax*sizeof(uint16),"Profile array");
  if (aprof == NULL)
    exit (1);
  mprof = aprof + pmax;
  pprof = mprof + pmax;

  for (i = 0; i < 2; i++)
    { char *asmx;
      if (ASM[i] == NULL)
        continue;
      asmx = Root(ASM[i],"");
      sprintf(name,"%s.prof",ASM[i]);
      AI = Open_Profiles(name);
      sprintf(name,"%s.%s.prof",asmx,MAT);
      MI = Open_Profiles(name);
      sprintf(name,"%s.%s.prof",asmx,PAT);
      PI = Open_Profiles(name);

      nct = 0;
      for (c = 0; c < AI->nreads; c++)
        { plen = Fetch_Profile(AI,i,pmax,aprof);
          if (plen > pmax)
            { pmax  = 1.2*plen + 1000;
              aprof = Realloc(aprof,3*pmax*sizeof(uint16),"Profile array");
              if (aprof == NULL)
                exit (1);
              mprof = aprof + pmax;
              pprof = mprof + pmax;
              Fetch_Profile(AI,c,pmax,aprof);
            }
          Fetch_Profile(MI,c,pmax,mprof);
          Fetch_Profile(PI,c,pmax,pprof);

          cln = 0;
          mct = pct = 0;
          for (x = 0; x < plen; x++)
            if (aprof[x] == 0)
              { if (cln > 0)
                  { nct += 1;
                    fprintf(f,"%s\t%d\t%d\t%d\t%d\n",asmx,nct,mct,pct,cln);
                    mct = pct = 0;
                    cln = 0;
                  }
              }
            else
              { mct += (mprof[x] > 0);
                pct += (pprof[x] > 0);
                cln += 1;
              }
          if (cln > 0)
            { nct += 1;
              fprintf(f,"%s\t%d\t%d\t%d\t%d\n",asmx,nct,mct,pct,cln);
            }
        }
      free(asmx);
    }

  free(aprof);

plot:
  fclose(f);

#ifdef KAMIL

  (void) XDIM; (void) YDIM; (void) PDF;

  sprintf(command,"mv %s.hpi %s.hpi",troot,OUT);
  system(command);

#else

  //  Generate the R plot script

  f = fopen(Catenate(troot,".R","",""),"w");
#ifdef DEBUG
  if (f == NULL)
    printf("Could not open %s\n",Catenate(troot,".R","",""));
  else
    printf("Generating %s\n",Catenate(troot,".R","",""));
  fflush(stdout);
#endif
  fwrite(hap_plot_script,strlen(hap_plot_script),1,f);
  fclose(f);

  //  Call the R plotter with arguments

  sprintf(command,"Rscript %s.R -f %s.hpi -o %s%s -x %g -y %g 2>/dev/null",
                  troot,troot,OUT,PDF?" -p":" ",XDIM,YDIM);
#ifdef DEBUG
  printf("%s\n",command);
  fflush(stdout);
#endif
  system(command);

  sprintf(command,"rm -f %s.hpi %s.R",troot,troot);
  system(command);

#endif // KAMIL

  free(PAT);
  free(MAT);
}
