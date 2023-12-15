/******************************************************************************************
 *
 *  smu2ascii:  Read a .smu file and output it in ascii to standard output.
 *
 *  Author:  Gene Myers
 *  Date  :  November, 2023
 *
 ********************************************************************************************/

#define SMAX  1000    //  Max. value of CovA+CovB
#define FMAX   500    //  Max. value of min(CovA,CovB)

#include <stdlib.h>
#include <stdio.h>

typedef long long int64;

int main(int argc, char *argv[])
{ FILE   *f;
  int64 **PLOT;
  int     i, a, nels;

  f = fopen(argv[1],"r");
  if (f == NULL)
    { fprintf(stderr,"cmu2ascii: Cannot open %s\n",argv[1]);
      exit (1);
    }

  PLOT = malloc(sizeof(int64 *)*(SMAX+1));
  if (PLOT == NULL)
    { fprintf(stderr,"cmu2ascii: Out of memory allocating smu array");
      exit (1);
    }
  PLOT[0] = malloc(sizeof(int64)*(SMAX+1)*(FMAX+1));
  if (PLOT[0] == NULL)
    { fprintf(stderr,"cmu2ascii: Out of memory allocating smu array");
      exit (1);
    }
  for (a = 1; a <= SMAX; a++)
    PLOT[a] = PLOT[a-1] + (FMAX+1);
  nels = (SMAX+1)*(FMAX+1);
  if (fread(PLOT[0],sizeof(int64),nels,f) != nels)
    { fprintf(stderr,"cmu2ascii: Not correct amount of data, really a .smu file?\n");
      exit (1);
    }

  fclose(f);

  printf("// %dx%d matrix, the i'th number in the j'th row give the number of hetmer pairs (a,b)",
         SMAX,FMAX);
  printf("//                     s.t. count(a)+count(b) = j+1 and min(count(a),count(b)) = i+1.");
  for (a = 0; a <= SMAX; a++)
    { for (i = 0; i < FMAX; i++)
        printf(" %lld",PLOT[a][i]);
      printf("\n");
    }

  exit (0);
}

