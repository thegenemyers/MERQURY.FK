/********************************************************************************************
 *
 *  Refactoring of Merqury scripts as a command line tool using FastK
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

static char *Usage = " [-v] <read> [<mat> <pat>] <asm1> [<asm2>] <out>";

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

static int   VERBOSE;

/****************************************************************************************
 *
 *  Main Routine
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ char  *READS, *MAT, *PAT, *ASM[2], *OUT;
  
  //  Command line processing

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) eptr;

    ARG_INIT("Merqury");

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc < 4 || argc > 7)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        exit (1);
      }
    else if (VERBOSE)
      { switch (argc)
        { case 4:
            fprintf(stderr,"\n No trio data, single diploid assembly\n");
            break;
          case 5:
            fprintf(stderr,"\n No trio data, two haploid assemblies\n");
            break;
          case 6:
            fprintf(stderr,"\n Trio data, single diploid assembly\n");
            break;
          case 7:
            fprintf(stderr,"\n Trio data, two haploid assemblies\n");
            break;
          default:
            ;
        }
      }

    READS = argv[1];
    if (argc >= 6)
      { MAT = argv[2];
        PAT = argv[3];
      }
    else
      MAT = PAT = NULL;
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

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
