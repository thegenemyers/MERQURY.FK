/********************************************************************************************
 *
 *  Command line utility to produce Assembly-spectra plots
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

void asm_plot(char  *OUT, int KEEP, char  **ASM, char **ATABLE, char  *RTABLE, 
              double XDIM, double YDIM,
              double XREL, double YREL,
              int    XMAX, int64  YMAX,
              int    PDF, int ZGRAM, int LINE, int FILL, int STACK,
              int NTHREADS);
