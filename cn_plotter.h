/********************************************************************************************
 *
 *  Routine to produce CN-spectra plots
 *
 *  Author:  Gene Myers
 *  Date  :  March, 2021
 *
 *********************************************************************************************/
 
int64 cnplot(char  *OUT, char  *ASM, char  *READS,
             double XDIM, double YDIM,
             double XREL, double YREL,
             int    XMAX, int64  YMAX,
             int    PDF, int ZGRAM, int LINE, int FILL, int STACK,
             char  *troot, int NTHREADS);
