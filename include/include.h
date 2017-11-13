/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
//#include <plasma.h>
#define USAGE(name, args, details)                  \
  printf(" Proper Usage is : ./exe " args " with\n" \
         "  " name "\n"   \
         details);
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif


static inline double cWtime(void)
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

int pdgeqdwh( int verbose, int prof, int M, int N, int optcond, 
	      double *U, int descU[9], double *A, int descA[9], 
              double *B, int descB[9], long int ldw, 
              double *C, int descC[9], double *tau, 
              double *Work, int lWork,
              int *Wi, int lWi, 
              int ictxt, double *flops);

int pdgeqsvd( char *jobu, char *jobvt, char *eigtype, 
              int nprow, int npcol, int nb, int ictxt,
              int m, int n, 
              double *A, int iA, int jA, int *descA, 
              double *S, 
              double *U,     int iU,     int jU, int *descU,
              double *VT,    int iVT,    int jVT, int *descVT,
              double *Wglo,  int iW,     int jW, int *descWglo,
              double *Wloc,  int lwork,
              int    *iWloc, int liwork, int *info);

void pdgenm2( double *A, int M, int N, int descA[9], 
              double *W, int descW[9], double *Sx, int descSx[9], 
              double *e, double tol);
