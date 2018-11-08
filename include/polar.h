/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpi.h>
#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#include "myscalapack.h"
#include "flops.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif


int pdgeqdwh( char *jobh, int m, int n,
	      double *A, int iA, int jA, int *descA, 
              double *H, int iH, int jH, int *descH,
              double *Work1, int lWork1, 
              double *Work2, int lWork2, 
              int *info);

int pdgezolopd( char *jobh, int m, int n,
	      double *A, int iA, int jA, int *descA, 
              double *H, int iH, int jH, int *descH,
              double *Work1, int lWork1, 
              double *Work2, int lWork2, 
              int *info);

void pdgenm2( double *A, int M, int N, int descA[9], 
              double *W, int descW[9], double *Sx, int descSx[9], 
              double *e, double tol);

int mellipj( double u, double alpha, 
             double *sn, double *cn, double *dn, 
             double *work);

int mellipke( double alpha,  
              double *k, double *e);

