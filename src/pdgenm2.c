/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file pdgenm2.c
 *
 *  QDWH is a high performance software framework for computing 
 *  the polar decomposition on distributed-memory manycore systems provided by KAUST
 *
 * @version 2.0.0
 * @author Dalal Sukkari
 * @author Hatem Ltaief
 * @date 2017-11-13
 *
 **/

#include "common.h"

void pdgenm2( double *A, int M, int N, int descA[9], double *W, int descW[9], double *Sx, int descSx[9], double *e, double tol)
{
    /*
    *
    *   NORMEST Estimate the matrix 2-norm.
    *   NORMEST(S) is an estimate of the 2-norm of the matrix S.
    *   NORMEST(S,tol) uses relative error tol instead of 1.e-6.
    *   [nrm,cnt] = NORMEST(..) also gives the number of iterations used.
    *
    *   This function is intended primarily for sparse matrices,
    *   although it works correctly and may be useful for large, full
    *   matrices as well.  Use NORMEST when your problem is large
    *   enough that NORM takes too long to compute and an approximate
    *   norm is acceptable.
    *
    */

    int maxiter = 100; /* should never take this many iterations. */
    int i1 = 1;
    int cnt = 0;
    int info;
    double e0, alpha, beta, normx, normSx;
    //x = sum(abs(S),1)';
    /* Since in qdwh.c we are finding norm1 of the matrix A, 
     * then we already have the sum of the columns saved in W
     * otherwise we should call the following pdlange
     */
    //e0 = pdlange_ ( "1", &M, &N, A, &i1, &i1, descA, W); 
    double *w  = (double *)malloc(1*sizeof(double)) ;
    
    *e = pdlange_ ( "f", &N, &i1, W, &i1, &i1, descW, w);
    //pdnrm2_( &N, e, W, &i1 , &i1, descW , &i1);
    
    if (*e == 0){ return;}
    //x = x/e;
    alpha = 1.0;
    pdlascl_( "G", e, &alpha, &N, &i1, W, &i1, &i1, descW, &info);

    e0 = 0;
    while ( (cnt < maxiter) &&
           (fabs((*e) - e0) > (tol * (*e))) )
   {
        e0 = *e; alpha = 1.0; beta = 0.0;
        pdgemv_ ("N", &M, &N, &alpha, A, &i1, &i1, descA, W, &i1, &i1, descW, &i1, &beta, Sx, &i1, &i1, descSx, &i1);
        normSx = pdlange_ ( "f", &N, &i1, Sx, &i1, &i1, descSx, w);

        //if nnz(Sx) == 0
        //    Sx = rand(size(Sx),class(Sx));
        //end

        pdgemv_ ("N", &M, &N, &alpha, A, &i1, &i1, descA, Sx, &i1, &i1, descSx, &i1, &beta, W, &i1, &i1, descW, &i1);
        normx = pdlange_ ( "f", &N, &i1, W, &i1, &i1, descW, w);
   
        *e = normx/normSx;
        pdlascl_( "G", &normx, &alpha, &N, &i1, W, &i1, &i1, descW, &info);
        cnt = cnt+1;
        if ( (cnt >= maxiter) &&
            (fabs((*e) - e0) > (tol * (*e))) ) {
            fprintf(stderr, "normest: didn't converge\n");
        }
    }
    return;
}
