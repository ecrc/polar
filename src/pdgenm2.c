/**
 QDWH-SVD
 *
 * (C) Copyright 2016 King Abdullah University of Science and Technology
 * Authors:
 * Dalal Sukkari (dalal.sukkari@kaust.edu.sa)
 * David Keyes (david.keyes@kaust.edu.sa)
 * Hatem Ltaief (hatem.ltaief@kaust.edu.sa)
 *  
 * Redistribution  and  use  in  source and binary forms, with or without
 * modification,  are  permitted  provided  that the following conditions
 * are met:
 * 
 * Redistributions  of  source  code  must  retain  the above copyright
 * notice,  this  list  of  conditions  and  the  following  disclaimer.
 * Redistributions  in  binary  form must reproduce the above copyright
 * notice,  this list of conditions and the following disclaimer in the
 * documentation  and/or other materials provided with the distribution.
 * Neither  the  name of the King Abdullah University of Science and
 * Technology nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior
 * written permission.
 *
 *
 THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
 LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
