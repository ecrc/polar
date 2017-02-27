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

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif


 /*******************************************************************************
 *
 *  VERBOSE (global input) INT*1
 *
 *  PROF    (global input) INT*1 
 *
 *  M       (global input) INT*1
 *          Specifies the number of the rows
 *
 *  N       (global input) INT*1
 *          Specifies the number of the columns
 *
 *  OPTCOND (global input) INT*1
 *          = 1:  the condition number can be calculated using the QR 
 *          = 0:  the condition number can be calculated using the LU 
 *
 *  U       (local input/output) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *          On entry, this array contains the matrix to be factorized 
 *          On exit, it contain the orthogonal polar factor U_P
 *
 *  DESCU   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix U.
 *
 *  A       (local workspace) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *
 *  DESCA   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix A.
 *
 *  B       (local workspace) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *
 *  DESCB   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix B.
 *
 *  LDW     (global input) INT*1
 *          Specifies the size of the workspace
 *
 *  C       (local output) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *          On exit, this array contains the symmetric positive semidefinite polar factor H 
 *
 *  DESCC   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix C.
 *
 *  TAU    (local workspace/output) DOUBLE PRECISION   array, dimension
 *          (NLOC)
 *
 *  WORK    (local workspace) DOUBLE PRECISION   array, dimension
 *          (LWORK)
 *
 * LWORK    (local input) INTEGER
 *          The dimension of the array WORK.
 *
 *  IW   (local workspace) INTEGER array, dimension (LWI)
 *          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
 *
 *  LWI  (input) INTEGER
 *          The dimension of the array IW.
 *          PDGETRF_LIWORK = ( LOCr(M_A)+MB_A )
 *          PDGECON_LIWORK >= MAX( 1, LOCr(N+MOD(IA-1,MB_A)) ).
 *          LWI = MAX(PDGETRF_LIWORK, PDGECON_LIWORK ) 
 *
 *
 * ICTXT (global) DESCA( CTXT_ ) The BLACS context handle, indicating
 *                                 the BLACS process grid A is distribu-
 *                                 ted over. The context itself is glo-
 *                                 bal, but the handle (the integer
 *                                 value) may vary.
 *
 * FLOPS (input/output)
 *       The flop count of QDWH
 *
 ******************************************************************************/

int    init = 0;
double eps;
double tol1;
double tol3;

int pdgeqdwh( int verbose, int prof, int M, int N, int optcond, 
	      double *U, int descU[9], double *A, int descA[9], 
              double *B, int descB[9], long int ldw, 
              double *C, int descC[9], double *tau, 
              double *Work, int lWork,
              int *Wi, int lWi, 
              int ictxt, double *flops)
{

    complex dd, sqd, a1;
    double conv = 100.;
    double a, b, c, L2, Liconv, alpha, beta, Anorm, Ainvnorm, Li, norm_est;
    double tol = 3.e-1;
    double flops_dgeqrf, flops_dorgqr, flops_dgemm, flops_dpotrf, flops_dtrsm;
    long int matsize;
    int MB = 2*M;
    int it, itconv, info, facto = -1;
    int itqr = 0, itpo =0;
    int i1 =1, iM = M+1;
    int myrank_mpi;
    double qwtime, litime, nrmtime, potime, qrtime, Htime;
    double sync_time_elapsed, reduced_time_elapsed;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);

    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Entering QDWH\n");}
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Preparing workspace for QDWH\n");}
    qwtime = 0.0;
    if(prof) {qwtime -= MPI_Wtime();}

    if (!init) {
	eps  = pdlamch_( &ictxt, "E" ); 
	tol1 = 5. * eps;
	tol3 = pow(tol1, 1./3.);
	init = 1;
    }

    if ( M < N ){
	fprintf(stderr, "error(m >= n is required)") ;
	return -1;
    }

    /**
     * Create the required workspaces
     */
    matsize = M*N;
    if ( B == NULL ) {
	matsize *= sizeof(double);
	B  = (double *)malloc(2*matsize);
    }
    else {
	if( ldw < (2 * matsize) ) {
	    fprintf(stderr, "Providing workspace is too small: %ld instead of %ld elements\n",
		    ldw, 2*matsize );
	    exit(-1);
	}
    }

    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Finish preparing workspace for QDWH\n");}

    /*
     * Save copy of U in A ==> H = U'*A
     */
    pdlacpy_ ( "A", &M, &N, U, &i1, &i1, descU, A, &i1, &i1, descA );

    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Cond estimate starts\n");}
    /*
     * Calculate Li: reciprocal of condition number estimation
     */

    litime = 0.0;
    if(prof) {litime =- MPI_Wtime();}

    pdlacpy_ ( "A", &M, &N, U, &i1, &i1, descU, B, &i1, &i1, descB );
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "lacpy ends\n");}
    Anorm = pdlange_ ( "1", &M, &N, U, &i1, &i1, descU, Work);
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "dlange ends\n");}

    alpha = 1.0; 
    pdgenm2( U, M, N, descU, B, descB, C, descC, &norm_est, tol);
    pdlascl_( "G", &norm_est, &alpha, &M, &N, U, &i1, &i1, descU, &info);
    //pdlascl_( "G", &alpha, &norm_est, &M, &N, U, &i1, &i1, descU, &info);


    /* estimate condition number using QR */
    if (optcond){
        pdgeqrf_(&M, &N, B, &i1, &i1, descB, tau, Work, &lWork, &info);

        sync_time_elapsed =- MPI_Wtime();
        pdtrtri_( "U", "N", &N, B, &i1, &i1, descB, &info );
        sync_time_elapsed += MPI_Wtime();
        MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        Ainvnorm = pdlange_ ( "1", &M, &N, B, &i1, &i1, descB, Work);
        Li = ( 1.0 / Ainvnorm)/Anorm;    
        Li = norm_est/1.1*Li;    
        *flops += FLOPS_DGEQRF( M, N )
               + FLOPS_DTRTRI(  N );
    }
    /* estimate condition number using LU */
    else {
        pdgetrf_ ( &M, &N, B, &i1, &i1, descB, Wi, &info );
        if (verbose & myrank_mpi == 0) { fprintf(stderr, "LU ends\n");}
        pdgecon_ ("1", &M, B, &i1, &i1, descB, &Anorm, &Li, Work, &lWork, Wi, &lWi, &info);
        Li = norm_est/1.1*Li;    
        /**
         * WARNING: The cost of the gecon is estimated with only one iteration
         */
        *flops += FLOPS_DGETRF(N, N)
               + 2. * FLOPS_DTRSM( 'L', N, 1 );
    }

    if(prof) {litime += MPI_Wtime();}
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Cond estimate ends\n");}

    /*
     * Calculate norm_est
     * Scal the matrix by norm_est
     */

    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Normest starts\n");}
    nrmtime = 0.0;
    if(prof) {nrmtime =- MPI_Wtime();}


    if(prof) {nrmtime += MPI_Wtime();}
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Normest ends\n");}


    itconv = 0; Liconv = Li;

    int itcqr = 0, itcpo = 0;     
    while(itconv == 0 || fabs(1-Liconv) > tol1 ) {
	/* To find the minimum number of iterations to converge. 
         * itconv = number of iterations needed until |Li - 1| < tol1 
	 * This should have converged in less than 50 iterations
         */
	if (itconv > 100) {
	    exit(-1);
	    break;
	}
	itconv++;

	L2  = Liconv * Liconv;
	dd  = cpow( 4. * (1. - L2 ) / (L2 * L2), 1./3. );
	sqd = sqrt(1. + dd);
	a1  = sqd + sqrt( 8. - 4. * dd + 8. * (2. - L2) / (L2 * sqd) ) / 2.;
	a   = creal(a1);
	b   = (a - 1.) * (a - 1.) / 4.;
	c   = a + b - 1.;
        if (c > 100) {itcqr += 1;} else {itcpo += 1;}
	// Update Liconv
	Liconv  = Liconv * (a + b * L2) / (1. + c * L2);
    }
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "QDWH loop starts\n");}
    if (myrank_mpi == 0) { fprintf(stderr, "\nItConv %d itcqr %d itcpo %d norm_est %2.4e Li %2.4e \n", itconv, itcqr, itcpo, norm_est, Li); fprintf(stderr, "It Facto Conv\n");}
    it = 0;
    while(conv > tol3 || it < itconv ) {
	/* This should have converged in less than 50 iterations */
	if (it > 100) {
	    exit(-1);
	    break;
	}
	it++;

	/* Copy U into B1 */
        //pdlacpy_( "A", &M, &N, U, &i1, &i1, descU, C, &i1, &i1, descC );

	// Compute parameters L,a,b,c (second, equivalent way).
	L2  = Li * Li;
	dd  = cpow( 4. * (1. - L2 ) / (L2 * L2), 1./3. );
	sqd = sqrt(1. + dd);
	a1  = sqd + sqrt( 8. - 4. * dd + 8. * (2. - L2) / (L2 * sqd) ) / 2.;
	a   = creal(a1);
	b   = (a - 1.) * (a - 1.) / 4.;
	c   = a + b - 1.;
	// Update Li
	Li  = Li * (a + b * L2) / (1. + c * L2);

	if ( c > 100) {

            qrtime = 0.0;
            if(prof) {qrtime =- MPI_Wtime();}

	    /* Copy U into C to check the convergence of QDWH */
            if (it >= itconv ){
            pdlacpy_( "A", &M, &N, U, &i1, &i1, descU, C, &i1, &i1, descC );
            }

	    /**
	     * Generate the matrix B = [ B1 ] = [ sqrt(c) * U ]
	     *                         [ B2 ] = [ Id          ]
	     */
            pdlacpy_( "A", &M, &N, U, &i1, &i1, descU, B, &i1, &i1, descB );
            alpha = 1.0; beta = sqrt(c);
            pdlascl_( "G", &alpha, &beta, &M, &N, B, &i1, &i1, descB, &info);
            alpha = 0.; beta =1.; 
            pdlaset_( "G", &M, &N, &alpha, &beta, B, &iM, &i1, descB);

	    /**
	     * Factorize B = QR, and generate the associated Q
	     */
            sync_time_elapsed =- MPI_Wtime();

            pdgeqrf_(&MB, &N, B, &i1, &i1, descB, tau, Work, &lWork, &info);
            pdorgqr_(&MB, &N, &N, B, &i1, &i1, descB, tau, Work, &lWork, &info);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    /**
	     * Gemm to find the conv-norm
	     *  U = ( (a-b/c)/sqrt(c) ) * Q1 * Q2' + (b/c) * U
	     */
            alpha = (a-b/c)/sqrt(c); beta = (b/c);
            pdgemm_( "N", "T", &M, &N, &N, &alpha, B, &i1, &i1, descB, B, &iM, &i1, 
                     descB, &beta, U, &i1, &i1, descU);

            if(prof) {qrtime += MPI_Wtime();}

	    /* Main flops used in this step */
	    flops_dgeqrf = FLOPS_DGEQRF( 2*M, N );
	    flops_dorgqr = FLOPS_DORGQR( 2*M, N, N );
	    flops_dgemm  = FLOPS_DGEMM( M, N, N );
	    *flops += flops_dgeqrf + flops_dorgqr + flops_dgemm;

            itqr += 1;
	    facto = 0;
	}
	else {
	    /**
	     * Compute Q1 = c * U * U' + I
	     */

            potime = 0.0;
            if(prof) {potime =- MPI_Wtime();}

            alpha = 0.; beta =1.; 
            pdlaset_( "G", &M, &N, &alpha, &beta, C, &i1, &i1, descC);

            sync_time_elapsed =- MPI_Wtime();

            pdgemm_( "T", "N", &M, &N, &N, &c, U, &i1, &i1, descU, U, &i1, &i1, 
                     descU, &beta, C, &i1, &i1, descC);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    /**
	     * Solve Q1 x = Q2, with Q2 = U
	    */
            alpha = 1.0; beta = 0.0;
            pdgeadd_( "T", &M, &N, &alpha, U, &i1, &i1, descU, &beta, B, &i1, &i1, descB);

            sync_time_elapsed =- MPI_Wtime();

            pdposv_( "U", &M, &N, C, &i1, &i1, descC, B, &i1, &i1, descB, &info);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    /* Copy U into C to check the convergence of QDWH */
            if (it >= itconv ){
            pdlacpy_( "A", &M, &N, U, &i1, &i1, descU, C, &i1, &i1, descC );
            }

	    /**
	     * Compute U =  (a-b/c) * Q2' + (b/c) * U
	     */
            alpha = (a-b/c); beta = (b/c);
            pdgeadd_ ( "T", &M, &N, &alpha, B, &i1, &i1, descB, &beta, U, &i1, &i1, descU);

            if(prof) {potime += MPI_Wtime();}

	    /* Main flops used in this step */
	    flops_dgemm  = FLOPS_DGEMM( M, N, N );
	    flops_dpotrf = FLOPS_DPOTRF( M );
	    flops_dtrsm  = FLOPS_DTRSM( 'L', M, N );
	    *flops += flops_dgemm + flops_dpotrf + 2. * flops_dtrsm;

            itpo += 1;
	    facto = 1;
        }

	/**
	 * Compute the norm of the symmetric matrix U - B1
	 */
        conv = 10.;
        if(it >= itconv ){
            alpha = 1.0; beta = -1.0;
            pdgeadd_ ( "N", &M, &N, &alpha, U, &i1, &i1, descU, &beta, C, &i1, &i1, descC);

            sync_time_elapsed =- MPI_Wtime();

            conv = pdlange_( "F", &M, &N, C, &i1, &i1, descC, Work);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        if (verbose & myrank_mpi == 0) fprintf(stderr, "%02d %-5s %e\n", it,
	       facto == 0 ? "QR" : "PO", conv );
    }
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "QDWH loop ends\n");}

    /*
     * A = U*H ==> H = U'*A ==> H = 0.5*(H'+H)
     */

     Htime = 0.0;
     if(prof) {Htime =- MPI_Wtime();}

     alpha = 1.0; beta = 0.0;
     pdgemm_( "T", "N", &M, &N, &N, &alpha, U, &i1, &i1, descU, A, &i1, &i1, 
              descA, &beta, C, &i1, &i1, descC);
     pdlacpy_( "A", &M, &N, C, &i1, &i1, descC, B, &i1, &i1, descB );
     alpha = 0.5; 
     pdgeadd_ ( "T", &M, &N, &alpha, B, &i1, &i1, descB, &alpha, C, &i1, &i1, descC);

     if(prof) {Htime += MPI_Wtime();}

     flops_dgemm  = FLOPS_DGEMM( M, N, N );
     *flops += flops_dgemm;


    if ( B == NULL ) {
        free(B);
    }

    if(prof) {qwtime += MPI_Wtime();}

    if (prof && (myrank_mpi == 0)) {
        fprintf(stderr, "# QDWH Profiling \n"); 
        fprintf(stderr, "#\n");
        fprintf(stderr, "# \tN    \ttimeQDWH     \ttimeLi     \ttimeNrm    \ttime1itQR   \t#QR    \ttime1itPO   \t#PO    \ttimeFormH \n");
	fprintf(stderr, "  \t%d \t%2.4e \t%2.4e \t%2.4e \t%2.4e \t%d \t%2.4e \t%d \t%2.4e \n", M, qwtime, litime, nrmtime, qrtime, itqr, potime, itpo, Htime);
    }
    if (myrank_mpi == 0) {
        fprintf(stderr, "#\n");
        fprintf(stderr, "# \t#QR  \t#PO  \n");
	fprintf(stderr, "  \t%d  \t%d \n", itqr, itpo);
    }
    //if (myrank_mpi == 0) fprintf(stderr, "======================> FINISH QDWH\n");
    //free(Wi);
    //free(Work);
    //free(Work2);
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Exiting QDWH\n");}
    return 0;

}
