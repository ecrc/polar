/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file pdgeqdwh.c
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

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif


 /*******************************************************************************
 *  .. Scalar Arguments ..
 *  INTEGER            IA, INFO, JA, LWORK, M, N
 *   ..
 *  .. Array Arguments ..
 *  INTEGER            DESCA( * )
    DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
 *     ..
 *
 *  Purpose
 *  =======
 *  
 *  PDGDWH computes the polar decomposition of a real distributed M-by-N  
 *
 *  matrix A = U * H.
 *
 *  Notes
 *  =====
 *
 *  Each global data object is described by an associated description
 *  vector.  This vector stores the information required to establish
 *  the mapping between an object element and its corresponding process
 *  and memory location.
 *
 *  Let A be a generic term for any 2D block cyclicly distributed array.
 *  Such a global array has an associated description vector DESCA.
 *  In the following comments, the character _ should be read as
 *  "of the global array".
 *
 *  NOTATION        STORED IN      EXPLANATION
 *  --------------- -------------- --------------------------------------
 *  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
 *                                 DTYPE_A = 1.
 *  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
 *                                 the BLACS process grid A is distribu-
 *                                 ted over. The context itself is glo-
 *                                 bal, but the handle (the integer
 *                                 value) may vary.
 *  M_A    (global) DESCA( M_ )    The number of rows in the global
 *                                 array A.
 *  N_A    (global) DESCA( N_ )    The number of columns in the global
 *                                 array A.
 *  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
 *                                 the rows of the array.
 *  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
 *                                 the columns of the array.
 *  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
 *                                 row of the array A is distributed.
 *  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
 *                                 first column of the array A is
 *                                 distributed.
 *  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
 *                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
 *
 *  Let K be the number of rows or columns of a distributed matrix,
 *  and assume that its process grid has dimension p x q.
 *  LOCr( K ) denotes the number of elements of K that a process
 *  would receive if K were distributed over the p processes of its
 *  process column.
 *  Similarly, LOCc( K ) denotes the number of elements of K that a
 *  process would receive if K were distributed over the q processes of
 *  its process row.
 *  The values of LOCr() and LOCc() may be determined via a call to the
 *  ScaLAPACK tool function, NUMROC:
 *          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
 *          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
 *  An upper bound for these quantities may be computed by:
 *          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
 *          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
 *
 *  Arguments
 *************
 *  JOBH    (global input) CHARACTER*1
 *          Specifies options for computing H:
 *          = 'H':  the H (the symmetric positive
 *                  semidefinite polar factor) are returned in the array H;
 *          = 'N':  no columns of H (no symmetric positive semidefinite polar factor) are
 *                  computed.
 *
 *  M       (global input) INTEGER
 *          The number of rows of the input matrix A.  M >= 0.
 *
 *  N       (global input) INTEGER
 *          The number of columns of the input matrix A.  N >= 0.
 *
 *  A       (local input/output) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *          On entry, this array contains the matrix to be factorized 
 *          On exit, it contain the orthogonal polar factor A_P
 *
 *  IA      (global input) INTEGER
 *          The row index in the global array A indicating the first
 *          row of sub( A ).
 *
 *  JA      (global input) INTEGER
 *          The column index in the global array A indicating the
 *          first column of sub( A ).
 *
 *  DESCA   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix A.
 *
 *  H       (local output) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *          On exit, this array contains the symmetric positive semidefinite polar factor H 
 *
 *  IH      (global input) INTEGER
 *          The row index in the global array H indicating the first
 *          row of sub( H ).
 *
 *  JH      (global input) INTEGER
 *          The column index in the global array H indicating the
 *          first column of sub( H ).
 *
 *  DESCH   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix H.
 *
 *  WORK    (local workspace/output) DOUBLE PRECISION   array, dimension
 *          (LWORK* NQ)
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 * LWORK    (local input) INTEGER
 *          The dimension of the array WORK.
 *          LWORK must be at least LWORK >= 3*MP*NQ
 *          If LWORK = -1, then LWORK is global input and a workspace    
 *          query is assumed; the routine only calculates the minimum
 *          and optimal size for all work arrays. Each of these
 *          values is returned in the first entry of the corresponding
 *          work array, and no error message is issued by PXERBLA.         
 *
 *  INFO (global output) INTEGER
 *          = 0:  successful exit
 *          < 0:  If the i-th argument is an array and the j-entry had
 *                an illegal value, then INFO = -(i*100+j), if the i-th
 *                argument is a scalar and had an illegal value, then
 *                INFO = -i.
 *
 ******************************************************************************/

int    init = 0;
double eps;
double tol1;
double tol3;

int pdgeqdwh( char *jobh, int m, int n,
	      double *A, int iA, int jA, int *descA, 
              double *H, int iH, int jH, int *descH,
              double *Work1, int lWork1, 
              double *Work2, int lWork2, 
              int *info)
{

    complex dd, sqd, a1;
    double conv = 100.;
    double a, b, c, L2, Liconv, alpha, beta, Anorm, Ainvnorm, Li, norm_est;
    double tol = 3.e-1;
    double flops_dgeqrf, flops_dorgqr, flops_dgemm, flops_dpotrf, flops_dtrsm;
    long int matsize;
    int MB = 2*m;
    int it, itconv, facto = -1;
    int itqr = 0, itpo = 0, alloc_qr = 0;
    int i1 =1, i0 = 0, iM = m+1;
    int myrank_mpi;
    double qwtime, litime, nrmtime, potime, qrtime, Htime;
    double sync_time_elapsed, reduced_time_elapsed;

    int verbose = 0, prof = 0, optcond = 0;
    double flops;
    flops = 0.0;

    int mloc, nloc, mlocW, nb;   
    int myrow, mycol, nprow, npcol;   
    int ctxt_ = 1, nb_ = 5;
    int ictxt;
    int wantH;

    int lWi, lwork_qr, lwork_cn;
    int *Wi = (int *)calloc(1,sizeof(int)) ;
    double *W  = (double *)calloc(1,sizeof(double)) ;

    int iinfo;

    /*
     * Get the grid parameters
     */
    ictxt = descA[ctxt_];
    Cblacs_get( -1, 0, &ictxt );
    nb = descA[nb_];
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    mloc  = numroc_( &m, &nb, &myrow, &i0, &nprow );
    nloc  = numroc_( &n, &nb, &mycol, &i0, &npcol );
    mlocW = numroc_( &MB, &nb, &myrow, &i0, &nprow );
 
    int lmin1, lmin2, lquery;
    *info = 0; 
    lquery =  (lWork1 == -1 || lWork2 == -1); 
    wantH = 0;

   /*
    * Test the input parameters
    */
    if( nprow == -1 ){
       *info = -(700+ctxt_);
    }
    else { 
       if ( m < n ){
	   fprintf(stderr, "error(m >= n is required)") ;
	   return -1;
       }
       if (jobh[0] == 'H' || jobh[0] == 'h'){
           wantH = 1;
       }

       int i2 = 2, i3 = 3, i7 = 7, i11 = 11, i_1 = -1;
       int *idum1, *idum2;
       idum1 = (int *)malloc(2*sizeof(int)) ;
       idum2 = (int *)malloc(2*sizeof(int)) ;
       chk1mat_(&m, &i2, &n, &i3, &iA, &jA, descA, &i7, info);
       if (wantH){
          chk1mat_(&m, &i2, &n, &i3, &iH, &jH, descH, &i11, info);
       }
       //igamx2d_(descA[ctxt_], "A", " ", &i1, &i1, info, &i1, &i1, &i1, &i_1, &i_1, &i0);

       lquery =  (lWork1 == -1 || lWork2 == -1); 
       if (*info == 0){
          lmin1 = mloc;
          lmin2 = mlocW;
          Work1[0] = lmin1;
          Work2[0] = lmin2;
          lquery =  (lWork1 == -1 || lWork2 == -1); 
          if( (lWork1 < lmin1) & !lquery ){
              *info = -13;
          }
          if( (lWork2 < lmin2) & !lquery ){
              *info = -15;
          }
       }

       idum1[0] = wantH;
       if( lWork1 == -1 || lWork2 == -1) {
             idum1[1] = -1;
       }
       else {
             idum1[1] =  1;
       }
       idum2[0] =  1;
       idum2[1] =  15;
       pchk1mat_( &m, &i2, &n, &i3, &iA, &jA, descA, &i7, &i2, &idum1, &idum2,
                        info );
       if ((*info == 0) && wantH){
          pchk1mat_( &m, &i2, &n, &i3, &iH, &jH, descH, &i11, &i0, &idum1, &idum2,
                        info );
       }
    }
      
    if( *info != 0 ){
        pxerbla_( ictxt, "PDGEQDWH", -1*info[0] ); 
        return 0;
    }
    else if ( lquery ){
    //lquery =  (lWork1 == -1 || lWork2 == -1); 
    //if ( lquery ){
    /*
     * Find Workspace 
     */
      /*
       int lwork_qr = -1, lwork_cn = -1;
       pdgecon_ ("1", &m, H, &iH, &jH, descH, 
                 &Anorm, &Li, 
                 Work, &lwork_cn, Wi, &lWi, info);
       lwork_cn  = (int)Work[0];
       lWi = N;//(int)iWloc[0];

       pdgeqrf_(&MB, &n, H, &iH, &iH, descH, 
                tau, Work, &lwork_qr, info);
       lwork_qr  = Work[0];
       lWork  = max ( lwork_cn, lwork_qr);
       lWi = N;
       */
       //Work[0]  = mlocW_; //B = Work;
       //Work[0] = 3*mloc;
       //Work[0] = 3*m;
       //Work[0] = ((3*m+nb)/nb)*nb;
       Work1[0] = mloc;
       Work2[0] = mlocW;
       return 0;
    } 

    /* Quick return if possible */
    if ( m == 0 || n == 0 ){
         return 0;
    }

    /**
     * Create the required workspaces
     * Needed for debugging the code
     */

    double *U=NULL, *B=NULL;
    int descU[9], descB[9];

    //int MB3 = 3*m;
    //int mlocW3 = numroc_( &MB3, &nb, &myrow, &i0, &nprow );
    if ( Work1 == NULL ) {
	U  = (double *)malloc(mloc*nloc*sizeof(double));
	//B  = (double *)malloc(mlocW*nloc*sizeof(double));
    }
    if ( Work2 == NULL ) {
	//A  = (double *)malloc(mloc*nloc*sizeof(double));
	B  = (double *)malloc(mlocW*nloc*sizeof(double));
    }
    else {
        U = Work1;
        //B = A + mlocW3*nloc;
        B = Work2;
    }

    descinit_( descU, &m, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &iinfo );
    descinit_( descB, &MB, &n, &nb, &nb, &i0, &i0, &ictxt, &mlocW, &iinfo ); //B = A + mloc*nloc;

    //lWork = 3*m; MB = 2*m;
    //descinit_( descA, &lWork, &n, &nb, &nb, &MB, &i0, &ictxt, &mloc, &iinfo );
    //descinit_( descB, &lWork, &N, &nb, &nb, &i0, &i0, &ictxt, &mlocW, &iinfo ); //B = Work;

    double *tau   = (double *)malloc(nloc*sizeof(double)) ;

    if ( !optcond ){
        lWi = n; 
        Wi  = (int *)malloc(lWi*sizeof(int)) ;
    }

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


    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Finish preparing workspace for QDWH\n");}

    /*
     * Save copy of A ==> H = U'*A
     */
    pdlacpy_ ( "A", &m, &n, A, &i1, &i1, descA, U, &i1, &i1, descU );

    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Cond estimate starts\n");}
    /*
     * Calculate Li: reciprocal of condition number estimation
     */

    litime = 0.0;
    if(prof) {litime =- MPI_Wtime();}

    pdlacpy_ ( "A", &m, &n, A, &i1, &i1, descA, B, &i1, &i1, descB );
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "lacpy ends\n");}
    Anorm = pdlange_ ( "1", &m, &n, U, &i1, &i1, descU, H);
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "dlange ends\n");}

    alpha = 1.0; 
    pdgenm2( U, m, n, descU, B, descB, H, descH, &norm_est, tol);
    pdlascl_( "G", &norm_est, &alpha, &m, &n, A, &i1, &i1, descA, &iinfo);
    //pdlascl_( "G", &alpha, &norm_est, &m, &n, A, &i1, &i1, descA, &iinfo);


    /* estimate condition number using QR */
    if ( optcond ){
        pdgeqrf_(&m, &n, B, &i1, &i1, descB, tau, H, &lWork1, &iinfo);

        sync_time_elapsed =- MPI_Wtime();
        pdtrtri_( "U", "N", &n, B, &i1, &i1, descB, &iinfo );
        sync_time_elapsed += MPI_Wtime();
        MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        Ainvnorm = pdlange_ ( "1", &m, &n, B, &i1, &i1, descB, H);
        Li = ( 1.0 / Ainvnorm)/Anorm;    
        Li = norm_est/1.1*Li;    
        flops += FLOPS_DGEQRF( m, n )
               + FLOPS_DTRTRI(  n );
    }
    /* estimate condition number using LU */
    else {
        pdgetrf_ ( &m, &n, B, &i1, &i1, descB, Wi, &iinfo );
        if (verbose & myrank_mpi == 0) { fprintf(stderr, "LU ends\n");}

        int lwork_cn = -1;
        pdgecon_ ("1", &m, B, &i1, &i1, descB, &Anorm, &Li, H, &lwork_cn, Wi, &lWi, &iinfo);
        lwork_cn = H[0];

        pdgecon_ ("1", &m, B, &i1, &i1, descB, &Anorm, &Li, H, &lwork_cn, Wi, &lWi, &iinfo);
        Li = norm_est/1.1*Li;    
        /**
         * WARNING: The cost of the gecon is estimated with only one iteration
         */
        flops += FLOPS_DGETRF(n, n)
               + 2. * FLOPS_DTRSM( 'L', n, 1 );
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
        if (c > 100) {itcqr += 1; alloc_qr = 1;} else {itcpo += 1;}
	// Update Liconv
	Liconv  = Liconv * (a + b * L2) / (1. + c * L2);
    }
    if (verbose & myrank_mpi == 0) { fprintf(stderr, "QDWH loop starts\n");}
    if (myrank_mpi == 0) { fprintf(stderr, "\nItConv %d itcqr %d itcpo %d norm_est %2.4e Li %2.4e \n", itconv, itcqr, itcpo, norm_est, Li); fprintf(stderr, "It Facto Conv\n");}
    it = 0;

    if ( alloc_qr ){
         lwork_qr = -1;
         pdgeqrf_(&MB, &n, B, &i1, &i1, descB, 
                  tau, W, &lwork_qr, &iinfo);
         lwork_qr  = W[0];
         W  = (double *)malloc((lwork_qr)*sizeof(double)) ;
    }


    while(conv > tol3 || it < itconv ) {
	/* This should have converged in less than 50 iterations */
	if (it > 100) {
	    exit(-1);
	    break;
	}
	it++;

	/* Copy U into B1 */
        //pdlacpy_( "A", &m, &n, U, &i1, &i1, descU, C, &i1, &i1, descC );

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
                pdlacpy_( "A", &m, &n, A, &i1, &i1, descA, H, &i1, &i1, descH );
            }

	    /**
	     * Generate the matrix B = [ B1 ] = [ sqrt(c) * U ]
	     *                         [ B2 ] = [ Id          ]
	     */
            pdlacpy_( "A", &m, &n, A, &i1, &i1, descA, B, &i1, &i1, descB );
            alpha = 1.0; beta = sqrt(c);
            pdlascl_( "G", &alpha, &beta, &m, &n, B, &i1, &i1, descB, &iinfo);
            alpha = 0.; beta =1.; 
            pdlaset_( "G", &m, &n, &alpha, &beta, B, &iM, &i1, descB);

	    /**
	     * Factorize B = QR, and generate the associated Q
	     */
            sync_time_elapsed =- MPI_Wtime();

            pdgeqrf_(&MB, &n, B, &i1, &i1, descB, tau, W, &lwork_qr, &iinfo);
            pdorgqr_(&MB, &n, &n, B, &i1, &i1, descB, tau, W, &lwork_qr, &iinfo);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    /**
	     * Gemm to find the conv-norm
	     *  U = ( (a-b/c)/sqrt(c) ) * Q1 * Q2' + (b/c) * U
	     */
            alpha = (a-b/c)/sqrt(c); beta = (b/c);
            pdgemm_( "N", "T", &m, &n, &n, &alpha, B, &i1, &i1, descB, B, &iM, &i1, 
                     descB, &beta, A, &i1, &i1, descA);

            if(prof) {qrtime += MPI_Wtime();}

	    /* Main flops used in this step */
	    flops_dgeqrf = FLOPS_DGEQRF( 2*m, n );
	    flops_dorgqr = FLOPS_DORGQR( 2*m, n, n );
	    flops_dgemm  = FLOPS_DGEMM( m, n, n );
	    flops += flops_dgeqrf + flops_dorgqr + flops_dgemm;

            itqr += 1;
	    facto = 0;
	}
	else {
	    /**
	     * Compute Q1 = c * U * U' + I
	     */
            alloc_qr = 0;
            potime = 0.0;
            if(prof) {potime =- MPI_Wtime();}

            alpha = 0.; beta =1.; 
            pdlaset_( "G", &m, &n, &alpha, &beta, H, &i1, &i1, descH);

            sync_time_elapsed =- MPI_Wtime();

            pdgemm_( "T", "N", &m, &n, &n, &c, A, &i1, &i1, descA, A, &i1, &i1, 
                     descA, &beta, H, &i1, &i1, descH);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    /**
	     * Solve Q1 x = Q2, with Q2 = U
	     */
            alpha = 1.0; beta = 0.0;
            pdgeadd_( "T", &m, &n, &alpha, A, &i1, &i1, descA, &beta, B, &i1, &i1, descB);

            sync_time_elapsed =- MPI_Wtime();

            pdposv_( "U", &m, &n, H, &i1, &i1, descH, B, &i1, &i1, descB, &iinfo);

            sync_time_elapsed += MPI_Wtime();
            MPI_Allreduce( &sync_time_elapsed, &reduced_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    /* Copy U into H to check the convergence of QDWH */
            if (it >= itconv ){
                pdlacpy_( "A", &m, &n, A, &i1, &i1, descA, H, &i1, &i1, descH );
            }

	    /**
	     * Compute U =  (a-b/c) * Q2' + (b/c) * U
	     */
            alpha = (a-b/c); beta = (b/c);
            pdgeadd_ ( "T", &m, &n, &alpha, B, &i1, &i1, descB, &beta, A, &i1, &i1, descA);

            if(prof) {potime += MPI_Wtime();}

	    /* Main flops used in this step */
	    flops_dgemm  = FLOPS_DGEMM( m, n, n );
	    flops_dpotrf = FLOPS_DPOTRF( m );
	    flops_dtrsm  = FLOPS_DTRSM( 'L', m, n );
	    flops += flops_dgemm + flops_dpotrf + 2. * flops_dtrsm;

            itpo += 1;
	    facto = 1;
        }

	/**
	 * Compute the norm of the symmetric matrix U - B1
	 */
        conv = 10.;
        if(it >= itconv ){
            alpha = 1.0; beta = -1.0;
            pdgeadd_ ( "N", &m, &n, &alpha, A, &i1, &i1, descA, &beta, H, &i1, &i1, descH);

            sync_time_elapsed =- MPI_Wtime();

            conv = pdlange_( "F", &m, &n, H, &i1, &i1, descH, W);

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
     pdgemm_( "T", "N", &m, &n, &n, &alpha, A, &i1, &i1, descA, U, &i1, &i1, 
              descU, &beta, H, &i1, &i1, descH);
     pdlacpy_( "A", &m, &n, H, &i1, &i1, descH, B, &i1, &i1, descB );
     alpha = 0.5; 
     pdgeadd_ ( "T", &m, &n, &alpha, B, &i1, &i1, descB, &alpha, H, &i1, &i1, descH);

     if(prof) {Htime += MPI_Wtime();}

     flops_dgemm  = FLOPS_DGEMM( m, n, n );
     flops += flops_dgemm;


    if(prof) {qwtime += MPI_Wtime();}

    if (prof && (myrank_mpi == 0)) {
        fprintf(stderr, "# QDWH Profiling \n"); 
        fprintf(stderr, "#\n");
        fprintf(stderr, "# \tn    \ttimeQDWH     \ttimeLi     \ttimeNrm    \ttime1itQR   \t#QR    \ttime1itPO   \t#PO    \ttimeFormH \n");
	fprintf(stderr, "  \t%d \t%2.4e \t%2.4e \t%2.4e \t%2.4e \t%d \t%2.4e \t%d \t%2.4e \n", m, qwtime, litime, nrmtime, qrtime, itqr, potime, itpo, Htime);
    }
    if (myrank_mpi == 0) {
        fprintf(stderr, "#\n");
        fprintf(stderr, "# \t#QR  \t#PO  \n");
	fprintf(stderr, "  \t%d  \t%d \n", itqr, itpo);
    }

    free( tau );
    if ( !optcond ){
        free( Wi );
    }
    if ( Work1 == NULL ) {
        free( U );
        //free( B );
    }
    if ( Work2 == NULL ) {
        //free( A );
        free( B );
    }

    if (verbose & myrank_mpi == 0) { fprintf(stderr, "Exiting QDWH\n");}
    return 0;

}
