/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file pdgezolopd.c
 *
 *  ZOLOPD is a high performance software framework for computing 
 *  the polar decomposition on distributed-memory manycore systems provided by KAUST
 *
 * @version 3.0.0
 * @author Dalal Sukkari
 * @author Hatem Ltaief
 * @date 2018-11-08
 *
 **/

#include "polar.h"

extern void pdgenm2( double *A, int M, int N, int descA[9], double *W, int descW[9], double *Sx, int descSx[9], double *e, double tol);

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
 *  PDGZOLOPD computes the polar decomposition of a real distributed M-by-N  
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
 *          If the symmetric polar factor is not needed, then H will be used as a workspace
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
 *  WORK1   (local workspace/output) DOUBLE PRECISION   array, dimension
 *          (LWORK* NQ)
 *          On exit, if INFO = 0, WORK1(1) returns the optimal LWORK.
 *
 *  LWORK1  (local input) INTEGER
 *          The dimension of the array WORK.
 *          LWORK must be at least LWORK >= MP*NQ
 *          If LWORK = -1, then LWORK is global input and a workspace    
 *          query is assumed; the routine only calculates the minimum
 *          and optimal size for all work arrays. Each of these
 *          values is returned in the first entry of the corresponding
 *          work array, and no error message is issued by PXERBLA.         
 *
 *  WORK2   (local workspace/output) DOUBLE PRECISION   array, dimension
 *          (LWORK2)
 *          On exit, if INFO = 0, WORK2(1) returns the optimal LWORK.
 *
 *  LWORK2  (local input) INTEGER
 *          The dimension of the array WORK.
 *          LWORK must be at least LWORK >= MP*NQ
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


int pdgezolopd( char *jobh, 
          int M, int N,
          double *A_all, int iA_all, int jA_all, int descA_all[9], 
          double *B_all, int iB_all, int jB_all, int descB_all[9],
          double *Work1, int lWork1, 
          double *Work2, int lWork2, 
          int *info)
{

    int    init = 0;
    double eps;
    double tol1;
    double tol3;

    complex dd, sqd, a1;
    double conv = 100.;
    double alpha, beta, Anorm, norm_est, Li;
    double tol = 1.e-1;
    double flops_dgeqrf, flops_dorgqr, flops_dgemm, flops_dpotrf, flops_dtrsm, flops_dtrtri;
    int MB = 2*M;
    int it, itconv, facto = -1;
    int itqr = 0, itpo =0;
    int i0 = 0, i1 =1, iM = M+1;
    int myrank_mpi, nprocs_mpi;

    double qwtime = 0.0, nrmtime = 0.0, litime = 0.0;
    double potime = 0.0, qrtime = 0.0;
    double Htime = 0.0, nstime = 0.0;
    double reduced_qwtime = 0.0, reduced_nrmtime = 0.0, reduced_litime = 0.0; 
    double reduced_potime = 0.0, reduced_qrtime = 0.0; 
    double reduced_Htime = 0.0, reduced_nstime = 0.0;

    int verbose = 0, prof = 0, optcond = 0, symm = 0, ns = 0;
    double *tau_all   = (double *)malloc(N*sizeof(double)) ;
    double *flops;   
    *flops = 0.;

    int mloc, nloc, mlocW;   
    int myrow, mycol, nprow, npcol;   
    int mloc_all, nloc_all, mlocW_all, nb;   
    int ctxt_ = 1, nb_ = 5;
    int ictxt; 
    int wantH;

    int myrow_all, mycol_all, nprow_all, npcol_all;   
    int ictxt_all;

    int k, i, j;
    int lWork3 = -1; 

    double con, kp, K, sn, cn, tn, ff_1, ff_con, enu, den;
    int nbprob;
    double *c, Rnorm, Rinvnorm, condest_R, maxc;
    int m, m_zol, itmax, ii, jj;
    int r;
    int howqr = 1;

    int *Wi_all = (int *)malloc(1*sizeof(int));
    int lWi_all = -1;
    int iinfo;

    /*
     * Get the grid parameters
     */
    if ( verbose ) fprintf(stderr, "Getting the grid parameters \n");
    ictxt_all = descA_all[ctxt_];
    Cblacs_get( -1, 0, &ictxt_all );
    nb = descA_all[nb_];
    Cblacs_gridinfo( ictxt_all, &nprow_all, &npcol_all, &myrow_all, &mycol_all );
    mloc_all  = numroc_( &M, &nb, &myrow_all, &i0, &nprow_all );
    nloc_all  = numroc_( &N, &nb, &mycol_all, &i0, &npcol_all );

    int lmin1, lmin2, lquery;
    *info = 0; 
    lquery =  (lWork1 == -1 || lWork2 == -1); 
    wantH = 0;

    MPI_Comm_rank( MPI_COMM_WORLD, &myrank_mpi );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs_mpi );
   /*
    * Test the input parameters
    */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, "Testing the input parameters \n");
    if ( nprow_all == -1 ){
        *info = -(700+ctxt_);
    }
    else { 
        if ( M < N ){
	   fprintf(stderr, "error(m >= n is required)") ;
	   return -1;
        }
        if ( jobh[0] == 'H' || jobh[0] == 'h' ){
           wantH = 1;
        }

        if ( verbose && myrank_mpi == 0 ) fprintf(stderr, "Testing the input parameters 280  \n");
        int i2 = 2, i3 = 3, i7 = 7, i11 = 11;
        int *idum1, *idum2;
        idum1 = (int *)malloc(2*sizeof(int)) ;
        idum2 = (int *)malloc(2*sizeof(int)) ;
        chk1mat_(&M, &i2, &N, &i3, &iA_all, &jA_all, descA_all, &i7, info);
        //if (wantH){
        chk1mat_(&M, &i2, &N, &i3, &iB_all, &jB_all, descB_all, &i7, info);
        //}
        //igamx2d_(descA[ctxt_], "A", " ", &i1, &i1, info, &i1, &i1, &i1, &i_1, &i_1, &i0);

        lquery =  (lWork1 == -1 || lWork2 == -1); 
        if ( *info == 0 ){
           lmin1 = mloc_all;
           //lmin2 = mloc_all; 
           lmin2 = -1;
           if ( optcond ){
              pdgeqrf_( &nb, &M, 
                        A_all, &i1, &i1, descA_all, 
                        tau_all, 
                        Work2, &lmin2, 
                        &iinfo );
           }
           else {
              pdgecon_( "1", &M, 
                        A_all, &i1, &i1, descA_all, 
                        &Anorm, &Li, 
                        Work2, &lmin2, Wi_all, &lWi_all, 
                        &iinfo );
           }
           lmin2 = (int)Work2[0]; 
           lWi_all = Wi_all[0];
           //Work1[0] = lmin1;
           //Work2[0] = lmin2;
           lquery =  (lWork1 == -1 || lWork2 == -1); 
           if ( (lWork1 < lmin1) & !lquery ){
              *info = -13;
           }
           if ( (lWork2 < lmin2) & !lquery ){
              *info = -15;
           }
        }

        idum1[0] = wantH;
        if ( lWork1 == -1 || lWork2 == -1) {
           idum1[1] = -1;
        }
        else {
           idum1[1] =  1;
        }
        idum2[0] =  1;
        idum2[1] =  15;
        pchk1mat_( &M, &i2, &N, &i3, &iA_all, &jA_all, descA_all, &i7, &i2, idum1, idum2,
                   info );
        //if ((*info == 0) && wantH){
        if (( *info == 0 )){
           pchk1mat_( &M, &i2, &N, &i3, &iB_all, &jB_all, descB_all, &i11, &i0, idum1, idum2,
                      info );
        }
    }

    if ( *info != 0 ){
        pxerbla_( ictxt_all, "PDGEZOLOPD", &(int){-1*info[0]} ); 
        return 0;
    }
    else if ( lquery ){
        lWork1 = -1; lWork2 = -1; lWi_all = -1;
        lWork1 = mloc_all;
        Work1[0] = lWork1;       
        if ( optcond ){
           pdgeqrf_( &nb, &M, 
                     A_all, &i1, &i1, descA_all, 
                     tau_all, 
                     Work2, &lWork2, 
                     &iinfo );
        }
        else {
           pdgecon_( "1", &M, 
                     A_all, &i1, &i1, descA_all, 
                     &Anorm, &Li, 
                     Work2, &lWork2, Wi_all, &lWi_all, 
                     &iinfo );
        }
        lWork2 = (int)Work2[0]; 
        Work2[0] = lWork2;
        return 0;
        if ( verbose && myrank_mpi == 0 ) fprintf(stderr, "Testing the input parameters 368  \n");
    } 

    /* Quick return if possible */
    if ( M == 0 || N == 0 ){
        return 0;
    }

    /**
     * Needed for debugging the code
    if ( Work1 == NULL ) {
	Work1  = (double *)malloc(mloc_all*nloc_all*sizeof(double));
    }
    if ( Work2 == NULL ) {
        lWork2 = -1;
        if(optcond){
           pdgeqrf_( &nb, &M, 
                    A_all, &i1, &i1, descA_all, 
                    tau_all, 
                    Work2, &lWork2, 
                    &iinfo );
        }
        else {
           pdgecon_( "1", &M, 
                     A_all, &i1, &i1, descA_all, 
                     &Anorm, &Li, 
                     Work2, &lWork2, NULL, &i0, 
                     &iinfo );
        }
        lWork2= (int)Work2[0];
	Work2  = (double *)malloc(lWork2*sizeof(double));
    }
    */

    /*
     * Use Work1 for Acpy_all
     */
    int descAcpy_all[9];
    double *Acpy_all = Work1;
    descinit_( descAcpy_all, &M, &N, &nb, &nb, &i0, &i0, &ictxt_all, &mloc_all, &iinfo );


    /*
     * Save copy of A_all in Acpy_all ==> for other computations and for H = U'*A
     */
    pdlacpy_( "All", &M, &N, 
              A_all,    &i1, &i1, descA_all, 
              Acpy_all, &i1, &i1, descAcpy_all ); 

    if ( prof ) {qwtime =- MPI_Wtime();}
   /*
    * Estimate the condition number and find the number of subproblems using all mpi processes
    * Li: reciprocal of condition number estimation
    */
    if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Cond estimate starts\n");}
    if ( prof ) {litime =- MPI_Wtime();}

    int lWork21 = -1;
    if ( !optcond ) {
        pdgecon_( "1", &M, 
                  A_all, &i1, &i1, descA_all, 
                  &Anorm, &Li, 
                  Work2, &lWork21, Wi_all, &lWi_all, 
                  &iinfo );
        lWork21 = (int)Work2[0]; lWi_all = Wi_all[0];
        Wi_all  = (int *)malloc(N*sizeof(int)) ;
    }
    //Work2 = (double *)malloc(lWork2*sizeof(double));


    if ( optcond ){
        // estimate condition number using QR 
        Anorm = pdlange_( "1", &M, &N, 
                          A_all, &i1, &i1, descA_all, 
                          Work2 );
        pdgeqrf_( &M, &N, 
                  A_all, &i1, &i1, descA_all, 
                  tau_all, Work2, &lWork2, 
                  &iinfo );
        pdtrtri_( "U", "N", &N, 
                  A_all, &i1, &i1, descA_all, 
                  &iinfo );
        double Ainvnorm = pdlange_( "1", &M, &N, 
                                    A_all, &i1, &i1, descA_all, 
                                    Work2 );
        Li = Ainvnorm*Anorm;    
        Li = Anorm/Li;    
        //Li = Li/sqrt(N); // This reduces the accuracy    
        *flops += FLOPS_DGEQRF( M, N )
               + FLOPS_DTRTRI(  N );
    }
    else{
       // estimate condition number using LU 
       Anorm = pdlange_( "1", &M, &N, 
                         A_all, &i1, &i1, descA_all, 
                         Work2 );
       pdgetrf_( &M, &N, 
                 A_all, &i1, &i1, descA_all, 
                 Wi_all, 
                 &iinfo );
       pdgecon_( "1", &M, 
                 A_all, &i1, &i1, descA_all, 
                 &Anorm, &Li, 
                 Work2, &lWork2, Wi_all, &lWi_all, 
                 &iinfo );
       Li = Anorm*Li;    
       //Li = Li/sqrt(N); // This reduces the accuracy    
       *flops += FLOPS_DGETRF(M, N);
       //*flops += FLOPS_DGETRF(M, N) + 2. * FLOPS_DTRSM( 'L', N, 1 );
    }
    if( prof ) {litime += MPI_Wtime();
        MPI_Allreduce( &litime, &reduced_litime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    }

   /*
    * Computing the second norm estimate
    * Scale the matrix by norm_est
    */
    if ( verbose && myrank_mpi == 0 ) {fprintf(stderr, " Computing the second norm estimate\n");}
    if ( prof ) {nrmtime =- MPI_Wtime();}
    alpha = 1.0; 
    pdgenm2( Acpy_all, M, N, descAcpy_all, A_all, descA_all, B_all, descB_all, &norm_est, tol );
    if ( prof ) {nrmtime += MPI_Wtime();
        MPI_Allreduce( &nrmtime, &reduced_nrmtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    }

    pdlacpy_( "All", &M, &N, 
              Acpy_all, &i1, &i1, descAcpy_all, 
              A_all,    &i1, &i1, descA_all ); 
    pdlascl_( "G", &norm_est, &alpha, &M, &N, 
              A_all, &i1, &i1, descA_all, 
              &iinfo );

    /*
     * scale the matrix by Li
     */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, " Scale by Li\n");
    alpha = 1.0; 
    pdlascl_( "G", &Li, &alpha, &M, &N, 
              A_all, &i1, &i1, descA_all, 
              &iinfo );

    con = 1.0/Li;
    choosem( con, &m_zol );
    nbprob = m_zol;

   /*
    * Computing the number of processor per subproblem 
    */
    if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, " Computing the number of processor per subproblem\n");}
    //int nbproc = nprow*npcol;
    //nprow = 2; npcol =2;
    int nbproc = (int)(nprow_all * npcol_all) / m_zol;


    /* Check if the nbproc is prime */
    int flag = 0;
    for( i = 2; i <= nbproc/2; ++i )
    {
        // condition for nonprime number
        if(nbproc%i == 0) { flag = 1; break;}
    }
    {
        if (flag == 0)
          /* nbproc is a prime number */
          { nbproc = nbproc - 1;}
    }

    int max_rank = nbproc * m_zol ;
    //nprow = ceil(nprow_all/2); 
    //nprow = 0;
    //npcol = 0;
    //npcol = ceil(npcol_all/4);
  
    int ndim[2] = {0, 0};
    MPI_Dims_create( nbproc, 2, ndim );
    nprow = (myrank_mpi < max_rank)? min(ndim[0], ndim[1]) : 0;
    npcol = (myrank_mpi < max_rank)? max(ndim[0], ndim[1]) : 0;
    nprow = min(ndim[0], ndim[1]);
    npcol = max(ndim[0], ndim[1]);
    
    if ( myrank_mpi == 0 ){ 
        fprintf(stderr, " The number of subproblems to be solved independently is %d\n", m_zol);

        if (flag == 0){                            
           fprintf(stderr, " The number of processors per subproblem is %d, for a better grid configuration we will use %d with grid configuration %d x %d\n", nbproc + 1, nbproc, nprow, npcol);}
        else if (flag == 1){ 
           fprintf(stderr, " The number of processors per subproblem is %d, with grid configuration %d x %d\n", nbproc, nprow, npcol);}

        fprintf(stderr, " There will be %d unused processors\n", (nprow_all * npcol_all - nbproc * m_zol));
        fprintf(stderr, " To use all the processors, use number of processors = nonprime x %d \n", m_zol);
    }

    int color = (myrank_mpi < max_rank)? 0 : 1;
    MPI_Comm  new_comm; 
    MPI_Comm_split(MPI_COMM_WORLD, color, myrank_mpi, &new_comm);

    if ( con < 2 ) 
        {itmax = 1;}
    else 
        {itmax = 2;} 
    //else 
    //    {itmax = 3;} // need this for ill-cond matrix with cond=1e16

    /*
     * Map processes to different context 
     */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, " Map processes to different context\n");
    int *imap = (int *)malloc(nprow*npcol*sizeof(int));
    int *ictxt_id = (int *)malloc(nbprob*sizeof(int));
    memset(ictxt_id, -1, nbprob*sizeof(int));

    k = 0;
    //printf("me %d nprow %d npcol %d \n", myrank_mpi, nprow, npcol);
    //nprow npcol: per group
    for (i = 0; i < nprow; i++){
        for (j = 0; j < npcol; j++){
           *(imap + i + j * nprow) = nprow*npcol*(int)(myrank_mpi/(nprow*npcol)) + k;
           k = k + 1;
        }
        //printf("\n");
    }

    /**
     * Create the required workspaces for independent problem solving
     * Initial the ctxt per independent problem
     */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, "Desc Init and ctxt create of subproblems\n");
    int descU[9], descA[9], descB[9], descUcpy[9];
    double *U=NULL,  *A=NULL, *B =NULL, *Ucpy=NULL, *tau;

if ( myrank_mpi < max_rank ) 
{
 
    int bhandle = Csys2blacs_handle(new_comm);
    ictxt = bhandle;

    //Cblacs_get( 0, 0, &ictxt );
    Cblacs_gridmap( &ictxt, imap, nprow, nprow, npcol );
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    *(ictxt_id + (int)(myrank_mpi/(nprow*npcol))) = ictxt;

    mloc  = numroc_( &M, &nb, &myrow, &i0, &nprow );
    nloc  = numroc_( &N, &nb, &mycol, &i0, &npcol );
    mlocW = numroc_( &MB, &nb, &myrow, &i0, &nprow );

    descinit_( descU, &M, &N, &nb, &nb, &i0, &i0, &ictxt, &mloc, &iinfo );
    descinit_( descA, &M, &N, &nb, &nb, &i0, &i0, &ictxt, &mloc, &iinfo );
    descinit_( descUcpy, &M, &N, &nb, &nb, &i0, &i0, &ictxt, &mloc, &iinfo );
    descinit_( descB, &MB, &N, &nb, &nb, &i0, &i0, &ictxt, &mlocW, &iinfo );

    descU[1] = (myrank_mpi < max_rank)? descU[1]:-1;
    descA[1] = (myrank_mpi < max_rank)? descA[1]:-1;
    descUcpy[1] = (myrank_mpi < max_rank)? descUcpy[1]:-1;
    descB[1] =  (myrank_mpi < max_rank)? descB[1]:-1;

    assert( mloc >= 0);
    U     = (double *)malloc(mloc*nloc*sizeof(double)) ;
    A     = (double *)malloc(mloc*nloc*sizeof(double)) ;
    Ucpy  = (double *)malloc(mloc*nloc*sizeof(double)) ;
    B     = (double *)malloc(mlocW*nloc*sizeof(double)) ;
    tau   = (double *)malloc(nloc*sizeof(double)) ;

    /**
     * Copy the input matrix to U in the different ctxt
     */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, "Start copy the input matrix to U in the different ctxt \n" );
    //descU[1] = ictxt;
    int group_id = (int)(myrank_mpi/(nprow*npcol));
    for ( j = 0; j < nbprob; j++ ){
        int temp = descU[1];
        if ( group_id != j ){
           descU[1] = -1;
        }
        pdgemr2d_( &M, &N, 
                   A_all, &i1, &i1, descA_all, 
                   U,     &i1, &i1, descU, 
                   &ictxt_all );
        if ( group_id != j ){
           descU[1] = temp;
        }
    }

    double *Work = (double *)malloc(1*sizeof(double));
    lWork3 = -1;

    pdgeqrf_( &MB, &N, B, &i1, &i1, descB, tau, Work, &lWork3, &iinfo );
    lWork3 = (int)Work[0];
    Work  = (double *)malloc(lWork3*sizeof(double)) ;
 
    if ( !init ) {
	eps  = pdlamch_( &ictxt, "E" ); 
	tol1 = 5. * eps;
	tol3 = pow(tol1, 1./3.);
	init = 1;
    }

    if ( M < N ){
	fprintf(stderr, "error(m >= n is required)") ;
	return -1;
    }

    con = 1.0/Li;
    it = 0;
    choosem(con, &m);
    m_zol = m; 
    if ( con < 2 ) 
        {itmax = 1;}
    else 
        {itmax = 2;} 
    //else 
    //    {itmax = 3;} // need this for ill-cond matrix with cond=1e16

    /*
     * Allocate U_ac on group0 so that all groups copy to it
     */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, " Allocate U_ac on group0\n" );
    group_id= (int)(myrank_mpi/(nprow*npcol));
    int descU_ac[9]; double *U_ac; descU_ac[1] = -1;
    if ( group_id == 0 && myrank_mpi < max_rank ){
        descinit_( descU_ac, &M, &N, &nb, &nb, &i0, &i0, &ictxt, &mloc, &iinfo );
        U_ac = (double *)malloc(mloc*nloc*sizeof(double)) ;
    }

    double time_gather, reduced_time_gather;
    double time_bcast, reduced_time_bcast;
    m_zol = nbprob;
    c  = (double *)malloc(2*m_zol*sizeof(double)) ;
    /*
     * Start the iterations to compute the orthogonal polar factor
     */
    if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Start computing the orthogonal polar factor on all the ctxt\n");}

    while ( myrank_mpi < max_rank && it < itmax ){
        it = it+1;
    
        kp = 1/con;
        alpha = acos(kp);
        /**
         * modified elliptic functions
         */
        mellipke( alpha, &K, &beta); 
        memset(c, 0, (2*m_zol)* sizeof(double));
        for ( ii = 1; ii <= (2*m_zol); ii++ ){
           mellipj(ii*K/(2*m_zol+1),alpha, &sn, &cn, &tn, Work);
           c[ii-1] = (sn*sn)/(cn*cn);
        }
    
        /*
         * U = computeAA(U,c,it);
         * the iterations (QR, Chol) on the matrix
         */
        pdlacpy_( "A", &M, &N, 
                  U,    &i1, &i1, descU, 
                  Ucpy, &i1, &i1, descUcpy );

        r = m_zol;//nbprob;
        ii = (int)(myrank_mpi/(nprow*npcol)) + 1;
        //r = m_zol;
        //for ( ii = 1; ii <= r; ii++ ){
        enu = 1;
        for ( jj = 1; jj <= r; jj++ ){
           enu = enu*(c[2*ii-2]-c[2*jj-1]);
        }
        den = 1;
        for ( jj = 1; jj <= r; jj++ ){
           if ( ii != jj ){
              den = den*(c[2*ii-2]-c[2*jj-2]);
           }
        }
    
        /*
         * max(c(1:end-1))
         */
        maxc = c[0]; int i;
        for ( i = 0; i < 2*m_zol-1; i++ ){
           if ( c[i] >= maxc )
              maxc = c[i];
        } 
        if ( it <= itmax && maxc > 1e2 ){    
           if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Start the ZOLOPD iterations\n");}
           //qrtime = 0.0;
           if ( prof ) {qrtime =- MPI_Wtime();}
           /*
            * QR-based
            */
           if ( howqr ){ 

	      /**
               * [AAA,~]=qr([A;sqrt(c(2*ii-1))*eye(n)],0);   
               * AA=AA-enu/den/sqrt(c(2*ii-1))*AAA(1:m,:)*AAA(m+1:end,:)' ;       
	       * Generate the matrix B = [ B1 ] = [ U                   ]
	       *                         [ B2 ] = [ sqrt(c(2*ii-2)) *Id ]
	       */
               if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Start the QR based first iteration\n");}
               pdlacpy_( "A", &M, &N, 
                         Ucpy, &i1, &i1, descUcpy, 
                         B,    &i1, &i1, descB );
               alpha = 0.; beta = sqrt(c[2*ii-2]);
               pdlaset_( "G", &M, &N, &alpha, &beta, 
                         B, &iM, &i1, descB );

	       /**
	        * Factorize B = QR, and generate the associated Q
	        */
               pdgeqrf_( &MB, &N, 
                         B, &i1, &i1, descB, 
                         tau, Work, &lWork3, 
                         &iinfo );
               pdorgqr_( &MB, &N, &N, 
                         B, &i1, &i1, descB, 
                         tau, Work, &lWork3, 
                         &iinfo );

	      /**
     	       * Gemm to find the conv-norm
	       *  U = ( (a-b/c)/sqrt(c) ) * Q1 * Q2' + (b/c) * U
               * alpha = -1.*enu/den/sqrt(c[2*ii-2]); beta = 1.;
	       */
              alpha = -1.*enu/den/sqrt(c[2*ii-2]); beta = 0.;
              if ( ii == 1 ) {beta = 1.;}
              pdgemm_( "N", "T", &M, &N, &N, 
                       &alpha, B, &i1, &i1, descB, 
                               B, &iM, &i1, descB, 
                       &beta,  U, &i1, &i1, descU );

	      /* Main flops used in this step */
	      flops_dgeqrf = FLOPS_DGEQRF( MB, N );
	      flops_dorgqr = FLOPS_DORGQR( MB, N, N );
	      flops_dgemm  = FLOPS_DGEMM( M, N, N );
	      *flops += flops_dgeqrf + flops_dorgqr + flops_dgemm;

              itqr += 1;
	      facto = 0;
           }//end of QR-based iteration
           /*
            * Chol-based first iteration
            * Still not tested yet
            */
           else {
	      /*
	       * Compute C = U' * U + c[2*ii-1] * I
               * This part is not tested yet
	       */
              if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Start the Chol based first iteration\n");}

              alpha = 0.; beta = 1.; 
              alpha = 1.; beta = c[2*ii-2];
              //pdgemm_( "T", "N", &M, &N, &N, &alpha, Ucpy, &i1, &i1, descUcpy, Ucpy, &i1, &i1, 
              //         descUcpy, &beta, C, &i1, &i1, descC);

	      /**
	       * dpotrf on C
	       */
              alpha = 1.0; beta = 0.0;
              //pdposv_( "U", &M, &N, C, &i1, &i1, descC, B, &i1, &i1, descB, &iinfo);
              //pdpotrf_( "U", &M, C, &i1, &i1, descC, &iinfo);
              /*
               * condest(R)
               */
              alpha = 0.; 
              pdlaset_( "G", &M, &N, 
                        &alpha, &alpha, 
                        B, &i1, &i1, descB );
              //pdlacpy_( "U", &M, &N, C, &i1, &i1, descC, B, &i1, &i1, descB );
              Rnorm    = pdlange_( "1", &M, &N, 
                                   B, &i1, &i1, descB, 
                                   Work );
              pdtrtri_( "U", "N", &N, 
                        B, &i1, &i1, descB, 
                        &iinfo );
              Rinvnorm = pdlange_( "1", &M, &N, 
                                   B, &i1, &i1, descB, 
                                   Work );
              condest_R = ( 1.0 / Rinvnorm)/Rnorm; 

              /* the following two dtrsm can be done in parallel */

              /*
               * solve Q = A/R; 
               */
              pdlacpy_( "A", &M, &N, 
                        Ucpy, &i1, &i1, descUcpy, 
                        B,    &i1, &i1, descB );
              //pdtrsm_ ("R", "U", "N", "N", &M, &N, &alpha, C, &i1, &i1, descC, B, &i1, &i1, descB);
 
              /*
               * solve II = II/R;  
               */
              alpha = 0.; beta = sqrt(c[2*ii-2]); 
              pdlaset_( "G", &M, &N, &alpha, &beta, 
                        B, &iM, &i1, descB );
              //pdtrsm_ ("R", "U", "N", "N", &M, &N, &alpha, C, &i1, &i1, descC, B, &iM, &i1, descB);

              if ( condest_R > 10 ) {  // repeat
                 /* the following two dgemm can be done in parallel */
                
                 /*
                  * II'*II
                  */ 
                 alpha = 1.; beta = 0.;
                 //pdgemm_( "T", "N", &M, &N, &N, &alpha, B, &iM, &i1, descB, B, &iM, &i1, 
                 //         descB, &beta, C, &i1, &i1, descC);
                 /*
                  * R = Q'*Q + II'*II
                  */ 
                 //pdgemm_( "T", "N", &M, &N, &N, &alpha, B, &i1, &i1, descB, B, &i1, &i1, 
                 //         descB, &alpha, C, &i1, &i1, descC);

	         /**
	          * dpotrf on R
	          */
                 //pdpotrf_( "U", &M, C, &i1, &i1, descC, &iinfo);

                 /* the following two dtrsm can be done in parallel */

                 /*
                  * solve Q = Q/R; 
                  */
                 //pdtrsm_ ("L", "U", "N", "N", &M, &N, &alpha, C, &i1, &i1, descC, B, &i1, &i1, descB);

                 /*
                  * solve II = II/R; 
                  */
                 //pdtrsm_ ("L", "U", "N", "N", &M, &N, &alpha, C, &i1, &i1, descC, B, &iM, &i1, descB);

	         /* Main flops used in this step */
	         flops_dgemm  = 2. * FLOPS_DGEMM( M, N, N );
	         flops_dpotrf = FLOPS_DPOTRF( M );
	         flops_dtrsm  = 2. * FLOPS_DTRSM( 'L', M, N );
	         *flops += flops_dgemm + flops_dpotrf + flops_dtrsm;
              }    
              //AA = AA-enu/den/sqrt(c(2*ii-1))*(Q*II');
              alpha = -enu/den/sqrt(c[2*ii-2]); beta = 1.;
              pdgemm_( "N", "T", &M, &N, &N, 
                       &alpha, B, &i1, &i1, descB, 
                               B, &iM, &i1, descB, 
                       &beta,  U, &i1, &i1, descU );

	      /* Main flops used in this step */
	      flops_dgemm  = 2. * FLOPS_DGEMM( M, N, N );
	      flops_dpotrf = FLOPS_DPOTRF( M );
              flops_dtrtri = FLOPS_DTRTRI(  N );
	      flops_dtrsm  = 2. * FLOPS_DTRSM( 'R', M, N );
	      *flops += flops_dgemm + flops_dpotrf + flops_dtrtri + flops_dtrsm;

           }//end Chol-based first iteration 
           if ( prof ) {qrtime += MPI_Wtime();
              MPI_Allreduce( &qrtime, &reduced_qrtime, 1, MPI_DOUBLE, MPI_MAX, new_comm );
              if ( myrank_mpi == 0 ) {
	         fprintf(stderr, "#  \treduced_qrtime    \n");
	         fprintf(stderr, "  \t%2.4e \n", reduced_qrtime);
              }
           }
        }// end of (it<=1 && maxc>1e2)
        else {// Chol-based second iteration
	   /**
	    * Compute C = U * U' + c[2*ii-1] * I
	    */
           if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Start the Chol based second iteration\n");}
           //potime = 0.0;
           if( prof ) {potime =- MPI_Wtime();}

           alpha = 0.; beta =1.; 
           pdlaset_( "G", &M, &N, &alpha, &beta, 
                     B, &i1, &i1, descB);
           alpha = 1.; beta = c[2*ii-2];
           pdgemm_( "T", "N", &M, &N, &N, 
                    &alpha, Ucpy, &i1, &i1, descUcpy, 
                            Ucpy, &i1, &i1, descUcpy, 
                    &beta,  B,    &i1, &i1, descB );
	   /**
	    * dpotrf on R
	    */
           pdpotrf_( "U", &M, 
                     B, &i1, &i1, descB, 
                     &iinfo );

           /*
            * solve Qtmp = A/Cinv; 
            */
           alpha = 1.; beta = 0.;
           pdlacpy_( "A", &M, &N, 
                     Ucpy, &i1, &i1, descUcpy, 
                     A,    &i1, &i1, descA );
           pdtrsm_( "R", "U", "N", "N", &M, &N, 
                    &alpha, B, &i1, &i1, descB, 
                            A, &i1, &i1, descA );
           /*
            * solve Qtmp=Qtmp/Cinv'; 
            */
           pdtrsm_( "R", "U", "T", "N", &M, &N, 
                    &alpha, B, &i1, &i1, descB, 
                            A, &i1, &i1, descA );
           alpha = -enu/den; beta = 0.;
           if ( ii == 1 ){
              beta = 1.;}
           pdgeadd_( "N", &M, &N, 
                     &alpha, A, &i1, &i1, descA, 
                     &beta,  U, &i1, &i1, descU );
	   /* Main flops used in this step */
	   flops_dgemm  = FLOPS_DGEMM( M, N, N );
	   flops_dpotrf = FLOPS_DPOTRF( M );
	   flops_dtrsm  = 2. * FLOPS_DTRSM( 'R', M, N );
	   *flops += flops_dgemm + flops_dpotrf + flops_dtrsm;

           if ( prof ) {potime += MPI_Wtime();
              MPI_Allreduce( &potime, &reduced_potime, 1, MPI_DOUBLE, MPI_MAX, new_comm );
              if ( myrank_mpi == 0 ) {
	                fprintf(stderr, "#  \treduced_potime    \n");
	                fprintf(stderr, "  \t%2.4e \n", reduced_potime);
              }
           }

           itpo += 1;
	   facto = 1;
        }//end Chol-based second iteration
        //}//end for-loop over r 
        /*
         * end of computeAA
         */
    
        /*
         * scale the matrix by ff(1)
         */
        ff_1 = 1.; ff_con = con;
        for ( i=1; i <= m; i++ ){
           ff_1 = ff_1 * (1.+c[2*i-1])/(1.+c[2*i-2]); 
           ff_con = ff_con * (con*con+c[2*i-1])/(con*con+c[2*i-2]); 
        }
        if( myrank_mpi < nprocs_mpi - max_rank ) 
           MPI_Send(&ff_1, 1, MPI_DOUBLE, max_rank + myrank_mpi, 0, MPI_COMM_WORLD);
        
        //if ( con < 2 ) { con = max(ff_con/ff_1,1);break;} // if k(A) small, one step is enough
        con = max(ff_con/ff_1,1);
        //if symm;    U=(U'+U)/2;end    // force symmetry if A symmetric

        /*
         * Joint here to accumulate the matrices form different mysim  
         * All the different contexts can copy the resulting matrix into a matrix with the global context  
         * U = U_s(mysim=0) + U_s(mysim=1) + ... + U_s(mysim=r-1) 
         * Note: we can copy a matrix with the global contxt to partial contxt, but we can't do the otherway!!
         */
         
        /*
         * Gather + Accumulation = reduction
         */
        MPI_Barrier( new_comm );
        if ( prof ){ time_gather =- MPI_Wtime();}
        if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Gather the result matrix from different ctxt\n");}
        alpha = 1.; 
        int j;

        if( it < itmax ) {         
           if ( group_id == 0 )
              descU[1] = -1;

           for ( j = 1; j < nbprob; j++ ){
              if ( group_id != 0 && group_id != j ){
                 descU[1] = -1;
              }

	      MPI_Barrier( new_comm );
              pdgemr2d_( &M, &N, 
                         U,    &i1, &i1, descU, 
                         U_ac, &i1, &i1, descU_ac, 
                         &ictxt_all );
	      MPI_Barrier( new_comm );

              if (group_id != 0 && group_id != j ){
                 descU[1] = ictxt;
              }

              if ( group_id == 0 ) {
                 descU[1] = ictxt;
                 pdgeadd_( "N", &M, &N, 
                           &alpha, U_ac, &i1, &i1, descU_ac, 
                           &alpha, U,    &i1, &i1, descU );
                 descU[1] = -1;
              }
           }

           /**
            * Scale/Symmetrize the matrix on group id 0
            * U = U/ff(1);
            */
           if ( group_id == 0 ) {
              descU[1] = ictxt;
              alpha = 1.0; 
              /**
               * normalize so that min(svd(A))=1
               * U = U/ff(1);
               */
              pdlascl_( "G", &ff_1, &alpha, &M, &N, 
                        U, &i1, &i1, descU, 
                        &iinfo ); 
              /**
               * if symm;    U=(U'+U)/2;end    // force symmetry if A symmetric
               */
              if ( symm ){
                 pdlacpy_( "A", &M, &N, 
                           U, &i1, &i1, descU, 
                           B, &i1, &i1, descB );
                 alpha = 0.5; 
                 pdgeadd_( "T", &M, &N, 
                           &alpha, B, &i1, &i1, descB, 
                           &alpha, U, &i1, &i1, descU );
              }
           }
        }
        else {
           alpha = 1.; 
           for ( j = 0; j < nbprob; j++ ){
              if ( group_id != j ){
                 descU[1] = -1;
              }

              MPI_Barrier( new_comm );
              pdgemr2d_( &M, &N, 
                         U,     &i1, &i1, descU, 
                         B_all, &i1, &i1, descB_all, 
                         &ictxt_all );
              beta = 1.; if (j ==0 ) {beta = 0.;}
              pdgeadd_( "N", &M, &N, 
                        &alpha, B_all, &i1, &i1, descB_all, 
                        &beta,  A_all, &i1, &i1, descA_all );

              if ( group_id != j ){
                 descU[1] = ictxt;
              }
           }
           alpha = 1.0; 
           pdlascl_( "G", &ff_1, &alpha, &M, &N, 
                     A_all, &i1, &i1, descA_all, 
                     &iinfo ); //U = U/ff(1);// normalize so that min(svd(A))=1
           if ( symm ){
              pdlacpy_( "A", &M, &N, 
                        A_all, &i1, &i1, descA_all, 
                        B_all, &i1, &i1, descB_all );
              alpha = 0.5; 
              pdgeadd_( "T", &M, &N, 
                        &alpha, B_all, &i1, &i1, descB_all, 
                        &alpha, A_all, &i1, &i1, descA_all );
           }
        }
        MPI_Barrier( new_comm );
        if ( prof ){ time_gather += MPI_Wtime();
           MPI_Allreduce( &time_gather, &reduced_time_gather, 1, MPI_DOUBLE, MPI_MAX, new_comm );
        }

        /*
         * Bcast 
         */
        MPI_Barrier( new_comm );
        if ( prof ){ time_bcast =- MPI_Wtime();}
        if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Bcast the result matrix to all the ctxt's\n");}

        if ( it < itmax ){ 
         
           if ( group_id == 0 ){
              pdlacpy_( "A", &M, &N, 
                        U,    &i1, &i1, descU, 
                        U_ac, &i1, &i1, descU_ac );
              descU[1] = -1;
           }

           descU[1] = ictxt;
           for ( j = 1; j < nbprob; j++ ){
              if ( group_id != j ){
                 descU[1] = -1;
              }
              MPI_Barrier( new_comm );
              pdgemr2d_( &M, &N, 
                         U_ac, &i1, &i1, descU_ac, 
                         U,    &i1, &i1, descU, 
                         &ictxt_all );
              MPI_Barrier( new_comm );
              if ( group_id != j ){
                 descU[1] = ictxt;
              }
           }
           descU[1] = ictxt;
           MPI_Barrier( new_comm );
           if ( prof ){ time_bcast += MPI_Wtime();
              MPI_Allreduce( &time_bcast, &reduced_time_bcast, 1, MPI_DOUBLE, MPI_MAX, new_comm );
           }
        }
         
        if ( prof && myrank_mpi == 0 ) {
           fprintf(stderr, "# \treduced_time_gather \treduced_time_bcast     \n");
	   fprintf(stderr, "  \t%2.4e \t\t%2.4e \n", reduced_time_gather, reduced_time_bcast);
        }
    } // end of while-loop
    if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "Done computing the orthogonal polar factor\n");}

    if ( ns ){
        assert(0);
        // U = (3/2)*U-U*(U'*U)/2; % Newton-Schulz refinement.
        if ( verbose && myrank_mpi == 0 ) { fprintf(stderr, "# \tNewton-Schulz refinement \n");}
        //nstime = 0.0;
        if ( prof ) {nstime =- MPI_Wtime();}
        alpha = 1./2.; beta = 0.;
        pdgemm_( "T", "N", &M, &N, &N, 
                 &alpha, A_all, &i1, &i1, descA_all, 
                         A_all, &i1, &i1, descA_all, 
                 &beta,  B_all, &i1, &i1, descB_all );
        alpha = -1.; beta = 3./2.;
        /* Save a copy of Acpy_all in U */
        if ( group_id != 1 ){
           descB[1] = -1;
        }
        pdgemr2d_( &M, &N, 
                   Acpy_all, &i1, &i1, descAcpy_all, 
                   B,        &i1, &i1, descB, 
                   &ictxt_all );
        if ( group_id != j){
           descB[1] = ictxt;
        }
        pdlacpy_( "A", &M, &N, 
                  A_all,    &i1, &i1, descA_all, 
                  Acpy_all, &i1, &i1, descAcpy_all );
        pdgemm_( "N", "N", &M, &N, &N, 
                 &alpha, Acpy_all, &i1, &i1, descAcpy_all, 
                         B_all,    &i1, &i1, descB_all, 
                 &beta,  A_all,    &i1, &i1, descA_all );
        /* Copy of  U to Acpy_all*/
        if ( group_id != 1 ){
           descB[1] = -1;
        }
        pdgemr2d_( &M, &N, 
                   B,        &i1, &i1, descB, 
                   Acpy_all, &i1, &i1, descAcpy_all, 
                   &ictxt_all );
        if ( group_id != j ){
           descB[1] = ictxt;
        }
        if ( prof ) {nstime += MPI_Wtime();
           MPI_Allreduce( &nstime, &reduced_nstime, 1, MPI_DOUBLE, MPI_MAX, new_comm );
           if ( myrank_mpi == 0 ) {
	             fprintf(stderr, "# \tNewton-Schulz-time     \n");
	             fprintf(stderr, "  \t%2.4e\n", reduced_nstime);
           }
        }
        /* Main flops used in this step: Newton-Schulz refinement*/
        flops_dgemm  = 2. * FLOPS_DGEMM( M, N, N );
        *flops += flops_dgemm;
    }

    /*
     * A = U*H ==> H = U'*A ==> H = 0.5*(H'+H)
     */
    if ( verbose && myrank_mpi == 0 ) fprintf(stderr, "Computing the symmetric positive semidefinite matrix H \n");
    if ( wantH ){ 
        if ( prof ) {Htime =- MPI_Wtime();}
        alpha = 1.0; beta = 0.0;
        pdgemm_( "T", "N", &M, &N, &N, 
                 &alpha, A_all,    &i1, &i1, descA_all, 
                         Acpy_all, &i1, &i1, descAcpy_all, 
                 &beta,  B_all, &i1, &i1, descB_all );
        pdlacpy_( "A", &M, &N, 
                  B_all,    &i1, &i1, descB_all, 
                  Acpy_all, &i1, &i1, descAcpy_all );
        alpha = 0.5; 
        pdgeadd_( "T", &M, &N,  
                  &alpha, Acpy_all, &i1, &i1, descAcpy_all, 
                  &alpha, B_all,    &i1, &i1, descB_all );
    
        /* Main flops used in this step: Calculating H */
        flops_dgemm  = FLOPS_DGEMM( M, N, N );
        *flops += flops_dgemm;
        if ( prof ) {Htime += MPI_Wtime();
           MPI_Allreduce( &Htime, &reduced_Htime, 1, MPI_DOUBLE, MPI_MAX, new_comm );
           if ( myrank_mpi == 0 ) {
	          fprintf(stderr, "# \treduced_Htime     \n");
	          fprintf(stderr, "  \t%2.4e\n", reduced_Htime);
           }
        }
    }

    if ( prof ) {qwtime += MPI_Wtime();
        MPI_Allreduce( &qwtime, &reduced_qwtime, 1, MPI_DOUBLE, MPI_MAX, new_comm );
        if ( myrank_mpi == 0 ) { printf("\n reduced_qwtime_end %2.4e \n", reduced_qwtime);}
    }

    if ( prof && myrank_mpi == 0  ) {
        fprintf(stderr, "# ZOLOPD Profiling \n"); 
        fprintf(stderr, "#\n");
        fprintf(stderr, "# \tN    \ttimeZOLOPD     \ttimeLi   \ttime2nr     #itmax  \ttime1itQR   \t#QR  \treduced_time_gather  \treduced_time_bcast  \ttime1itPO   \t#PO    \ttimeFormH \n");
	fprintf(stderr, "  \t%d \t%2.4e \t%2.4e \t%2.4e  \t%d  \t%2.4e  \t%d \t%2.4e  \t\t%2.4e  \t\t%2.4e \t%d \t%2.4e \n", M, reduced_qwtime, reduced_litime, reduced_nrmtime, itmax, reduced_qrtime, itqr*nbprob, reduced_time_gather, reduced_time_bcast, reduced_potime, itpo*nbprob, reduced_Htime);
    }
    if ( myrank_mpi == 0 ) {
        fprintf(stderr, "#\n");
        fprintf(stderr, "# #itmax  \t#QR  \t#PO  \n");
	fprintf(stderr, "  \t%d  \t%d  \t%d \n", itmax, nbprob*itqr, nbprob*itpo);
    }

    free( U ); free( A ); free( B ); free( Ucpy ); free( tau ); 
    free( imap ); free( ictxt_id );
    free( c );
    free( Work );
    if ( group_id == 0)
         free( U_ac );
    if ( Work1 == NULL ) 
         free( Work1 );
    if(Work2 == NULL) 
         free( Work2 );

} else {

    int group_id= (int)(myrank_mpi/(nprow*npcol));
    int bhandle = Csys2blacs_handle(new_comm);
    double* Work11 = (double *)malloc(1*sizeof(double));

    int descU_ac[9]; double *U_ac; descU_ac[1] = -1;
    descU_ac[1] = -1;
    kp = 1/con;
    alpha = acos(kp);
    lWork3 = -1;
    descU[1] = -1;
        for ( j = 0; j < nbprob ; j++ ){
           pdgemr2d_( &M, &N, 
                      A_all, &i1, &i1, descA_all, 
                      U,     &i1, &i1, descU, 
                      &ictxt_all );
        }
        int it1 = 0;

    while ( it1 < itmax ){
    
        it1 = it1 + 1;

        MPI_Status stat;
        MPI_Recv( &ff_1, 1, MPI_DOUBLE, myrank_mpi - max_rank, 0,
                  MPI_COMM_WORLD, &stat );

        if ( it1 < itmax ) {         
           for ( j = 1; j < nbprob; j++ ){
              MPI_Barrier( new_comm );
              pdgemr2d_( &M, &N, 
                         U,    &i1, &i1, descU, 
                         U_ac, &i1, &i1, descU_ac, 
                         &ictxt_all );
              MPI_Barrier( new_comm );

              if ( group_id == 0 ) {
                 pdgeadd_( "N", &M, &N, 
                           &alpha, U_ac, &i1, &i1, descU_ac, 
                           &alpha, U,    &i1, &i1, descU );
              }
           }
           /**
            * Scale/Symmetrize the matrix on group id 0
            * U = U/ff(1);
            */
           if ( group_id == 0 ) {
              alpha = 1.0; 
              /**
               * normalize so that min(svd(A))=1
               * U = U/ff(1);
               */
              pdlascl_( "G", &ff_1, &alpha, &M, &N, 
                        U, &i1, &i1, descU, 
                        &iinfo ); 
              /**
               * if symm;    U=(U'+U)/2;end    // force symmetry if A symmetric
               */
              if ( symm ){
                 pdlacpy_( "A", &M, &N, 
                           U, &i1, &i1, descU, 
                           B, &i1, &i1, descB );
                 alpha = 0.5; 
                 pdgeadd_( "T", &M, &N, 
                           &alpha, B, &i1, &i1, descB, 
                           &alpha, U, &i1, &i1, descU );
              }
           }
        }
        else {
           alpha = 1.; 
           for ( j = 0; j < nbprob; j++ ){
                 MPI_Barrier( new_comm );
                 pdgemr2d_( &M, &N, 
                            U,     &i1, &i1, descU, 
                            B_all, &i1, &i1, descB_all, 
                            &ictxt_all );
                 beta = 1.; if (j ==0 ) {beta = 0.;}
                 pdgeadd_( "N", &M, &N, 
                           &alpha, B_all, &i1, &i1, descB_all, 
                           &beta,  A_all, &i1, &i1, descA_all );
           }
           alpha = 1.0;
           //descA_all[1]=-1; 

           pdlascl_( "G", &ff_1, &alpha, &M, &N, 
                     A_all, &i1, &i1, descA_all, 
                     &iinfo ); //U = U/ff(1);// normalize so that min(svd(A))=1
           if ( symm ){
              pdlacpy_( "A", &M, &N, 
                        A_all, &i1, &i1, descA_all, 
                        B_all, &i1, &i1, descB_all );
                        alpha = 0.5; 
              pdgeadd_( "T", &M, &N, 
                        &alpha, B_all, &i1, &i1, descB_all, 
                        &alpha, A_all, &i1, &i1, descA_all );
           }
        }
        MPI_Barrier( new_comm );


        /*
         * Bcast 
         */
        MPI_Barrier( new_comm );
        if ( it1 < itmax ){ 
           if ( group_id == 0 ){
              pdlacpy_( "A", &M, &N, 
                        U,    &i1, &i1, descU, 
                        U_ac, &i1, &i1, descU_ac );
           }

           for ( j = 1; j < nbprob; j++ ){
              MPI_Barrier( new_comm );
              pdgemr2d_( &M, &N, 
                         U_ac, &i1, &i1, descU_ac, 
                         U,    &i1, &i1, descU, 
                         &ictxt_all );
              MPI_Barrier( new_comm );
           }
           MPI_Barrier( new_comm );
        }
    }    
    if ( wantH ){ 
        if ( prof ) {Htime =- MPI_Wtime();}
        alpha = 1.0; beta = 0.0;
        pdgemm_( "T", "N", &M, &N, &N, 
                 &alpha, A_all,    &i1, &i1, descA_all, 
                 Acpy_all, &i1, &i1, descAcpy_all, 
                 &beta,  B_all, &i1, &i1, descB_all );
        pdlacpy_( "A", &M, &N, 
                  B_all,    &i1, &i1, descB_all, 
                  Acpy_all, &i1, &i1, descAcpy_all );
        alpha = 0.5; 
        pdgeadd_( "T", &M, &N,  
                  &alpha, Acpy_all, &i1, &i1, descAcpy_all, 
                  &alpha, B_all,    &i1, &i1, descB_all );
      
        /* Main flops used in this step: Calculating H */
        flops_dgemm  = FLOPS_DGEMM( M, N, N );
        *flops += flops_dgemm;
        if ( prof ) {Htime += MPI_Wtime();
           //MPI_Allreduce( &Htime, &reduced_Htime, 1, MPI_DOUBLE, MPI_MAX, new_comm );
           if ( myrank_mpi == 0 ) {
              fprintf(stderr, "# \treduced_Htime     \n");
              fprintf(stderr, "  \t%2.4e\n", reduced_Htime);
           }
        }
    }
} //end else



    MPI_Barrier( MPI_COMM_WORLD );
    return 0;
}




















