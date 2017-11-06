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
 *  PDGSVD computes the SVD of a real distributed M-by-N based
 *  on the polar decomposition QDWH
 *
 *  matrix A = U * Sigma * VT.
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
 * JOBU    (global input) CHARACTER*1
 *          Specifies options for computing U:
 *          = 'V':  the first SIZE columns of U (the left singular
 *                  vectors) are returned in the array U;
 *          = 'N':  no columns of U (no left singular vectors) are
 *                  computed.
 *
 * JOBVT   (global input) CHARACTER*1
 *          Specifies options for computing U:
 *          = 'V':  the first SIZE columns of V (the right singular
 *                  vectors) are returned in the array VT;
 *          = 'N':  no columns of VT (no right singular vectors) are
 *                  computed.
 *
 * eigtype  (global input) CHARACTER*1
 *          Specifies the eigensolver to be used after the polar decomposition:
 *          = 'r': Use PDSYEVR 
 *          = 'd': Use PDSYEVD 
 *          = 'r': Use ELPA-2stage 
 *
 *  M       (global input) INTEGER
 *          The number of rows of the input matrix A.  M >= 0.
 *
 *  N       (global input) INTEGER
 *          The number of columns of the input matrix A.  N >= 0.
 *
 *  A       (local input/workspace) block cyclic DOUBLE PRECISION
 *          array,
 *          global dimension (M, N), local dimension (MP, NQ)
 *          On exit, the contents of A are destroyed.
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
 *  S       (global output) DOUBLE PRECISION   array, dimension SIZE
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 *  U       (local output) DOUBLE PRECISION   array, local dimension
 *          (MP, SIZEQ), global dimension (M, SIZE)
 *          if JOBU = 'V', U contains the first min(m,n) columns of U
 *          if JOBU = 'N', U is not referenced.
 *
 *  IU      (global input) INTEGER
 *          The row index in the global array U indicating the first
 *          row of sub( U ).
 *
 *  JU      (global input) INTEGER
 *          The column index in the global array U indicating the
 *          first column of sub( U ).
 *
 *  DESCU   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix U.
 *
 *  VT      (local output) DOUBLE PRECISION   array, local dimension
 *          (SIZEP, NQ), global dimension (SIZE, N).
 *          If JOBVT = 'V', VT contains the first SIZE rows of
 *          V**T. If JOBVT = 'N', VT is not referenced.
 *
 *  IVT     (global input) INTEGER
 *          The row index in the global array VT indicating the first
 *          row of sub( VT ).
 *
 *  JVT     (global input) INTEGER
 *          The column index in the global array VT indicating the
 *          first column of sub( VT ).
 *
 *  DESCVT   (global input) INTEGER array of dimension DLEN_
 *          The array descriptor for the distributed matrix VT.
 *
 *  WORK    (local workspace/output) DOUBLE PRECISION   array, dimension
 *          (LWORK)
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (local input) INTEGER
 *          The dimension of the array WORK.
 *  if eigtype == 'r' then
 *  LWORK   (local input) INTEGER
 *          Size of WORK, must be at least 3.
 *          See below for definitions of variables used to define LWORK.
 *          If no eigenvectors are requested (JOBZ = 'N') then
 *             LWORK >= max(2 + 5*N + MAX( 12 * NN, NB * ( NP0 + 1 ) ), LOCrW*LOCc)
 *          If eigenvectors are requested (JOBZ = 'V' ) then
 *             the amount of workspace required is:
 *             LWORK >= max(2 + 5*N + MAX( 18*NN, NP0 * MQ0 + 2 * NB * NB ) +
 *               (2 + ICEIL( NEIG, NPROW*NPCOL))*NN, LOCrW*LOCc)
 *
 *  if eigtype == 'd' then
 *  LWORK   (local input) INTEGER
 *          LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
 *          TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
 *          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
 *          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
 *          LWORK >= max(1 + 6*SIZEB + MAX(WATOBD, WBDTOSVD), LOCrW*LOCc)
 *
 *          If LWORK = -1, then LWORK is global input and a workspace
 *          query is assumed; the routine only calculates the minimum
 *          size for the work array. The required workspace is returned
 *          as the first element of WORK and no error message is issued
 *          by PXERBLA.
 *          Where, MLOCW is 
 *          LOCrW( M ) = NUMROC( 2*M, MB_A, MYROW, RSRC_A, NPROW ),
 *
 *
 *  if eigtype == 'd' then
 *  The WORK is NULL and not referenced
 *
 *  IWLOC   (local workspace/output) INTEGER array, dimension (LIWORK)
 *          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
 *
 *  LIWORK  (input) INTEGER
 *          The dimension of the array IWORK.
 *          PDGETRF_LIWORK = ( LOCr(M_A)+MB_A )
 *          PDGECON_LIWORK >= MAX( 1, LOCr(N+MOD(IA-1,MB_A)) ).
 *          Let  NNP = MAX( N, NPROW*NPCOL + 1, 4 ). Then:
 *          PDSYEVR_LIWORK >= 12*NNP + 2*N when the eigenvectors are desired
 *          PDSYEVR_LIWORK >= 10*NNP + 2*N when only the eigenvalues have to be computed
 *          PDSYEVD_LIWORK = 7*N + 8*NPCOL + 2. 
 *          LIWORK = MAX(PDGETRF_LIWORK, PDGECON_LIWORK , EIGTYPE_LIWORK ) 
 *
 *          If LIWORK = -1, then LIWORK is global input and a workspace
 *          query is assumed; the routine only calculates the minimum
 *          and optimal size for all work arrays. Each of these
 *          values is returned in the first entry of the corresponding
 *          work array, and no error message is issued by PXERBLA.
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *         != 0:  the eigensolver fail
 *
 ******************************************************************************/


int pdgeqsvd( char *jobu, char *jobvt, char *eigtype, 
              int m, int n, 
              double *A, int iA, int jA, int *descA, 
              double *S, 
              double *U,     int iU,     int jU, int *descU,
              double *VT,    int iVT,    int jVT, int *descVT,
              double *Work,  int lWork,
              int    *iWork, int liWork, int *info)
{

    int verbose = 0; int profqw = 0; int optcond = 0;
    int vl, vu, il, iu, nbeigvals, nbeigvecs;
    double flops, GFLOPS;
    flops = 0.0;

    int i0 = 0;
    int i1 = 1;

    int mloc, nloc, mlocW, nb;   
    int myrow, mycol, nprow, npcol;   
    int ctxt_ = 1, nb_ = 5;
    int ictxt;
    int MB = 2*n; 
           
    /*
     * Get the grid parameters
     */
    ictxt = descU[ctxt_];
    Cblacs_get( -1, 0, &ictxt );
    nb = descU[nb_];
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    mloc  = numroc_( &m, &nb, &myrow, &i0, &nprow );
    nloc  = numroc_( &n, &nb, &mycol, &i0, &npcol );
    mlocW = numroc_( &MB, &nb, &myrow, &i0, &nprow );


    double alpha = 1.0, beta = 0.0;
    int lwork_cn, liwork_cn;

   /*
    * Find Workspace 
    */
    if (lWork  == -1 && liWork == -1){
        double Anorm = 1., Li = 1.;
        lwork_cn = -1; liwork_cn = -1;
        //pdgecon_ ("1", &n, U, &iU, &jU, descU, 
        //          &Anorm, &Li, 
        //          Work, &lWork, iWork, &liWork, info);
        //liwork_cn = n;//(int)iWork[0];
        lwork_cn  = Work[0];
        liwork_cn = (int)iWork[0];


        if (eigtype[0] == 'r') {
            pdsyevr_( "V", "A", "L", &n, 
                  U, &iU, &jU, descU, 
                  &vl, &vu, &il, &iu, &nbeigvals, &nbeigvecs,
                  S, 
                  VT, &iVT, &jVT, descVT, 
                  Work, &lWork, 
                  iWork, &liWork, info );
        }   
        else if (eigtype[0] == 'd') {
            pdsyevd_( jobvt, "L", &n, 
                  U, &iU, &jU, descU, 
                  S, 
                  VT, &iVT, &jVT, descVT, 
                  Work, &lWork, 
                  iWork, &liWork, info );
        }   
        //lWork  = max ( Work[0], lwork_cn);
        lWork  = max ( Work[0], mlocW*nloc);
        liWork = max ( (int)iWork[0], liwork_cn);
        Work[0]  = lWork;
        iWork[0] = liWork;
        return 0;
    }

    pdgeqdwh( n, n,
              A, iA, jA, descA, // UP 
              U, iU, jU, descU, // H 
              VT, mloc,
              Work, mlocW,
              &info);

    if (eigtype[0] == 'r'){
        pdsyevr_( "V", "A", "L", &n, 
                   U, &iU, &jU, descU, 
                   &vl, &vu, &il, &iu, &nbeigvals, &nbeigvecs,
                   S, 
                   VT, &iVT, &jVT, descVT, 
                   Work, &lWork, 
                   iWork, &liWork, info );
              
    }
    else if(eigtype[0] == 'd'){
        pdsyevd_( jobvt, "L", &n, 
                  U, &iU, &jU, descU, 
                  S, 
                  VT, &iVT, &jVT, descVT, 
                  Work, &lWork, 
                  iWork, &liWork, info );
    }
    else if(eigtype[0] == 'e'){
        Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
	mloc  = numroc_( &n, &nb, &myrow, &i0, &nprow );
	nloc  = numroc_( &n, &nb, &mycol, &i0, &npcol );
        int useQr, THIS_REAL_ELPA_KERNEL_API;
        int mpi_comm_rows, mpi_comm_cols;
        int mpierr = elpa_get_communicators(MPI_Comm_c2f(MPI_COMM_WORLD), myrow, mycol, &mpi_comm_rows, &mpi_comm_cols);
        useQr = 0;
        THIS_REAL_ELPA_KERNEL_API = ELPA2_REAL_KERNEL_AVX_BLOCK6;
        *info = elpa_solve_evp_real_2stage( n, n, U, mloc, 
                                            S, VT, 
                                            mloc, nb, nloc, 
                                            mpi_comm_rows, mpi_comm_cols, MPI_Comm_c2f(MPI_COMM_WORLD),
                                            THIS_REAL_ELPA_KERNEL_API, useQr);
    }

    if(jobu = "V") {
        pdgemm_( "N", "N", &n, &n, &n, 
                 &alpha, 
                 A, &iA, &jA, descA, 
                 VT, &iVT, &jVT, descVT, 
                 &beta, 
                 U, &iU, &jU, descU);
    }
    return info[0];
}
