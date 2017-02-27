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


/* Default values of parameters */
int nprow         = 1;
int npcol         = 1;
int lvec          = 1;
int rvec          = 1;
int n             = 5120;
int nb            = 128;
int mode          = 4;
double cond       = 9.0072e+15;
int optcond       = 0;
int start         = 5120;
int stop          = 5120;
int step          = 1;
int niter         = 1;
int check         = 0;
int verbose       = 0;

static inline double cWtime(void)
{
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

void print_usage(void)
{
    fprintf(stderr,
            "======= QDWHsvd timing using ScaLAPACK\n"
            " -p      --nprow         : Number of MPI process rows\n"
            " -q      --npcol         : Number of MPI process cols\n"
            " -jl     --lvec          : Compute left singular vectors\n"
            " -jr     --rvec          : Compute right singular vectors\n"
            " -n      --N             : Dimension of the matrix\n"
            " -b      --nb            : Block size\n"
            " -m      --mode          : [1:6] Mode from pdlatms used to generate the matrix\n"
            " -k      --cond          : Condition number used to generate the matrix\n"
            " -o      --optcond       : Estimate Condition number using QR\n"
            " -i      --niter         : Number of iterations\n"
            " -r      --n_range       : Range for matrix sizes Start:Stop:Step\n"
            " -c      --check         : Check the solution\n"
            " -v      --verbose       : Verbose\n"
            " -h      --help          : Print this help\n" );
}

#define GETOPT_STRING "p:q:x:y:n:b:m:i:o:r:Q,S:s:w:e:c:f:t:v:h"

static struct option long_options[] =
    {
        /* PaRSEC specific options */
        {"nprow",      required_argument,  0, 'p'},
        {"npcol",      required_argument,  0, 'q'},
        {"jl",         no_argument,        0, 'x'},
        {"lvec",       no_argument,        0, 'x'},
        {"jr",         no_argument,        0, 'y'},
        {"rvec",       no_argument,        0, 'y'},
        {"N",          required_argument,  0, 'n'},
        {"n",          required_argument,  0, 'n'},
        {"nb",         required_argument,  0, 'b'},
        {"b",          required_argument,  0, 'b'},
        {"mode",       required_argument,  0, 'm'},
        {"m",          required_argument,  0, 'm'},
        {"cond",       required_argument,  0, 'k'},
        {"k",          required_argument,  0, 'k'},
        {"optcond",    required_argument,  0, 'o'},
        {"o",          required_argument,  0, 'o'},
        {"i",          required_argument,  0, 'i'},
        {"niter",      required_argument,  0, 'i'},
        {"r",          required_argument,  0, 'r'},
        {"n_range",    required_argument,  0, 'r'},
        {"check",      no_argument,        0, 'c'},
        {"verbose",    no_argument,        0, 'v'},
        {"help",       no_argument,        0, 'h'},
        {"h",          no_argument,        0, 'h'},
        {0, 0, 0, 0}
    };

static void parse_arguments(int argc, char** argv)
{
    int opt = 0;
    int c;
    int myrank_mpi;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);

    do {
#if defined(HAVE_GETOPT_LONG)
        c = getopt_long_only(argc, argv, "",
                        long_options, &opt);
#else
        c = getopt(argc, argv, GETOPT_STRING);
        (void) opt;
#endif  /* defined(HAVE_GETOPT_LONG) */

        switch(c) {
        case 'p': nprow     = atoi(optarg); break;
        case 'q': npcol     = atoi(optarg); break;
        case 'n': n         = atoi(optarg); start = n; stop = n; step = 1; break;
        case 'b': nb        = atoi(optarg); break;
        case 'm': mode      = atoi(optarg); break;
        case 'k': cond      = atof(optarg); break;
        case 'o': optcond   = atof(optarg); break;
        case 'i': niter     = atoi(optarg); break;
        case 'r': get_range( optarg, &start, &stop, &step ); break;
        case 'c': check     = 1; break;
        case 'v': verbose   = 1; break;
        case 'h':
            if (myrank_mpi == 0) print_usage(); MPI_Finalize(); exit(0);
            break;
        default:
            break;
        }
    } while(-1 != c);
}

int main(int argc, char **argv) {

	int myrank_mpi, nprocs_mpi;
	int ictxt, myrow, mycol;
	int mloc, nloc, mlocW;
        int mpi_comm_rows, mpi_comm_cols;
        int useQr, THIS_REAL_ELPA_KERNEL_API;
	int i, j, k, iter, size, info_facto_mr, info_facto_dc, info_facto_qw, info_facto_el, info_facto_sl, info, iseed;
	int my_info_facto;
        int i0 = 0, i1 = 1;
        int lwork, liwork, *iWloc=NULL, ldw;
        long int LDW;

	int descA[9], descAcpy[9], descU[9], descVT[9], descWglo[9], descC[9], descH[9], descSigma[9];
        int descCcpy[9], descAu[9];

	double *A=NULL, *Acpy=NULL, *Wloc=NULL, *H=NULL, *D=NULL, *C=NULL;
	double *tau=NULL, *Wglo=NULL;
        double *Sigma=NULL;

        double eps = LAPACKE_dlamch_work('e');
        int iprepad, ipostpad, sizemqrleft, sizemqrright, sizeqrf, sizeqtq,
                     sizechk, sizesyevx, isizesyevx,
                     sizesubtst, isizesubtst, sizetst,
                     isizetst;
/**/
        double my_elapsed_polarqdwh = 0.0, elapsed_polarqdwh = 0.0, sumtime_polarqdwh = 0.0;
        double max_time_polarqdwh  = 0.0, min_time_polarqdwh  = 1e20;

        double flops, GFLOPS;

        double orth_Uqw, berr_UHqw;
        double orth_Usvd, berr_UHsvd, frobA;

        int success;
        double alpha, beta;
        char *jobu, *jobvt;


/**/

        if (verbose & myrank_mpi == 0) fprintf(stderr, "Program starts... \n");

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

        if (verbose & myrank_mpi == 0) fprintf(stderr, "MPI Init done\n");
        parse_arguments(argc, argv);
        if (verbose & myrank_mpi == 0) fprintf(stderr, "Checking arguments done\n");

        Cblacs_get( -1, 0, &ictxt );
	Cblacs_gridinit( &ictxt, "R", nprow, npcol );
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
        if (myrank_mpi == 0) printf("\n =================== nprow %d npcol %d \n", nprow, npcol);
        if (verbose & myrank_mpi == 0) fprintf(stderr, "BLACS Init done\n");

	if (myrank_mpi == 0) {
           fprintf(stderr, "# \n");
           fprintf(stderr, "# NPROCS %d P %d Q %d\n", nprocs_mpi, nprow, npcol);
           fprintf(stderr, "# niter %d\n", niter);
           fprintf(stderr, "# n_range %d:%d:%d mode: %d cond: %2.4e \n", start, stop, step, mode, cond);
           fprintf(stderr, "# \n");
        }

        if (verbose & myrank_mpi == 0) fprintf(stderr, "Range loop starts\n");


        // Begin loop over range
        for (size = start; size <= stop; size += step) {
            while ( (int)((double)size / (double)nb) < ( max(nprow , npcol) )){
               if (myrank_mpi == 0) fprintf(stderr, " Matrix size is small to be facrorized using this number of processors \n");
               size += step;
            }
            n = size; ldw = 2*n, LDW = ldw*n; long int matsize = n*n;
	    mloc  = numroc_( &n, &nb, &myrow, &i0, &nprow );
	    nloc  = numroc_( &n, &nb, &mycol, &i0, &npcol );
	    mlocW = numroc_( &ldw, &nb, &myrow, &i0, &nprow );
            if (verbose & myrank_mpi == 0) fprintf(stderr, "Desc Init starts %d\n", mloc);

	    descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    descinit_( descAcpy, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    descinit_( descC, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    descinit_( descH, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    descinit_( descWglo, &ldw, &n, &nb, &nb, &i0, &i0, &ictxt, &mlocW, &info );
            if (verbose & myrank_mpi == 0) fprintf(stderr, "Desc Init ends %d\n", mloc);

	    A     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    Acpy  = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    C     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    H     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    Wglo  = (double *)malloc(mlocW*nloc*sizeof(double)) ;
	    tau   = (double *)malloc(nloc*sizeof(double)) ;
	    D     = (double *)malloc(n*sizeof(double)) ;

	    descinit_( descSigma, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    Sigma = (double *)calloc(mloc*nloc,sizeof(double)) ;

            /* Initialize the timing counters */
            my_elapsed_polarqdwh = 0.0, elapsed_polarqdwh = 0.0, sumtime_polarqdwh = 0.0; 
            max_time_polarqdwh  = 0.0, min_time_polarqdwh  = 1e20;

        
            /* Generate matrix by pdlatms */
            {
	       //D  = (double *)malloc(n*sizeof(double)) ;
               char   *dist = "N"; /* NORMAL( 0, 1 )  ( 'N' for normal ) */
               int    iseed[4] = {1, 0, 0, 1};
               char   *sym = "P"; /* The generated matrix is symmetric, with
                                    eigenvalues (= singular values) specified by D, COND,
                                    MODE, and DMAX; they will not be negative.
                                    "N" not supported. */
               //int    mode = 4; /* sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND) */
               //double cond = 1.0/eps;
               double dmax = 1.0;
               int    kl   = n;
               int    ku   = n;
               char   *pack = "N"; /* no packing */
               int    order = n;
               int    info;
         
               pdlasizesep_( descA, 
                             &iprepad, &ipostpad, &sizemqrleft, &sizemqrright, &sizeqrf, 
                             &lwork, 
                             &sizeqtq, &sizechk, &sizesyevx, &isizesyevx, &sizesubtst, 
                             &isizesubtst, &sizetst, &isizetst );
               if (verbose & myrank_mpi == 0) fprintf(stderr, "Setting lwork done\n");
               Wloc = (double *)calloc(lwork,sizeof(double)) ;
               pdlatms_(&n, &n, dist,
                        iseed, sym, D, &mode, &cond, &dmax,
                        &kl, &ku, pack, 
                        A, &i1, &i1, descA, &order, 
                        Wloc, &lwork, &info);
               if (verbose & myrank_mpi == 0) fprintf(stderr, "MatGen done\n");
               if (info != 0) {
                   fprintf(stderr, "An error occured during matrix generation: %d\n", info );
                   return EXIT_FAILURE;
               }
               pdlacpy_( "All", &n, &n, 
                         A, &i1, &i1, descA, 
                         Acpy, &i1, &i1, descAcpy ); 
               frobA  = pdlange_ ( "f", &n, &n, A, &i1, &i1, descA, Wloc);
               beta = 0.0;
               pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);

               for (i = 1; i <= n; i++) {
                   int idum1, idum2, iloc, jloc;
                   if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                         &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                               iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                               jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                               Sigma[ (jloc-1)*mloc + (iloc-1) ] = D[i-1];
                   }
               } 

               //norm_sv    = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc);
               if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy to Acpy done\n");

               free( Wloc );
            }

            if (myrank_mpi == 0) fprintf(stderr, "\n\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

            // QDWH
	    //if ( qwmr || qwdc || qwel || polarqdwh) {
               jobu   = lvec ? "V" : "N";
               Wloc  = (double *)calloc(1,sizeof(double)) ;
               iWloc = (int *)calloc(1,sizeof(int)) ;
               int lwork_qr;
               int lwork_cn, liwork_cn;

               lwork_cn = -1; liwork_cn = -1;
               double Anorm = 1., Li = 1.;
               pdgecon_ ("1", &n, A, &i1, &i1, descA, 
                         &Anorm, &Li, 
                         Wloc, &lwork_cn, iWloc, &liwork_cn, &my_info_facto);
               lwork_cn  = Wloc[0];
               liwork_cn = (int)iWloc[0];

               lwork_qr = -1; 
               pdgeqrf_(&ldw, &n, Wglo, &i1, &i1, descWglo, 
                         tau, Wloc, &lwork_qr, &info);
               lwork_qr  =  Wloc[0];

               lwork     = max ( lwork_qr, lwork_cn);
	       Wloc  = (double *)malloc(lwork*sizeof(double)) ;
	       iWloc = (int *)malloc(n*sizeof(int)) ;

               for (iter = 0; iter < niter; iter++) {
                  pdlacpy_( "All", &n, &n, 
                            Acpy, &i1, &i1, descAcpy, 
                            A,    &i1, &i1, descA ); 
                  flops = 0.0;
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "\nQDWH + ScaLAPACK EIG done\n");
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "QDWH starts...\n");

                  /*
                   * Find polar decomposition using QDWH. 
                   * C contains the positive-definite factor. 
                   * A contains the orthogonal polar factor.
                   */
                  my_elapsed_polarqdwh   = 0.0;
                  my_elapsed_polarqdwh   =- MPI_Wtime();

                  my_info_facto = pdgeqdwh( verbose, 0, n, n, optcond, 
                                        A, descA, 
                                        H, descH, 
                                        Wglo, descWglo, LDW, 
                                        C, descC, 
                                        tau,
                                        Wloc, lwork,
                                        iWloc, liwork_cn, 
                                        ictxt, &flops);

                  my_elapsed_polarqdwh   += MPI_Wtime();
                  MPI_Allreduce( &my_elapsed_polarqdwh, &elapsed_polarqdwh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                  sumtime_polarqdwh += elapsed_polarqdwh;
                  if ( elapsed_polarqdwh >= max_time_polarqdwh ) { max_time_polarqdwh = elapsed_polarqdwh;} 
                  if ( elapsed_polarqdwh <= min_time_polarqdwh ) { min_time_polarqdwh = elapsed_polarqdwh;} 

	           MPI_Allreduce( &my_info_facto, &info_facto_qw, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

                  if (verbose & myrank_mpi == 0) fprintf(stderr, "QDWH ends...\n");
                  if (verbose & check & myrank_mpi == 0) fprintf(stderr, "Testing QDWH starts...\n");
                  /*
                   * Checking the polar factorization
                   */
                  if(check && iter == 0){
                       /*
                        * checking orthogonality of Up
                        */ 
                        alpha = 0.0; beta = 1.0;
                        pdlaset_( "G", &n, &n, &alpha, &beta, H, &i1, &i1, descH);
                        alpha = 1.0; beta = -1.0;
                        pdgemm_( "T", "N", &n, &n, &n, 
                                 &alpha, 
                                 A, &i1, &i1, descA, 
                                 A, &i1, &i1, descA, 
                                 &beta, 
                                 H, &i1, &i1, descH);
                        orth_Uqw  = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc)/frobA;
                       /*
                        * checking the factorization |A-Up*H|
                        */ 
                        pdlacpy_( "A", &n, &n, Acpy, &i1, &i1, descAcpy, H, &i1, &i1, descH );
                        pdgemm_( "N", "N", &n, &n, &n, 
                                 &alpha, 
                                 A, &i1, &i1, descA, 
                                 C, &i1, &i1, descC, 
                                 &beta, 
                                 H, &i1, &i1, descH);
                        berr_UHqw  = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc)/frobA;
                  }
               }
               if (  myrank_mpi == 0) {
                  fprintf(stderr, "# QDWH \n"); 
                  fprintf(stderr, "#\n");
	          fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ   \tGflop/s \tAvg-Time     \tMax-Time    \tMin-Time    \tBerr_UpH  \tOrth_Up  \tinfo     \n");
	          fprintf(stderr, "   %6d \t%4d \t%4d \t%3d \t%3d \t%8.2f", n, nb, nprocs_mpi, nprow, npcol, flops/1e9/min_time_polarqdwh);
	          fprintf(stderr, "\t%6.2f \t\t%6.2f \t\t%6.2f \t\t%2.4e \t%2.4e \t%d \n", sumtime_polarqdwh/niter, max_time_polarqdwh, min_time_polarqdwh, berr_UHqw, orth_Uqw,info_facto_qw);
                  fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
               }


            //}

            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

	    free( Sigma );
	    free( A );
	    free( Acpy );
	    free( C );
	    free( H );
	    free( Wglo );
	    free( tau );
	    free( D );
            if (verbose & myrank_mpi == 0) fprintf(stderr, "Free matrices done\n");
        } // End loop over range


        if (verbose & myrank_mpi == 0) fprintf(stderr, "Range loop ends\n");

	blacs_gridexit_( &i0 );
	MPI_Finalize();
        if (verbose & myrank_mpi == 0) fprintf(stderr, "Program ends...\n");
	return 0;
}
