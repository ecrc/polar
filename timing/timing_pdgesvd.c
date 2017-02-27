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

// TO RUN THE POLAR DECOMPOSITION:
// make; mpirun -np 2 ./main_svd --nprow 2 --npcol 1 --b 128 --niter 2 --n_range 1024:1024:1024  --check --polarqdwh --polarsvd

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
int polarsvd     = 0;
int slsvd         = 0;
int qwmr          = 0;
int qwdc          = 0;
int qwel          = 0;
int check         = 0;
int profqw        = 0;
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
            " -polarqdwh --polarqdwh  : Find polar decomposition using QDWH A=UH \n"
            " -polarsvd  --polarsvd   : Find the polar decomposition using scalapack-svd \n"
            " -s      --slsvd         : Run reference ScaLAPACK SVD\n"
            " -w      --qwmr          : Run QDWH SVD with ScaLAPACK MRRR EIG\n"
            " -e      --qwdc          : Run QDWH SVD with ScaLAPACK DC EIG\n"
            " -l      --qwel          : Run QDWH SVD with ScaLAPACK DC EIG\n"
            " -c      --check         : Check the solution\n"
            " -fqwsvd --profqwsvd     : Enable profiling QDWHsvd\n"
            " -fqw    --profqw        : Enable profiling QDWH\n"
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
        //{"polar",      no_argument,        0, 'u'},
        {"polarqdwh",  no_argument,        0, 'Q'},
        {"polarsvd",   no_argument,        0, 'S'},
        {"slsvd",      no_argument,        0, 's'},
        {"qwmr",       no_argument,        0, 'w'},
        {"qwdc",       no_argument,        0, 'e'},
        {"qwel",       no_argument,        0, 'l'},
        {"e",          no_argument,        0, 'c'},
        {"check",      no_argument,        0, 'c'},
        {"profqwsvd",  no_argument,        0, 'f'},
        {"fqwsvd",     no_argument,        0, 'f'},
        {"profqw",     no_argument,        0, 't'},
        {"fqw",        no_argument,        0, 't'},
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
        case 'S': polarsvd  = 1; break;
        case 's': slsvd     = 1; break;
        case 'w': qwmr      = 1; break;
        case 'e': qwdc      = 1; break;
        case 'l': qwel      = 1; break;
        case 'i': niter     = atoi(optarg); break;
        case 'r': get_range( optarg, &start, &stop, &step ); break;
        case 'c': check     = 1; break;
        case 't': profqw    = 1; break;
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
	int i, j, k, iter, size, info_facto, info, iseed;
	int my_info_facto;
        int i0 = 0, i1 = 1;
        int lwork, liwork, *iWloc=NULL, ldw;
        long int LDW;
	int descA[9], descAcpy[9], descU[9], descVT[9], descWglo[9], descH[9], descSigma[9];
	double *A=NULL, *Acpy=NULL, *U=NULL, *VT=NULL, *S=NULL, *Wloc=NULL, *D=NULL, *C=NULL;
	double *H=NULL, *tau=NULL, *Wglo=NULL;
        double *Sigma=NULL;


        double eps = LAPACKE_dlamch_work('e');
        int iprepad, ipostpad, sizemqrleft, sizemqrright, sizeqrf, sizeqtq,
                     sizechk, sizesyevx, isizesyevx,
                     sizesubtst, isizesubtst, sizetst,
                     isizetst;

        double my_elapsed_qwsvd = 0.0, elapsed_qwsvd = 0.0, sumtime_qwsvd = 0.0;
        double my_elapsed_slsvd = 0.0, elapsed_slsvd = 0.0, sumtime_slsvd = 0.0;
        double my_elapsed_polarsvd  = 0.0, elapsed_polarsvd  = 0.0, sumtime_polarsvd  = 0.0;

        double max_time_qwsvd = 0.0, min_time_qwsvd = 1e20;
        double max_time_slsvd = 0.0, min_time_slsvd = 1e20;
        double max_time_polarsvd = 0.0, min_time_polarsvd = 1e20;

        double flops, GFLOPS;

        double berr = 0.0, my_berr = 0.0, norm_sv;
        double my_berr_qwmr = 0.0, my_berr_qwdc = 0.0, my_berr_qwel = 0.0;
        double my_acc_qwmr = 0.0, my_acc_qwdc = 0.0, my_acc_qwel = 0.0, my_acc_slsvd = 0.0;
        double my_orthR_qwmr = 0.0, my_orthR_qwdc = 0.0, my_orthR_qwel = 0.0, my_orthR_slsvd = 0.0;
        double my_orthL_qwmr = 0.0, my_orthL_qwdc = 0.0, my_orthL_qwel = 0.0, my_orthL_slsvd = 0.0;

        double orth_Uqw, berr_UHqw;
        double orth_Usvd, berr_UHsvd, frobA;

        int success;
        double alpha, beta;
        char *jobu, *jobvt, *eigtype;
        int vl, vu, il, iu, nbeigvals, nbeigvecs;

        jobu  = lvec ? "V" : "N";
        jobvt = rvec ? "V" : "N";


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

        /* to run only the polar decompsition */
        if(polarsvd )  {slsvd = 0;} 

        if ( qwmr )
           eigtype = "r";
        else if (qwdc)
           eigtype = "d";
        else if (qwel)
           eigtype = "e";

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
	    descinit_( descH, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    descinit_( descWglo, &ldw, &n, &nb, &nb, &i0, &i0, &ictxt, &mlocW, &info );
	    descinit_( descSigma, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );

	    descinit_( descU, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	    descinit_( descVT, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
            if (verbose & myrank_mpi == 0) fprintf(stderr, "Desc Init ends %d\n", mloc);

	    A     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    Acpy  = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    H     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    Wglo  = (double *)malloc(mlocW*nloc*sizeof(double)) ;
	    tau   = (double *)malloc(nloc*sizeof(double)) ;
	    D  = (double *)malloc(n*sizeof(double)) ;
	    Sigma = (double *)calloc(mloc*nloc,sizeof(double)) ;


	    U      = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    VT      = (double *)malloc(mloc*nloc*sizeof(double)) ;
	    S      = (double *)malloc(n*sizeof(double)) ;

            /* Initialize the timing counters */
            my_elapsed_qwsvd = 0.0, elapsed_qwsvd = 0.0, sumtime_qwsvd = 0.0;
            my_elapsed_slsvd = 0.0, elapsed_slsvd = 0.0, sumtime_slsvd = 0.0;
            my_elapsed_polarsvd  = 0.0, elapsed_polarsvd  = 0.0, sumtime_polarsvd  = 0.0;

            max_time_qwsvd = 0.0, min_time_qwsvd = 1e20;
            max_time_slsvd = 0.0, min_time_slsvd = 1e20;
            max_time_polarsvd = 0.0, min_time_polarsvd = 1e20;
        
            /* Generate matrix by pdlatms */
            {
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

               norm_sv    = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc);
               if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy to Acpy done\n");

               free( Wloc );
            }

            if (myrank_mpi == 0) fprintf(stderr, "\n\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

            // QDWH + EIG
	    if ( qwmr || qwdc || qwel ) {
               /*
                * SVD decomposition is (A *U)* S * U' = C * S * U'.
                */
               for (iter = 0; iter < niter; iter++) {
                  flops = 0.0;

                  if( (qwmr || qwdc || qwel) ){
                      pdlacpy_( "A", &n, &n, A, &i1, &i1, descA, Acpy, &i1, &i1, descAcpy );
                  }

                  /*
                   * Find the SVD using QDWH + EIG 
                   */
	          if ( qwmr || qwdc || qwel ) {
                     if (verbose & myrank_mpi == 0) fprintf(stderr, "EIG starts...\n");
                     /*
                      * Find Workspace 
                      */
                     lwork  = -1; liwork = -1;
                     Wloc   = (double *)calloc(1,sizeof(double));
                     iWloc  = (int *)calloc(1,sizeof(int));
                     pdgeqsvd( jobu, jobvt, eigtype,  
                                nprow, npcol, nb, ictxt, 
                                n, n, 
                                A, i1, i1, descA, 
                                S, 
                                U,     i1,     i1,  descU,
                                VT,    i1,     i1,  descVT,
                                Wglo,  i1,     i1,  descWglo,
                                Wloc,  lwork,
                                iWloc, liwork, &my_info_facto);
                     lwork  = (int)Wloc[0];
                     liwork = (int)iWloc[0];
	             Wloc   = (double *)calloc(lwork,sizeof(double)) ;
	             iWloc  = (int *)calloc(liwork,sizeof(int)) ;
                     /*
                      * QDWH + EIG
                      */

                     my_elapsed_qwsvd   = 0.0;
                     my_elapsed_qwsvd   =- MPI_Wtime();
                     pdgeqsvd( jobu, jobvt, eigtype,
                                nprow, npcol, nb, ictxt, 
                                n, n, 
                                A, i1, i1, descA, 
                                S, 
                                U,     i1,     i1,  descU,
                                VT,    i1,     i1,  descVT,
                                Wglo,  i1,     i1,  descWglo,
                                Wloc,  lwork,
                                iWloc, liwork, &my_info_facto);
                     my_elapsed_qwsvd   += MPI_Wtime();
                     MPI_Allreduce( &my_elapsed_qwsvd, &elapsed_qwsvd, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                     sumtime_qwsvd += elapsed_qwsvd;
                     if ( elapsed_qwsvd >= max_time_qwsvd ) { max_time_qwsvd = elapsed_qwsvd;} 
                     if ( elapsed_qwsvd <= min_time_qwsvd ) { min_time_qwsvd = elapsed_qwsvd;} 

                     if (verbose & myrank_mpi == 0) fprintf(stderr, "Compute left singular vectors end...\n");
   
	             MPI_Allreduce( &my_info_facto, &info_facto, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

                     if (verbose & myrank_mpi == 0) fprintf(stderr, "\nQDWH + ScaLAPACK EIG done\n");
         
                     /*
                      * Checking the SVD decomposition
                      */
                     if (check ) {
                        if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing starts...\n");
                        alpha = 1.0; beta = 0.0;
                        /*
                         * Set the singular values on the main diagonal 
                         * |A - U*Sigma*V'|
                         */
                        pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);

                        for (i = 1; i <= n; i++) {
                                int idum1, idum2, iloc, jloc;
                                if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                                &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                        iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                        jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                        Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                                }
                        } 


                        pdlacpy_( "All", &n, &n, 
                                 Acpy, &i1, &i1, descAcpy, 
                                 H, &i1, &i1, descH ); 
                        alpha = 1.0; beta = 0.0;
                        pdgemm_( "N", "N", &n, &n, &n, 
                                 &alpha, 
                                 U   , &i1, &i1, descU, 
                                 Sigma, &i1, &i1, descSigma, 
                                 &beta, 
                                 A, &i1, &i1, descA);
                        beta = -1.0;
                        pdgemm_( "N", "T", &n, &n, &n, 
                                 &alpha, 
                                 A, &i1, &i1, descA, 
                                 VT, &i1, &i1, descVT, 
                                 &beta, 
                                 H, &i1, &i1, descH);
                        my_berr_qwmr = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc) / (frobA * n);
                        
                        /* 
                         * Accuracy of singular values 
                         */
                        //for(i=0; i < n ; i++ )
                        //    D[i] = fabs(D[i]);
                        dlasrt_( "D", &n, S, &info );
                        dlasrt_( "D", &n, D, &info );
                        for(i=0; i < n ; i++ )
                            S[i] = S[i] - D[i];
                        for (i = 1; i <= n; i++) {
                             int idum1, idum2, iloc, jloc;
                             if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                             &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                     iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                     jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                     Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                             }
                        } 
                        my_acc_qwmr = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc) / norm_sv;

                        /* 
                         * Orthogonality of Left singular vectors 
                         */
                        alpha = 0.0; beta = 1.0;
                        pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                        alpha = 1.0; beta = -1.0;
                        pdgemm_( "T", "N", &n, &n, &n, 
                                 &alpha, 
                                 U, &i1, &i1, descU, 
                                 U, &i1, &i1, descU, 
                                 &beta, 
                                 Sigma, &i1, &i1, descSigma);
                        my_orthL_qwmr = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                        /* 
                         * Orthogonality of Right singular vectors 
                         */
                        alpha = 0.0; beta = 1.0;
                        pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                        alpha = 1.0; beta = -1.0;
                        pdgemm_( "T", "N", &n, &n, &n, 
                                 &alpha, 
                                 VT, &i1, &i1, descVT, 
                                 VT, &i1, &i1, descVT, 
                                 &beta, 
                                 Sigma, &i1, &i1, descSigma);
                        my_orthR_qwmr = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                        if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing ends...\n");
                     }
                     free(Wloc); free(iWloc);
                  }
                  if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

                  /*
                   * Save copy of A in Acpy  
                   */
                  pdlacpy_( "All", &n, &n, 
                            Acpy, &i1, &i1, descAcpy, 
                            A, &i1, &i1, descA ); 
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy back to A done\n");
               }
               if ( (qwmr || qwdc || qwel)  && myrank_mpi == 0) {
                  fprintf(stderr, "# QDWH + EIG\n"); 
                  fprintf(stderr, "#\n");
	          fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ   \tAvg-Time     \tMax-Time    \tMin-Time  \tinfo     \tAcc-sv    \tOrth-Rsv    \tOrth-Lsv     \tBerr\n");
	          fprintf(stderr, "   %6d \t%4d \t%3d \t%3d \t%3d", n, nb, nprocs_mpi, nprow, npcol);
	          fprintf(stderr, "\t%6.2f \t\t%6.2f \t\t%6.2f \t\t%d \t\t%2.4e \t%2.4e \t%2.4e \t%2.4e\n", 
                                  sumtime_qwsvd/niter, max_time_qwsvd, min_time_qwsvd, info_facto, my_acc_qwmr, my_orthR_qwmr, my_orthL_qwmr, my_berr_qwmr);
               }

	       my_berr = 0.0;
	       my_berr_qwmr = 0.0;
	       my_berr_qwdc = 0.0;
            }

            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

            /*
             * ScaLAPACK SVD
             */
	    if ( slsvd || polarsvd) {
               lwork = -1;
               Wloc  = (double *)calloc(1,sizeof(double));
               if( polarsvd){
                  jobu  = "V";
                  jobvt = "V";
               }

               pdgesvd_( jobu, jobvt, &n, &n, A, &i1, &i1, descA, 
                         S, 
                         U, &i1, &i1, descU, 
                         VT, &i1, &i1, descVT, 
                         Wloc, &lwork, &my_info_facto );
               lwork = (int)Wloc[0];
	       Wloc  = (double *)calloc(lwork,sizeof(double)) ;
 
               for (iter = 0; iter < niter; iter++) {
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "\nScaLAPACK dgesvd starts...\n");
                  if( polarsvd){
                      jobu  = "V";
                      jobvt = "V";
                  }

                  pdlacpy_( "All", &n, &n, 
                            Acpy, &i1, &i1, descAcpy, 
                            A, &i1, &i1, descA ); 

                  if (slsvd || polarsvd) {
                      my_elapsed_slsvd    = 0.0;
                      my_elapsed_slsvd    =- MPI_Wtime();
                      my_elapsed_polarsvd = 0.0;
                      my_elapsed_polarsvd =- MPI_Wtime();
                  }
                  pdgesvd_( jobu, jobvt, &n, &n, A, &i1, &i1, descA, 
                            S, 
                            U, &i1, &i1, descU, 
                            VT, &i1, &i1, descVT, 
                            Wloc, &lwork, &my_info_facto );
                  if (slsvd) {
                      my_elapsed_slsvd    += MPI_Wtime();
                      MPI_Allreduce( &my_elapsed_slsvd, &elapsed_slsvd, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                      sumtime_slsvd += elapsed_slsvd;
                      if ( elapsed_slsvd >= max_time_slsvd ) { max_time_slsvd = elapsed_slsvd;} 
                      if ( elapsed_slsvd <= min_time_slsvd ) { min_time_slsvd = elapsed_slsvd;} 
                  }
	          MPI_Allreduce( &my_info_facto, &info_facto, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "\nScaLAPACK dgesvd ends\n");

                  if (polarsvd && !slsvd){
                      /*
                       */ 
                       alpha = 1.0; beta = 0.0;
                       pdgemm_( "N", "N", &n, &n, &n, 
                                &alpha, 
                                U   , &i1, &i1, descU, 
                                VT,   &i1, &i1, descVT, 
                                &beta, 
                                A, &i1, &i1, descA);
                       alpha = 1.0; beta = 0.0;
                       pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);
                       for (i = 1; i <= n; i++) {
                           int idum1, idum2, iloc, jloc;
                           if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                           &&     ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                    iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                    jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                    Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                           }
                       } 
                       pdgemm_( "T", "N", &n, &n, &n, 
                                &alpha, 
                                VT   , &i1, &i1, descVT, 
                                Sigma,   &i1, &i1, descSigma, 
                                &beta, 
                                U, &i1, &i1, descU);
                       pdgemm_( "N", "N", &n, &n, &n, 
                                &alpha, 
                                U   , &i1, &i1, descU, 
                                VT,   &i1, &i1, descVT, 
                                &beta, 
                                Sigma, &i1, &i1, descSigma);
                       if (polarsvd) {
                           my_elapsed_polarsvd    += MPI_Wtime();
                           MPI_Allreduce( &my_elapsed_polarsvd, &elapsed_polarsvd, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                           sumtime_polarsvd += elapsed_polarsvd;
                           if ( elapsed_polarsvd >= max_time_polarsvd ) { max_time_polarsvd = elapsed_polarsvd;} 
                           if ( elapsed_polarsvd <= min_time_polarsvd ) { min_time_polarsvd = elapsed_polarsvd;} 
                       }
	               MPI_Allreduce( &my_info_facto, &info_facto, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                  }
                  /*
                   * Checking the polar factorization
                   */

                  if(polarsvd && check && !slsvd && iter == 0){
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
                        orth_Usvd  = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc)/frobA;

                       /*
                        * checking the factorization |A-Up*H|
                        */ 
                        pdlacpy_( "A", &n, &n, Acpy, &i1, &i1, descAcpy, H, &i1, &i1, descH );
                        pdgemm_( "N", "N", &n, &n, &n, 
                                 &alpha, 
                                 A, &i1, &i1, descA, 
                                 Sigma, &i1, &i1, descSigma, 
                                 &beta, 
                                 H, &i1, &i1, descH);
                        berr_UHsvd  = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc)/frobA;
                  }

                  if (check && slsvd && !polarsvd && iter == 0) {
                     if ( verbose & myrank_mpi == 0) printf( "\n check Only on first iter %d \n", iter);
                     if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing starts...\n");
                     alpha = 1.0; beta = 0.0;
                     pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);

                     for (i = 1; i <= n; i++) {
                             int idum1, idum2, iloc, jloc;
                             if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                             &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                     iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                     jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                     Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                             }

                     } 
                     pdgemm_( "N", "N", &n, &n, &n, 
                              &alpha, 
                              U   , &i1, &i1, descU, 
                              Sigma, &i1, &i1, descSigma, 
                              &beta, 
                              A, &i1, &i1, descA);
                     pdlacpy_( "All", &n, &n, 
                              Acpy, &i1, &i1, descAcpy, 
                              Wglo, &i1, &i1, descWglo ); 
                     beta = -1.0;
                     pdgemm_( "N", "N", &n, &n, &n, 
                              &alpha, 
                              A, &i1, &i1, descA, 
                              VT, &i1, &i1, descVT, 
                              &beta, 
                              Wglo, &i1, &i1, descWglo);
                     my_berr = pdlange_ ( "f", &n, &n, Wglo, &i1, &i1, descWglo, Wloc) / (frobA * n);

                        /* Accuracy of singular values */
                        for(i=0; i < n ; i++ )
                            D[i] = fabs(D[i]);
                        dlasrt_( "D", &n, D, &info );
                        dlasrt_( "D", &n, S, &info );
                        for(i=0; i < n ; i++ )
                            S[i] = S[i] - D[i];
                        for (i = 1; i <= n; i++) {
                             int idum1, idum2, iloc, jloc;
                             if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                             &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                     iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                     jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                     Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                             }
                        } 
                        my_acc_slsvd = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc) / norm_sv;

                        /* Orthogonality of Left singular vectors */
                        alpha = 0.0; beta = 1.0;
                        pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                        alpha = 1.0; beta = -1.0;
                        pdgemm_( "T", "N", &n, &n, &n, 
                                 &alpha, 
                                 U, &i1, &i1, descU, 
                                 U, &i1, &i1, descU, 
                                 &beta, 
                                 Sigma, &i1, &i1, descSigma);
                        my_orthL_slsvd = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                        /* Orthogonality of Right singular vectors */
                        alpha = 0.0; beta = 1.0;
                        pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                        alpha = 1.0; beta = -1.0;
                        pdgemm_( "N", "T", &n, &n, &n, 
                                 &alpha, 
                                 VT, &i1, &i1, descVT, 
                                 VT, &i1, &i1, descVT, 
                                 &beta, 
                                 Sigma, &i1, &i1, descSigma);
                        my_orthR_slsvd = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                        if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing ends...\n");
                  } 
                  pdlacpy_( "All", &n, &n, 
                            Acpy, &i1, &i1, descAcpy, 
                            A, &i1, &i1, descA ); 
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy back to A done\n");
               }
	       free( Wloc );

            }
            if ( polarsvd && myrank_mpi == 0) {
                  fprintf(stderr, "# Polar decomposition using ScaLAPACK DGESVD \n"); 
                  fprintf(stderr, "#\n");
	          fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ   \tAvg-Time     \tMax-Time    \tMin-Time    \tBerr_UpH  \tOrth_Up  \tinfo     \n");
	          fprintf(stderr, "   %6d \t%4d \t%4d \t%3d \t%3d ", n, nb, nprocs_mpi, nprow, npcol);
	          fprintf(stderr, "\t%6.2f \t\t%6.2f \t\t%6.2f \t\t%2.4e \t%2.4e \t%d \n", sumtime_polarsvd/niter, max_time_polarsvd,  min_time_polarsvd, berr_UHsvd, orth_Usvd,info_facto);
                  fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            }
	    if (slsvd && myrank_mpi == 0){
                  fprintf(stderr, "# ScaLAPACK DGESVD\n"); 
                  fprintf(stderr, "#\n");
	          fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ   \tAvg-Time     \tMax-Time    \tMin-Time  \tinfo     \tAcc-sv    \tOrth-Rsv    \tOrth-Lsv     \tBerr \n");
	          fprintf(stderr, "   %6d \t%4d \t%3d \t%3d \t%3d", n, nb, nprocs_mpi, nprow, npcol);
                  fprintf(stderr, "\t%6.2f  \t%6.2f \t\t%6.2f \t\t%d \t\t%2.4e \t%2.4e \t%2.4e \t%2.4e\n", sumtime_slsvd/niter, max_time_slsvd, min_time_slsvd, info_facto, my_acc_slsvd, my_orthR_slsvd, my_orthL_slsvd, my_berr);
            }

            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
            if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

	    free( Sigma );
	    free( A );
	    free( Acpy );
	    free( H );
	    free( Wglo );
	    free( tau );
	    free( D );

	    free( U );
	    free( VT );
	    free( S );
            if (verbose & myrank_mpi == 0) fprintf(stderr, "Free matrices done\n");
        } // End loop over range


        if (verbose & myrank_mpi == 0) fprintf(stderr, "Range loop ends\n");

	blacs_gridexit_( &i0 );
	MPI_Finalize();
        if (verbose & myrank_mpi == 0) fprintf(stderr, "Program ends...\n");
	return 0;
}



