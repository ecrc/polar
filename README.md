
 The  main  purpose of QDWH is to compute the polar decomposition A=UH using ScaLAPACK numerical library.
 ScaLAPACK-QDWH is written in C and requires ScaLAPACK installation, as the main software dependency.
 The main driver compares the ScaLAPACK performance of QDWH against the SVD-based polar decomposition
 (using PDGESVD).

1. How to run the code:
   1.1 To run the test on N nodes, the number of tasks (nT) = N * (number_of_cores per node ). The programming model is pure MPI (no OpenMP, i.e., sequential BLAS).
   1.2 PxQ is the process grid configuration, where (nT - PxQ = 0)
   1.3 To find the polar decomposition: 
       a) Using ScaLAPACK-QDWH (--polarqdwh)
       b) Using ScaLAPACK-SVD-based (--polarsvd)
       c) On Cray systems, the launching command typically looks like: 
       srun --ntasks=nT --hint=nomultithread ./main --nprow p --npcol q  --b 64 --cond 1e16 --niter 1 --n_range start:stop:step  --check --polarqdwh --polarsvd
       This runs QDWH-based and SVD-based polar decomposition and will check on the numerical accuracy, for one single iteration, with ill-conditioned matrix, for PxQ processor grid size, block size 64, 
       no multithreading, on nT number of MPI processes, across a range of matrix sizes defined by start:stop:step.
   1.4 To Find the SVD using QDWH: 
       Using ScaLAPACK-QDWH: the calculation of the polar decomposition is followed by MRRR (--qwmr) or DC (--qwdc) as the symmetric eigensolver.
       Using ScaLAPCK-PDGESVD: --slsvd
       srun --ntasks=nT --hint=nomultithread ./main --nprow p --npcol q  --b 64 --cond 1e16 --niter 1 --n_range start:stop:step  --check --qwmr --qwdc --slsvd 
   1.5 To use LU to estimate the condition number of the matrix: add the option --optcond 0 (default:optcond 1 to use qr)

   1.6 The complete list of options is available below with -h option:
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
       " -c      --check         : Check the solution\n"
       " -fqwsvd --profqwsvd     : Enable profiling QDWHsvd\n"
       " -fqw    --profqw        : Enable profiling QDWH\n"
       " -v      --verbose       : Verbose\n"
       " -h      --help          : Print this help\n" );

