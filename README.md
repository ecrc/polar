The  main  purpose of QDWH is to compute the polar decomposition A=UH using ScaLAPACK numerical library.
ScaLAPACK-QDWH is written in C and requires ScaLAPACK installation, as the main software dependency.
The main driver compares the ScaLAPACK performance of QDWH against the SVD-based polar decomposition
(using PDGESVD).
1. How to run the code
  1. To run the test on N nodes, the number of tasks (nT) = N * (number_of_cores per node ). The programming model is pure MPI (no OpenMP, i.e., sequential BLAS).
  2. PxQ is the process grid configuration, where (nT - PxQ = 0)
  3. To find the polar decomposition:
    1. Using ScaLAPACK-QDWH (--polarqdwh)
    2. Using ScaLAPACK-SVD-based (--polarsvd)
    3. On Cray systems, the launching command typically looks like:
    
       srun --ntasks=nT --hint=nomultithread ./main --nprow p --npcol q  --b 64 --cond 1e16 --niter 1 --n_range start:stop:step  --check --polarqdwh --polarsvd
       
       This runs QDWH-based and SVD-based polar decomposition and will check on the numerical accuracy, for one single iteration, with ill-conditioned matrix, for PxQ processor grid size, block size 64, no multithreading, on nT number of MPI processes, across a range of matrix sizes defined by start:stop:step.
  4. To Find the SVD using QDWH:
       Using ScaLAPACK-QDWH: the calculation of the polar decomposition is followed by MRRR (--qwmr) or DC (--qwdc) as the symmetric eigensolver.
       
       Using ScaLAPCK-PDGESVD: --slsvd
       
       srun --ntasks=nT --hint=nomultithread ./main --nprow p --npcol q  --b 64 --cond 1e16 --niter 1 --n_range start:stop:step  --check --qwmr --qwdc --slsvd

  5. To use LU to estimate the condition number of the matrix: add the option --optcond 0 (default:optcond 1 to use qr)
  6. The complete list of options is available below with -h option:
  
  ```
       "======= QDWHsvd timing using ScaLAPACK"
       " -p      --nprow         : Number of MPI process rows"
       " -q      --npcol         : Number of MPI process cols"
       " -jl     --lvec          : Compute left singular vectors"
       " -jr     --rvec          : Compute right singular vectors"
       " -n      --N             : Dimension of the matrix"
       " -b      --nb            : Block size"
       " -m      --mode          : [1:6] Mode from pdlatms used to generate the matrix"
       " -k      --cond          : Condition number used to generate the matrix"
       " -o      --optcond       : Estimate Condition number using QR"
       " -i      --niter         : Number of iterations"
       " -r      --n_range       : Range for matrix sizes Start:Stop:Step"
       " -polarqdwh --polarqdwh  : Find polar decomposition using QDWH A=UH "
       " -polarsvd  --polarsvd   : Find the polar decomposition using scalapack-svd "
       " -s      --slsvd         : Run reference ScaLAPACK SVD"
       " -w      --qwmr          : Run QDWH SVD with ScaLAPACK MRRR EIG"
       " -e      --qwdc          : Run QDWH SVD with ScaLAPACK DC EIG"
       " -c      --check         : Check the solution"
       " -fqwsvd --profqwsvd     : Enable profiling QDWHsvd"
       " -fqw    --profqw        : Enable profiling QDWH"
       " -v      --verbose       : Verbose"
       " -h      --help          : Print this help"
```
