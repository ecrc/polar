QDWH
================

The **QR-based Dynamically Weighted Halley** (QDWH) package is a high performance open-source software
for computing the polar decomposition of a dense matrix A = UH. 
QDWH is written in C and requires ScaLAPACK installation, as the main software dependency.
QDWH provides the polar decomposition and the performance and the accuracy of the results. 
QDWH currently supports double precision arithmetics and run on shared and distributed-memory systems,
using MPI.

Current Features of QDWH
===========================

- Support double precision.
- Support dense two-dimensional block cyclic data distribution.
- ScaLAPACK Interface / Native Interface.
- ScaLAPACK-Compliant Error Handling.
- ScaLAPACK-Derived Testing Suite
- ScaLAPACK-Compliant Accuracy.
 
Programming models (backends):
1.  MPI
2.  ScaLAPACK


Installation
============

Installation requires at least **CMake** of version 3.2.3. To build QDWH,
follow these instructions:

1.  Get QDWH from git repository

        git clone git@github.com:ecrc/qdwh

2.  Go into QDWH folder

        cd qdwh

3.  Create build directory and go there

        mkdir build && cd build

4.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/ 

5.  To build the testing binaries (optional)

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/ -DQDWH_TESTING:BOOL=ON

5.  Build QDWH

        make -j

6.  Install QDWH

        make install

7. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
QDWH.

Testing and Timing
==================

The directories testing and timing contain an example 
to test the accuracy and the performance of QDWH using
ill/well-conditioned random matrices.

   The complete list of options is available below with -h option:
  
  ```
       "======= QDWH testing using ScaLAPACK\n"
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
```
     On Cray systems, the launching command typically looks like:
    
       srun --ntasks=nT --hint=nomultithread ./main --nprow p --npcol q  --b 64 --cond 1e16 --niter 1 --n_range start:stop:step  --check

     1. The number of the nodes is N, the number of tasks (nT) = N * (number_of_cores per node ). The programming model is pure MPI (no OpenMP, i.e., sequential BLAS).
     2. PxQ is the process grid configuration, where (nT - PxQ = 0)


TODO List
=========

1.  Add support for the other precisions 
2.  Provide full StarPU support (GPUs and distributed-memory systems)
3.  Port to other dynamic runtime systems


References
==========
1. D. Sukkari, H. Ltaief, A. Esposito and D. Keyes, A QDWH-Based SVD Software Framework on
Distributed-Memory Manycore Systems, *Submitted to ACM Transactions on Mathematical Software*, 
http://hdl.handle.net/10754/626212, 2017.
2. D. Sukkari, H. Ltaief, M. Faverge, and D. Keyes, Asynchronous Task-Based Polar
Decomposition on Massively Parallel Systems, *IEEE Transactions on Parallel and 
Distributed Systems*, 2017.
3. D. Sukkari, H. Ltaief and D. Keyes, A High Performance QDWH-SVD Solver using
Hardware Accelerators, *ACM Transactions on Mathematical Software*, vol. 43 (1), pp. 1-25, 2016.
4. D. Sukkari, H. Ltaief and D. Keyes, High Performance Polar Decomposition for SVD
Solvers on Distributed Memory Systems, Best Papers, Proceedings of the 22nd International Euro-Par Conference, 2016.
5. Y. Nakatsukasa and N. J. Higham, Stable and Efficient Spectral Divide and Conquer 
Algorithms for the Symmetric Eigenvalue Decomposition and the SVD, *SIAM Journal on Scientific Computing*,
vol. 35, no. 3, pp. A1325–A1349, http://epubs.siam.org/doi/abs/10.1137/120876605, 2013.


Questions?
==========
Please feel free to create an issue on Github for any questions and inquiries.

