POLAR
================
POLAR is a high performance open-source software package to compute the polar decomposition ((PD)) of a dense matrix A = UH
based on QR-based Dynamically Weighted Halley (QDWH) and ZOLO-PD algorithms.

QDWH
================

The **QR-based Dynamically Weighted Halley** (QDWH) is one of
the most popular algorithms to compute the polar decomposition. It is backward stable and converges
in at most six iterations. 

ZOLO-PD
================

ZOLO-PD relies on the **Zolotarev function which is the best rational approximant to the sign function**  
for computing the polar decomposition.
Building upon on the QR-based
Dynamically Weighted Halley (QDWH) algorithm, the key idea
lies in finding the best rational approximation for the scalar sign function,
which also corresponds to the polar factor for symmetric matrices,
to further accelerate the QDWH convergence.
Based on the Zolotarev rational functions---introduced by Zolotarev (ZOLO) in
1877--- this new PD algorithm ZOLO-PD converges within two iterations even for ill-conditioned matrices,
instead of the original six iterations needed for QDWH.
ZOLO-PD uses the property of Zolotarev functions that optimality is maintained when
two functions are composed in an appropriate manner.
The resulting ZOLO-PD has a convergence rate up to seventeen,
in contrast to the cubic convergence rate for QDWH.
This comes at the price of higher arithmetic costs and memory footprint. These
extra floating-point operations can, however, be processed in an
embarrassingly parallel fashion. 

Current Features of QDWH/ZOLOPD
===========================

- Written in C.
- Support for Double Precision.
- Support for Two-Dimensional Block Cyclic Data Distribution.
- ScaLAPACK Interface / Native Interface.
- ScaLAPACK-Compliant Error Handling.
- ScaLAPACK-Derived Testing Suite
- ScaLAPACK-Compliant Accuracy.
 
Programming models and backends:
1.  MPI
2.  ScaLAPACK


Installation
============

The installation requires at least **CMake** of version 3.2.3. To build the polar decomposition based on QDWH/ZOLOPD,
follow these instructions:

1.  Get polar from git repository

        git clone git@github.com:ecrc/polar

2.  Go into polar folder

        cd polar

3.  Create build directory and go there

        mkdir build && cd build

4.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/ 

5.  To use exist dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/ -DSCALAPACK_DIR=/path/to/scalapack/install/ -DSLTMG_LIBRARIES=/path/to/scalapack/install/lib/libsltmg.a

5.  To build the testing binaries (optional) add the following:

        -DPOLAR_TESTING:BOOL=ON

5.  Build polar

        make -j

6.  Install polar

        make install

7. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
polar based on QDWH/ZOLO-PD.

Testing and Timing
==================

The directories testing and timing contain an example 
to test the accuracy and the performance of QDWH/ZOLO-PD using
ill (with condition number less than 5.e12) and well-conditioned random matrices.

   The complete list of options is available below with -h option:
  
  ```
       "======= QDWH/ZOLOPD testing using ScaLAPACK\n"
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
2.  Extend task-based programming model
3.  Port to various dynamic runtime systems


References
==========
1. H. Ltaief, D. Sukkari, A. Esposito, Y. Nakatsukasa and D. Keyes, Massively Parallel 
Polar Decomposition on Distributed-Memory Systems, *Submitted to IEEE Transactions on 
Parallel Computing TOPC*, http://hdl.handle.net/10754/626359.1, 2018.
2. D. Sukkari, H. Ltaief, A. Esposito and D. Keyes, A QDWH-Based SVD Software Framework on
Distributed-Memory Manycore Systems, *Submitted to ACM Transactions on Mathematical Software TOMS*, 
http://hdl.handle.net/10754/626212, 2017.
3. D. Sukkari, H. Ltaief, M. Faverge, and D. Keyes, Asynchronous Task-Based Polar
Decomposition on Massively Parallel Systems, *IEEE Transactions on Parallel and 
Distributed Systems TPDS*, volume 29, pages 312–323, https://ieeexplore.ieee.org/document/8053812/, 2017.
4. D. Sukkari, H. Ltaief and D. Keyes, A High Performance QDWH-SVD Solver using
Hardware Accelerators, *ACM Transactions on Mathematical Software TOMS*, vol. 43 (1), pp. 1-25, 2016.
5. D. Sukkari, H. Ltaief and D. Keyes, High Performance Polar Decomposition for SVD
Solvers on Distributed Memory Systems, Best Papers, *Proceedings of the 22nd International 
Euro-Par Conference*, https://doi.org/10.1007/978-3-319-43659-3_44, 2016.
6. D.Sukkari, H. Ltaief and D. Keyes, A High Performance QDWH-SVD Solver using 
Hardware Accelerators, *ACM Transactions on Mathematical Software TOMS*, 
http://doi.acm. org/10.1145/2894747, volume 43, pages 6:1–6:25, 2016.
7. Y. Nakatsukasa and N. J. Higham, Stable and Efficient Spectral Divide and Conquer 
Algorithms for the Symmetric Eigenvalue Decomposition and the SVD, *SIAM Journal on Scientific Computing*,
vol. 35, no. 3, pp. A1325–A1349, http://epubs.siam.org/doi/abs/10.1137/120876605, 2013.
8. Y. Nakatsukasa, R. Freund, using Zolotarev's Rational Approximation for Computing the Polar, 
Symmetric Eigenvalue, and Singular Value Decompositions, *SIAM Review*, 
https://books.google.com.sa/books?id=a9d7rgEACAAJ, 2016.


Questions?
==========
Please feel free to create an issue on Github for any questions and inquiries.

