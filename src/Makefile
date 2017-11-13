#/**
# *
# * Copyright (c) 2017, King Abdullah University of Science and Technology
# * All rights reserved.
# *
# **/
#
#/**
# *
# * @file Makefile
# *
# *  QDWH is a high performance software framework for computing 
# *  the polar decompositon on distributed-memory manycore systems provided by KAUST
# *
# * @version 2.0.0
# * @author Dalal Sukkari
# * @author Hatem Ltaief
# * @date 2017-11-13
# *
# **/


ELPADIR = $(HOME)/codes/elpa-2015.11.001/install
SCALAPACKDIR = $(HOME)/codes/scalapack_installer-intel/install
PLASMADIR = $(HOME)/codes/plasma/plasma-installer_2.8.0-intel/install

INC := -I../include/ $(INC)

MPICC = mpicc -cc=icc
MPIF77 = mpif90 -fc=ifort
CFLAGS = -Wall -I$(ELPADIR)/include/elpa-2015.11.001 -I$(PLASMADIR)/include -qopenmp -mkl=parallel -DHAVE_GETOPT_LONG

CC = $(MPICC)
LD = mpif90 -fc=ifort
#LIBS = -nofor_main -mkl=parallel -L$(PLASMADIR)/lib -llapacke -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -L$(SCALAPACKDIR)/lib -lsltmg -lmkl_blacs_intelmpi_lp64 -lpthread -lm 
LIBS = -nofor_main -mkl=parallel -L$(PLASMADIR)/lib -llapacke -L$(ELPADIR)/lib -lelpa -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -L$(SCALAPACKDIR)/lib -lsltmg -lmkl_blacs_intelmpi_lp64 -lpthread -lm 

TARGETS = main

all: $(TARGETS)

main:  pdgeqdwh.o pdgeqsvd.o pdgenm2.o
#main_svd:  test-blacs.o
	$(LD) -o $@ $^ $(LIBS)
.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@


clean:
	rm -f *.o main $(TARGETS)
