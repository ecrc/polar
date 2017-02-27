#/**
# QDWH-SVD
# *
# * (C) Copyright 2016 King Abdullah University of Science and Technology
# * Authors:
# * Dalal Sukkari (dalal.sukkari@kaust.edu.sa)
# * David Keyes (david.keyes@kaust.edu.sa)
# * Hatem Ltaief (hatem.ltaief@kaust.edu.sa)
# *  
# * Redistribution  and  use  in  source and binary forms, with or without
# * modification,  are  permitted  provided  that the following conditions
# * are met:
# * 
# * Redistributions  of  source  code  must  retain  the above copyright
# * notice,  this  list  of  conditions  and  the  following  disclaimer.
# * Redistributions  in  binary  form must reproduce the above copyright
# * notice,  this list of conditions and the following disclaimer in the
# * documentation  and/or other materials provided with the distribution.
# * Neither  the  name of the King Abdullah University of Science and
# * Technology nor the names of its contributors may be used to endorse
# * or promote products derived from this software without specific prior
# * written permission.
# *
# *
# THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
# LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

TARGETS = testing_pdgesvd testing_pdgeqdwh

all: $(TARGETS)

testing_pdgesvd testing_pdgeqdwh:  testing_pdgesvd.o ../src/pdgeqsvd.o ../src/pdgeqdwh.o ../src/pdgenm2.o
	$(LD) -o $@ $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

clean:
	rm -f *.o testing_pdgesvd testing_pdgeqdwh
