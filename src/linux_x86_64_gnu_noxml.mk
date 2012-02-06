#-------------------------------------------------------------------------------
#
#  linux_x86_64_gnu_noxml.mk
#
#-------------------------------------------------------------------------------
# $Id: linux_x86_64_intel.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=LINUX_X86_64_GNU
#-------------------------------------------------------------------------------
 FFTWDIR=$(HOME)/software/fftw/fftw-linux_x86_64/fftw-2.1.3/fftw
 FFTWLIB=$(FFTWDIR)/libfftw.a
# BLASDIR=$(HOME)/software/blas/blas-linux_x86-64
# LAPACKDIR=$(HOME)/software/lapack/lapack-linux_x86-64
 SCALAPACK_DIR = $(HOME)/software/scalapack-2.0/scalapack-linux-gnu_x86-64
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 BLASDIR=/usr/local/tools/mkl-10.3.1/lib

 CXX=mpig++
 LD=$(CXX)

 DFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 
 INCLUDE = -I$(FFTWDIR)
 
 CXXFLAGS= -g -O3 -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)

# LIBPATH = -L$(FFTWDIR) -L$(BLASDIR) -L$(LAPACKDIR)
# LIBS =  $(SCALAPACKLIB) -lgfortran -lfftw -lblas -llapack 
 LIBPATH = -L$(FFTWDIR) -L$(BLASDIR)
 LIBS =  $(SCALAPACKLIB) -lfftw -omp -openmp -lmkl_core -lmkl_gnu_thread -lmkl_gf_lp64

# LDFLAGS = $(LIBPATH) $(LIBS) -Wl,-rpath,/usr/local/tools/mkl-8.1.1.004/lib
 LDFLAGS = $(LIBPATH) $(LIBS) 

#-------------------------------------------------------------------------------
