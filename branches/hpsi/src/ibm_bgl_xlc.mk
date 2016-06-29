#-------------------------------------------------------------------------------
#
#  ibm_bgl_xlc.mk
#
#-------------------------------------------------------------------------------
# $Id: ibm_bgl_xlc.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=BGL
#-------------------------------------------------------------------------------
 SCALAPACK_DIR = $(HOME)/software/scalapack-2.0/scalapack-bgl-xlc
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 BLASDIR=/bgl/local/lib
# LAPACKLIB = -llapack440
# BLASLIB = -lblas440
 LAPACKLIB = $(HOME)/software/lapack/lapack-bgl-xlc/liblapack.a
 BLASLIB = $(HOME)/software/blas/blas-bgl-xlc/libblas.a
 DGEMMLIB = $(HOME)/bglsoftware/blas/libdgemm.rts.a
 ZGEMMLIB = $(HOME)/bglsoftware/blas/libzgemm_dc_rel3.rts.a
 FFTWDIR=$(HOME)/bglsoftware/bglfftwgel-2.1.5.pre5/fftw
 FFTWLIBDIR=$(FFTWDIR)/.libs
 FFTWLIB=$(FFTWLIBDIR)/libfftw.a

 CXX=mpxlC
 LD=$(CXX)

 DFLAGS += -DUSE_FFTW -DUSE_MPI -DSCALAPACK -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DMPICH_IGNORE_CXX_SEEK -D$(PLT)
 INCLUDE =  -I$(FFTWDIR) 
 CXXFLAGS= -g -O3 -qarch=440 $(INCLUDE) $(DFLAGS)

 XLCMATHLIBS = -L/opt/ibmcmp/xlf/bg/11.1/blrts_lib -lxlf90 -lxlopt -lxlomp_ser -lxl -lxlfmath  -lmassv

 LIBPATH = -L$(FFTWLIBDIR) -L$(BLASDIR)
 LIBS =  $(SCALAPACKLIB) -lfftw $(LAPACKLIB) $(DGEMMLIB) $(ZGEMMLIB) $(BLASLIB) -lg2c $(XLCMATHLIBS)

 LDFLAGS = $(LIBPATH) $(LIBS)

#-------------------------------------------------------------------------------
