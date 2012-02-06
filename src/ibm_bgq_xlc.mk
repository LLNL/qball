#-------------------------------------------------------------------------------
#
#  BlueGene/Q (IBM Q32)
#
#-------------------------------------------------------------------------------
# $Id: linux_x86_64_intel.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=BGQ
#-------------------------------------------------------------------------------

 LIBHOME = /bgusr/draeger/software

# XERCESCDIR=$(LIBHOME)/xml/xml-bgq/xerces-c-src_2_5_0
# XERCESCLIBDIR=$(XERCESCDIR)/lib
# XERCESLIB=$(XERCESCLIBDIR)/libxerces-c.a
 FFTWDIR=$(LIBHOME)/fftw/fftw-ibm_bgq/fftw-2.1.3/fftw
 FFTWLIB=$(FFTWDIR)/libfftw.a
 BLASDIR=$(LIBHOME)/blas
 LAPACKDIR=$(LIBHOME)/lapack
 SCALAPACK_DIR = $(LIBHOME)/scalapack-2.0/scalapack-ibm_bgq
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

 CXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r
 LD=$(CXX)

 DFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 
 INCLUDE = -I$(FFTWDIR)
 
 CXXFLAGS= -g -O3 -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(FFTWDIR) -L$(BLASDIR) -L$(LAPACKDIR)
 LIBS =  $(SCALAPACKLIB) -lfftw -lblas -llapack
 LDFLAGS = $(LIBPATH) $(LIBS)

#-------------------------------------------------------------------------------
