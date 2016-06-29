#-------------------------------------------------------------------------------
#
#  zeus-chaos4.mk
#
#-------------------------------------------------------------------------------
# $Id: linux_x86_64_intel.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 XERCESCDIR=$(HOME)/software/xml/xml-zeus/xerces-c-src_2_5_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib
 XERCESLIB=$(XERCESCLIBDIR)/libxerces-c.a
 FFTWDIR=$(HOME)/software/fftw/fftw-linux_x86_64/fftw-2.1.3/fftw
 FFTWLIB=$(FFTWDIR)/libfftw.a
 BLASDIR=/usr/local/tools/mkl/lib

 CXX=mpig++
 LD=$(CXX)

 DFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DXML -DUSE_XERCES
 
 INCLUDE = -I$(FFTWDIR) -I$(XERCESCDIR)/include
 
 CXXFLAGS= -g -O3 -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(FFTWDIR) -L$(BLASDIR) -L$(XERCESCLIBDIR)
 LIBS =  $(PLIBS) -lgfortran -lfftw  -openmp -lmkl -lmkl_lapack64 -lxerces-c

 LDFLAGS = $(LIBPATH) $(LIBS) -Wl,-rpath,/usr/local/tools/mkl-8.1.1.004/lib

 SCALAPACK_DIR = $(HOME)/software/scalapack-2.0/scalapack-gnu-linux_x86-64
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 PLIBS = $(SCALAPACKLIB)

#-------------------------------------------------------------------------------
