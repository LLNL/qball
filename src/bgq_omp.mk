#-------------------------------------------------------------------------------
#
#  BlueGene/Q (IBM Q32)
#
#-------------------------------------------------------------------------------
# $Id: linux_x86_64_intel.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=BGQ
#-------------------------------------------------------------------------------

 LIBHOME = $(HOME)/software

# XERCESCDIR=$(LIBHOME)/xml/xml-bgq/xerces-c-src_2_5_0
# XERCESCLIBDIR=$(XERCESCDIR)/lib
# XERCESLIB=$(XERCESCLIBDIR)/libxerces-c.a
 FFTWDIR=$(LIBHOME)/fftw/fftw-bgq/fftw-2.1.3/fftw
 FFTWLIB=$(FFTWDIR)/libfftw.a
 BLASDIR=$(LIBHOME)/blas/blas-bgq-xlc
 LAPACKDIR=$(LIBHOME)/lapack/lapack-bgq-xlc
 SCALAPACK_DIR = $(LIBHOME)/scalapack-2.0/scalapack-bgq-xlc
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 ESSLDIR = /usr/local/tools/essl/5.1

 CXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r
 LD=$(CXX)

 DFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 
 INCLUDE = -I$(FFTWDIR) -I$(ESSLDIR)/include
 
# CXXFLAGS= -g -O3 -qarch=qp -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)
 CXXFLAGS= -g -O3 -qsmp=omp -qarch=qp -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

# LIBPATH = -L$(FFTWDIR) -L$(LAPACKDIR) -L$(BLASDIR) 
# LIBS =  $(SCALAPACKLIB) -lfftw -lblas -llapack -lm
# LIBPATH = -L$(FFTWDIR) -L$(LAPACKDIR) -L$(BLASDIR) -L$(ESSLDIR)/lib -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -L/opt/ibmcmp/xlf/bg/14.1/bglib64
# LIBS =  $(SCALAPACKLIB) -lfftw -lesslsmpbg -lblas -llapack -lxlf90_r -lxlsmp -lxlfmath
# LDFLAGS = $(LIBPATH) $(LIBS) -qarch=qp -lc -lnss_files -lnss_dns -lresolv

 LIBPATH = -L$(FFTWDIR) -L$(LAPACKDIR) -L$(BLASDIR) -L$(ESSLDIR)/lib -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -L/opt/ibmcmp/xlmass/bg/7.3/bglib64 -L/opt/ibmcmp/xlf/bg/14.1/bglib64
 LIBS =  $(SCALAPACKLIB) -lfftw -lesslbg -lblas -llapack -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath
 LDFLAGS = $(LIBPATH) $(LIBS) -qarch=qp -lc -lnss_files -lnss_dns -lresolv


#-------------------------------------------------------------------------------
