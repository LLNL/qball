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
 SCALAPACK_DIR = $(LIBHOME)/scalapack-2.0/scalapack-bgq-xlc-jaggemm
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 ESSLDIR = /usr/local/tools/essl/5.1
 HPMLIBS = -L/usr/local/tools/mpitrace/lib -lmpihpm -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm
 JAGGEMMLIB = $(LIBHOME)/jaggemm_opt/libjaggemm.a

# CXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r
BGQ_SDK_PATH = /bgsys/drivers/ppcfloor
#BGQ_SDK_PATH=/bgsys/drivers/DRV2012_0229_1625/ppc64-rhel60
CXX=$(BGQ_SDK_PATH)/comm/xl/bin/mpixlcxx_r
CC=$(BGQ_SDK_PATH)/comm/xl/bin/mpixlc_r

 LD=$(CXX)

 DFLAGS += -DPRINTALL -DUSE_JAGGEMM -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DTAU
 
 INCLUDE = -I$(FFTWDIR) -I$(ESSLDIR)/include
 
# CXXFLAGS= -g -O3 -qarch=qp -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)
# CXXFLAGS= -g -O3 -qsmp=omp -qarch=qp -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
 CXXFLAGS= -g -O3 -qarch=qp -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
 CFLAGS= -qhot=novector -qsimd=auto -g -O3 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

# LIBPATH = -L$(FFTWDIR) -L$(LAPACKDIR) -L$(BLASDIR) 
# LIBS =  $(SCALAPACKLIB) -lfftw -lblas -llapack -lm
 LIBPATH = -L$(FFTWDIR) -L$(LAPACKDIR) -L$(BLASDIR) -L$(ESSLDIR)/lib -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -L/opt/ibmcmp/xlf/bg/14.1/bglib64
 LIBS =  $(SCALAPACKLIB) -lfftw $(JAGGEMMLIB) -lesslsmpbg -lblas -llapack -lxlf90_r -lxlsmp -lxlfmath $(HPMLIBS)
 LDFLAGS = $(LIBPATH) $(LIBS) -qarch=qp -lc -lnss_files -lnss_dns -lresolv

TAUROOTDIR = $(LIBHOME)/tau/tau-2.21.2
ifneq (,$(findstring DTAU,$(DFLAGS)))
        include  $(TAUROOTDIR)/include/Makefile
        CXXFLAGS+=$(TAU_INCLUDE) $(TAU_DEFS)
#       LIBS+=$(TAU_MPI_LIBS) $(TAU_LIBS)                                                                                  
        LDFLAGS+= $(TAU_LIBS)
endif


#-------------------------------------------------------------------------------
