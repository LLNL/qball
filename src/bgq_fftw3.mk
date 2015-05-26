#-------------------------------------------------------------------------------
#
#  BlueGene/Q (IBM Q32)
#
#-------------------------------------------------------------------------------
# $Id: linux_x86_64_intel.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=BGQ
#-------------------------------------------------------------------------------

# do this to force manual compilation of key threaded objects
 OMPHACK = 1

 LIBHOME = $(HOME)/software
 BLASDIR=$(LIBHOME)/blas/blas-bgq-xlc
 LAPACKDIR=$(LIBHOME)/lapack/lapack-bgq-xlc
 SCALAPACK_DIR = $(LIBHOME)/scalapack-2.0/scalapack-bgq-xlc-jaggemm
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 ESSLDIR = /usr/local/tools/essl/5.1
 HPMLIBS = -L/usr/local/tools/mpitrace/lib -lmpihpm_smp -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm
 JAGGEMMLIB = $(LIBHOME)/jaggemm_opt/libjaggemm.a
 XERCESCDIR=$(HOME)/software/xml/xerces-c-3.1.1-bgq/src
 XERCESCLIBDIR=$(XERCESCDIR)/.libs
 XERCESLIB=$(XERCESCLIBDIR)/libxerces-c.a
# FFTW3DIR = /usr/local/tools/fftw-3.3.3
# FFTWLIB=$(FFTW3DIR)/lib/libfftw3.a
 FFTW3DIR = /usr/local/tools/fftw-3.3.3
 FFTWLIB=$(FFTW3DIR)/lib/libfftw3.a $(FFTW3DIR)/lib/libfftw3_threads.a

 BGQ_SDK_PATH = /bgsys/drivers/ppcfloor
 CXX=$(BGQ_SDK_PATH)/comm/xl/bin/mpixlcxx_r
 CC=$(BGQ_SDK_PATH)/comm/xl/bin/mpixlc_r
 AR=$(BGQ_SDK_PATH)/gnu-linux/powerpc64-bgq-linux/bin/ar
 RANLIB=$(BGQ_SDK_PATH)/gnu-linux/powerpc64-bgq-linux/bin/ranlib

 LD=$(CXX)

 DFLAGS += -DPRINTALL -DUSE_FFTW3 -DUSE_FFTW3_THREADS -DUSE_JAGGEMM -DUSE_CSTDIO_LFS \
	-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DHPM -DUSE_XERCES -DXERCESC_3
 
 INCLUDE = -I$(FFTW3DIR)/include -I$(ESSLDIR)/include -I$(XERCESCDIR)
 
 CXXFLAGS= -g -O3 -qarch=qp -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
 CFLAGS= -qhot=novector -qsimd=auto -g -O3 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(LAPACKDIR) -L$(BLASDIR) -L$(ESSLDIR)/lib -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 \
	-L/opt/ibmcmp/xlf/bg/14.1/bglib64 -L$(XERCESCLIBDIR)
 LIBS =  $(SCALAPACKLIB) $(JAGGEMMLIB) $(FFTWLIB) -lesslsmpbg -lblas -llapack -lxlf90_r -lxlsmp -lxlfmath $(HPMLIBS) -lxerces-c
 LDFLAGS = $(LIBPATH) $(LIBS) -qarch=qp -lc -lnss_files -lnss_dns -lresolv

#-------------------------------------------------------------------------------
