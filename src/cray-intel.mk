 PLT=CRAY

# modules to load before building:
#
# module swap PrgEnv-pgi PrgEnv-intel
# module load fftw/2.1.5.7
#   to get performance information:
# module load perftools-lite 

 LIBHOME = $(HOME)/software
 LIBSCIDIR = $(CRAY_LIBSCI_DIR)/intel/140/haswell/lib/libsci_intel.a
 SCALAPACK_DIR = $(LIBHOME)/scalapack-cray-intel
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# build ScaLAPACK with -DAdd_

 #XERCESCDIR=$(HOME)/software/xml/xerces-c-3.1.1-bgq/src
 #XERCESCLIBDIR=$(XERCESCDIR)/.libs
 #XERCESLIB=$(XERCESCLIBDIR)/libxerces-c.a

# execute module load PrgEnv-intel before building
 CXX=CC
 CC=cc
 AR=ar
 RANLIB=ranlib

 LD=$(CXX)

 DFLAGS += -DUSE_FFTW2 -DUSE_DFFTW -DADD_
 
 INCLUDE = -I$(FFTW_INC)
 
 CXXFLAGS= -O3 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
 CFLAGS= -O3 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBS =  $(SCALAPACKLIB) $(LIBSCIDIR) $(FFTW_DIR)/libdfftw.a
 LDFLAGS = $(LIBS)
