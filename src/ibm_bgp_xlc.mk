#-------------------------------------------------------------------------------
#
#  ibm_bgp_xlc.mk
#
#-------------------------------------------------------------------------------
# $Id: ibm_bgp_xlc.mk,v 1.1 2010/01/22 00:35:01 draeger1 Exp $
#
 PLT=BGP
#-------------------------------------------------------------------------------
 BLASDIR=/usr/local/tools/mkl-10.3.1/lib
 SCALAPACK_DIR = $(HOME)/software/scalapack-2.0/scalapack-bgp-xlc
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 ESSLDIR = /bgsys/ibm_essl/sles10/prod/opt/ibmmath

 CXX=mpixlcxx
 LD=$(CXX)

# DFLAGS += -DUSE_ESSL -DUSE_MPI -DSCALAPACK -DADD_ -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D$(PLT)
 DFLAGS += -DUSE_ESSL -DUSE_MPI -DSCALAPACK -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D$(PLT)
 
 INCLUDE =  -I$(ESSLDIR)/include
 
 CXXFLAGS= -g -O3 -qarch=450 -qtune=450 $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(ESSLDIR)/lib -L/opt/ibmcmp/xlsmp/bg/1.7/bglib

 XLCMATHLIBS = -L/opt/ibmcmp/xlf/bg/11.1/lib -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath  -lm
 BLASLIB = $(HOME)/dawn/software/blas/libblas.a
 LAPACKLIB = $(HOME)/dawn/software/blas/liblapack_bgp.a

 LIBS =  $(SCALAPACKLIB) -lesslbg $(LAPACKLIB) $(BLASLIB) $(XLCMATHLIBS) -lrt -lpthread
 LDFLAGS = $(LIBPATH) $(LIBS)

#-------------------------------------------------------------------------------
