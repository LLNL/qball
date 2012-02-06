#-------------------------------------------------------------------------------
#
#  linux-pc_mpi.mk
#
#-------------------------------------------------------------------------------
# $Id: linux-pc.mk,v 1.10 2010/05/12 20:05:25 draeger1 Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
# GCCDIR=/usr/apps/gcc/3.1
# GCCDIR=/usr/apps/gcc/3.2.1
 GCCDIR=/usr
 MPIDIR=/usr/apps/mpich/1.2.7p1
 XERCESCDIR=$(HOME)/software/xml/xerces-c-src_2_5_0
 FFTWDIR=$(HOME)/software/fftw/fftw-linux/fftw-2.1.3/fftw
 BLASDIR=/usr/lib

 XERCESLIB=$(XERCESCDIR)/lib/libxerces-c.a
 FFTWLIB=$(FFTWDIR)/libfftw.a
 
# CXX=/usr/apps/mpich/1.2.4/bin/mpiCC
# CXX=$(GCCDIR)/bin/g++

# CXX=g++
 CXX=/usr/apps/mpich/1.2.7p1/bin/mpicxx
 LD=$(CXX)

# PLTFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
#             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
#             -DAPP_NO_THREADS -DXML_USE_NO_THREADS
 PLTFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DSIMPLECHKPT

# INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR)
 
 CXXFLAGS= -O3 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 
# CXXFLAGS= -g -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 

# LIBPATH = -L$(FFTWDIR) -L/usr/X11R6/lib \
#           -L$(MPIDIR)/lib -L $(BLASDIR) -L $(GCCDIR)/lib -L$(XERCESCDIR)/lib
 LIBPATH = -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L $(BLASDIR) -L $(GCCDIR)/lib
  
# LIBS =  $(PLIBS) -lfftw -llapack -lblas -lm -lmpich -lpmpich -lmpich \
#         -lg2c -lxerces-c
 LIBS =  $(PLIBS) -lfftw -llapack -lblas -lm -lmpich -lpmpich -lmpich \
         -lg2c
 
 LDFLAGS = $(LIBPATH) $(LIBS) 

 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/software/blacs/blacs-linux/BLACS/LIB
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/software/scalapack/scalapack-linux/SCALAPACK
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a
 TOOLSLIB      = $(SCALAPACK_DIR)/tools_$(PLT).a
 REDISTLIB     = $(SCALAPACK_DIR)/redist_$(PLT).a

 LAPACKLIB = -llapack
 BLASLIB = -lblas

 # Parallel libraries
# PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
