#-------------------------------------------------------------------------------
#
#  linux-pc.mk
#
#-------------------------------------------------------------------------------
# $Id: linux-pc.mk,v 1.1.1.1 2005/08/18 17:23:33 draeger1 Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# GCCDIR=/usr/apps/gcc/3.1

 CXX=g++
 LD=$(CXX)

# FFTWDIR=$(HOME)/fftw/linux-pc/fftw-1.3/src
 BLASDIR=/usr/lib

# INCLUDE = -I$(FFTWDIR)
  
 CXXFLAGS= -O2 -D$(PLT) -DADD_ $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(BLASDIR) -L$(FFTWDIR) -L$(GCCDIR)/lib

# LIBS = -lfftw -llapack -lblas -lm -lg2c

 LIBS = -lfftw -llapack -lblas -lm
 
 LDFLAGS = $(LIBPATH) $(LIBS)
#-------------------------------------------------------------------------------
