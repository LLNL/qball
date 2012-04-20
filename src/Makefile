#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.46 2010/08/26 17:44:16 draeger1 Exp $
#------------------------------------------------------------------------------
#
# Get name of architecture-specific include file from hostname
include Makefile.arch

# If defined, use architecture file set with ARCH variable
ARCH ?= $(ARCHGUESS)
include $(ARCH).mk
#------------------------------------------------------------------------------
#
TESTEXECS=testMatrix testFTGrid testGridFunction testBasis \
testBlacsContext testSlaterDet testEnergyFunctional testSample \
testChargeDensity  testFourierTransform testSpecies testContext \
testEigenSolvers

QBSRC=AtomSet.C Atom.C SymmetrySet.C Symmetry.C Species.C Wavefunction.C SlaterDet.C \
      EnergyFunctional.C SampleStepper.C Basis.C FourierTransform.C Matrix.C Context.C \
      sinft.C spline.C UnitCell.C StructureFactor.C ChargeDensity.C UserInterface.C RunCmd.C \
      LoadCmd.C SaveCmd.C SavesysCmd.C SavedenCmd.C SpeciesCmd.C SpeciesReader.C \
      SpeciesHandler.C XCPotential.C LDAFunctional.C PBEFunctional.C PBESolFunctional.C \
      PBERevFunctional.C BLYPFunctional.C NonLocalPotential.C SampleReader.C SampleWriter.C \
      StructuredDocumentHandler.C SampleHandler.C AtomSetHandler.C WavefunctionHandler.C \
      XMLGFPreprocessor.C Base64Transcoder.C CPSampleStepper.C BOSampleStepper.C ParOptCmd.C \
      ParallelOptimizer.C WavefunctionStepper.C SDWavefunctionStepper.C \
      MDWavefunctionStepper.C SDIonicStepper.C MDIonicStepper.C BMDIonicStepper.C \
      PSDWavefunctionStepper.C PSDAWavefunctionStepper.C JDWavefunctionStepper.C \
      SDCellStepper.C ConfinementPotential.C Preconditioner.C SimpleConvergenceDetector.C \
      MMSpecies.C EmpiricalPotential.C EnthalpyFunctional.C Hugoniostat.C release.C qbox_xmlns.C \
      isodate.C ComputeMLWFCmd.C BasisMapping.C MLWFTransform.C jade.C ConstraintSet.C \
      Constraint.C DistanceConstraint.C PositionConstraint.C AngleConstraint.C TorsionConstraint.C \
      SaveESPCmd.C AndersonMixer.C SDAIonicStepper.C CGIonicStepper.C jacobi.C HubbardPotential.C \
      SphericalIntegration.C PrintMem.C
KPTSRC=SymOpSet.C SymOp.C

SIMDSRC = clooper.c

QBOBJS = $(QBSRC:%.C=$(OBJDIR)/%.o)
KPTOBJS = $(KPTSRC:%.C=$(OBJDIR)/%.o)
SIMDOBJS = $(SIMDSRC:%.c=$(OBJDIR)/%.o)
ifdef CFLAGS
	QBOBJS += $(SIMDOBJS)
endif

OBJDIR = objs-$(ARCH)
EXENAME = qb-$(ARCH)
CXXFLAGS += -DARCH='"$(ARCH)"'

INSTALLDIR = /usr/apps/qbox
EXEINSTALL = qb@ll-$(VERSION)
NRUNSINSTALL = qb@ll-nruns-$(VERSION)

DEFAULT: $(EXENAME)

$(OBJDIR)/%.o: %.C
	$(CXX) $(CXXFLAGS) -o $@ -c $^

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJDIR)/qb.o: qb.C
	$(CXX) $(CXXFLAGS) -o $@ -c $^

$(EXENAME): makedirs depend $(QBOBJS) BGQOMPhack $(OBJDIR)/qb.o
	$(LD) $(DFLAGS) -o $(EXENAME) $(OBJDIR)/qb.o $(QBOBJS) $(LDFLAGS)

BGQOMPhack:
ifdef OMPHACK
	$(CXX) -qsmp=omp $(CXXFLAGS) -o $(OBJDIR)/FourierTransform.o -c FourierTransform.C
	$(CXX) -qsmp=omp $(CXXFLAGS) -o $(OBJDIR)/ChargeDensity.o -c ChargeDensity.C
	$(CXX) -qsmp=omp $(CXXFLAGS) -o $(OBJDIR)/SlaterDet.o -c SlaterDet.C
	$(CXX) -qsmp=omp $(CXXFLAGS) -o $(OBJDIR)/NonLocalPotential.o -c NonLocalPotential.C
endif

install: $(EXENAME) qb-nruns
ifdef VERSION
	cp $(EXENAME) $(INSTALLDIR)/$(EXEINSTALL)
	chmod a+rx $(INSTALLDIR)/$(EXEINSTALL)
	cp qb-nruns-$(ARCH) $(INSTALLDIR)/$(NRUNSINSTALL)
	chmod a+rx $(INSTALLDIR)/$(NRUNSINSTALL)
	cp ChangeLog $(INSTALLDIR)
	chmod a+r $(INSTALLDIR)/ChangeLog
	cp qbLink.h $(INSTALLDIR)/qblink/qbLink-$(VERSION).h
	chmod a+r $(INSTALLDIR)/qblink/qbLink-$(VERSION).h
	cp $(OBJDIR)/libqblink-$(ARCH).a $(INSTALLDIR)/qblink/libqblink-$(VERSION).a
	chmod a+r $(INSTALLDIR)/qblink/libqblink-$(VERSION).a
else
	@echo "Need to define VERSION variable to make install, e.g. VERSION=1.59.2-zeus"
endif

link: $(OBJDIR)/qbLink.o $(QBOBJS)
ifeq ($(ARCH),ibm_bgl_xlc)
	mkdir $(OBJDIR)/qblink.tmp.$$$$ && cd $(OBJDIR)/qblink.tmp.$$$$ && ar -x $(SCALAPACKLIB) && ar -x $(FFTWLIB) && cd ../.. && ar -cr $(OBJDIR)/libqblink-$(ARCH).a $(OBJDIR)/qblink.tmp.$$$$/*.o $^ && rm -rf $(OBJDIR)/qblink.tmp.$$$$
else
ifeq ($(ARCH),ibm_bgp_xlc)
	mkdir $(OBJDIR)/qblink.tmp.$$$$ && cd $(OBJDIR)/qblink.tmp.$$$$ && ar -x $(SCALAPACKLIB) && cd ../.. && ar -cr $(OBJDIR)/libqblink-$(ARCH).a $(OBJDIR)/qblink.tmp.$$$$/*.o $^ && rm -rf $(OBJDIR)/qblink.tmp.$$$$
else
ifeq ($(ARCH),ibm_sp5_xlc)
	mkdir $(OBJDIR)/qblink.tmp.$$$$ && cd $(OBJDIR)/qblink.tmp.$$$$ && ar -x $(SCALAPACKLIB) && cd ../.. && ar -cr $(OBJDIR)/libqblink-$(ARCH).a $(OBJDIR)/qblink.tmp.$$$$/*.o $^ && rm -rf $(OBJDIR)/qblink.tmp.$$$$
else
	mkdir $(OBJDIR)/qblink.tmp.$$$$ && cd $(OBJDIR)/qblink.tmp.$$$$ && ar x $(SCALAPACKLIB) && ar x $(FFTWLIB) && ar x $(XERCESLIB) && cd ../.. && ar cr $(OBJDIR)/libqblink-$(ARCH).a $(OBJDIR)/qblink.tmp.$$$$/*.o $^ && rm -rf $(OBJDIR)/qblink.tmp.$$$$
endif
endif
endif

lib:   $(OBJDIR)/qb.o $(QBOBJS)
	ar cr $(OBJDIR)/libqb.a $^

qb-nruns: makedirs $(OBJDIR)/qb-nruns.o link
	$(LD) -o qb-nruns-$(ARCH) $(OBJDIR)/qb-nruns.o $(OBJDIR)/libqblink-$(ARCH).a $(LDFLAGS)

qb-setupkpts: makedirs $(OBJDIR)/qb-setupkpts.o $(QBOBJS) $(KPTOBJS)
	$(LD) -o qb-setupkpts-$(ARCH) $(OBJDIR)/qb-setupkpts.o $(QBOBJS) $(KPTOBJS) $(LDFLAGS)

makedirs:
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR) ; fi
	@if [ ! -e .depend_$(ARCH) ]; then touch .depend_$(ARCH) ; fi

testMatrix: $(OBJDIR)/testMatrix.o $(OBJDIR)/Matrix.o $(OBJDIR)/Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
testBasis: $(OBJDIR)/testBasis.o $(OBJDIR)/Basis.o $(OBJDIR)/UnitCell.o $(OBJDIR)/Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
testFourierTransform: $(OBJDIR)/testFourierTransform.o $(OBJDIR)/FourierTransform.o $(OBJDIR)/Basis.o \
	$(OBJDIR)/UnitCell.o $(OBJDIR)/Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
testBlas: $(OBJDIR)/testBlas.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
testEigenSolvers: $(OBJDIR)/testEigenSolvers.o $(OBJDIR)/Matrix.o $(OBJDIR)/Context.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
xmlSpecies: xmlSpecies.o qbox_xmlns.o
	$(LD) $(DFLAGS) -o $@ $^ $(LDFLAGS)
xmlget: xmlget.o
	$(LD) -o xmlget xmlget.C -I$(XERCESCDIR)/include \
        -L$(XERCESCDIR)/lib -lxerces-c 
xmlextract: xmlextract.o
	$(LD) -o xmlextract xmlextract.C -I$(XERCESCDIR)/include \
        -L$(XERCESCDIR)/lib -lxerces-c 
#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
 depend :
	@makedepend -f .depend_$(ARCH) -p'($(OBJDIR)/' -o '.o)' -Y -D$(PLT) *.[cCh] 2> /dev/null

.depend_$(ARCH): $(QBSRC)
	@touch .depend_$(ARCH)
	@$(MAKE) depend

#	makedepend -Y -D$(PLT) *.[cCh]
#makedepend -D$(PLT) -I/usr/include/CC *.[cC]
#------------------------------------------------------------------------------
#  Cleanup object files
 clean :
	rm -f $(OBJDIR)/*.o
	rm -f $(EXENAME) $(TESTEXECS)
	rm -rf ti_files
	rm -f $(OBJDIR)/libqb.a
	rm -f $(OBJDIR)/libqblink-$(ARCH).a
	rm -rf $(OBJDIR)/qblink.tmp.*
	rmdir $(OBJDIR)
	rm -f .depend_$(ARCH)

#------------------------------------------------------------------------------
 ctags :
	etags -o tags *.[Ch]
#------------------------------------------------------------------------------
 html :
	enscript -Ecpp --color -Whtml --toc -pqbsrc.html *.h *.C
#------------------------------------------------------------------------------
 ps :
	enscript -Ecpp --toc -2rGj -pqb.ps *.h *.C Makefile *.mk
#------------------------------------------------------------------------------
 tar :
	tar czf qbsrc.tgz *.[Chc] Makefile Makefile.arch ChangeLog *.mk
#------------------------------------------------------------------------------

-include .depend_$(ARCH)
