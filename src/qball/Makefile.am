#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.46 2010/08/26 17:44:16 draeger1 Exp $
#------------------------------------------------------------------------------

include $(srcdir)/../common-rules.make

bin_PROGRAMS  =                             \
	qball                               \
	qball-nruns                         \
	qball-setupkpts                     \
	qballdiff                           \
	xmlSpecies                          \
	xmlget                              \
	xmlextract

lib_LIBRARIES =                             \
	libqbLink.a

include_HEADERS =                           \
	qbLink.h

noinst_HEADERS =                            \
	AndersonMixer.h                     \
	AngleConstraint.h                   \
	Atom.h                              \
	AtomSet.h                           \
	AtomSetHandler.h                    \
	Base64Transcoder.h                  \
	Basis.h                             \
	BasisMapping.h                      \
	blacs.h                             \
	BMDIonicStepper.h                   \
	BOSampleStepper.h                   \
	CellStepper.h                       \
	CGIonicStepper.h                    \
	ChargeDensity.h                     \
	ClebschGordan.h                     \
	clooper.h                           \
	ConfinementPotential.h              \
	Constraint.h                        \
	ConstraintSet.h                     \
	Context.h                           \
	Control.h                           \
	ConvergenceDetector.h               \
	CoordinateConstraint.h              \
	CPSampleStepper.h                   \
	CurrentDensity.h                    \
	dftd3.h                             \
	D3vector.h                          \
	Dimensions.h                        \
	DistanceConstraint.h                \
	EhrenSampleStepper.h                \
	EmpiricalPotential.h                \
	EnergyFunctional.h                  \
	EnthalpyFunctional.h                \
	ExponentialWavefunctionStepper.h    \
	Extrapolator.h                      \
	ExtrapolatorASP.h                   \
	ExtrapolatorNTC.h                   \
	ExtrapolatorSimple.h                \
	FCPStepper.h                        \
	FcpThemostat.h                      \
	FORKTDWavefunctionStepper.h         \
	FOTDWavefunctionStepper.h           \
	FourierTransform.h                  \
	HubbardPotential.h                  \
	Hugoniostat.h                       \
	IonicStepper.h                      \
	isodate.h                           \
	jacobi.h                            \
	jade.h                              \
	JDWavefunctionStepper.h             \
	MDIonicStepper.h                    \
	MDWavefunctionStepper.h             \
	Messages.h                          \
	MLWFTransform.h                     \
	MMSpecies.h                         \
	NameOf.h                            \
	NonLocalPotential.h                 \
	NonLocalPotential.save.h            \
	ParallelOptimizer.h                 \
	PositionConstraint.h                \
	Preconditioner.h                    \
	PrintMem.h                          \
	profile.h                           \
	PSDAWavefunctionStepper.h           \
	PSDWavefunctionStepper.h            \
	qbox_xmlns.h                        \
	RMMDIISWavefunctionStepper.h        \
	readTSC.h                           \
	release.h                           \
	Sample.h                            \
	SampleHandler.h                     \
	SampleReader.h                      \
	SampleStepper.h                     \
	SampleWriter.h                      \
	SDAIonicStepper.h                   \
	SDCellStepper.h                     \
	SDIonicStepper.h                    \
	SDWavefunctionStepper.h             \
	SelfConsistentPotential.h           \
	SharedFilePtr.h                     \
	SimpleConvergenceDetector.h         \
	sinft.h                             \
	SlaterDet.h                         \
	SORKTDWavefunctionStepper.h         \
	SOTDWavefunctionStepper.h           \
	Species.h                           \
	SpeciesHandler.h                    \
	SpeciesReader.h                     \
	SphericalIntegration.h              \
	StructuredDocumentHandler.h         \
	StructureFactor.h                   \
	StructureHandler.h                  \
	StrX.h                              \
	Symmetry.h                          \
	SymmetrySet.h                       \
	SymOp.h                             \
	SymOpSet.h                          \
	TDEULERWavefunctionStepper.h        \
	TestHandler.h                       \
	Timer.h                             \
	TorsionConstraint.h                 \
	Unit.h                              \
	UnitCell.h                          \
	uuid_str.h                          \
	Wavefunction.h                      \
	WavefunctionHandler.h               \
	WavefunctionStepper.h               \
	XCPotential.h                       \
	XMLFile.h                           \
	XMLGFPreprocessor.h

libqbLink_a_SOURCES =                       \
	AtomSet.cc                           \
	Atom.cc                              \
	CoordinateConstraint.cc              \
	SymmetrySet.cc                       \
	Symmetry.cc                          \
	Species.cc                           \
	Wavefunction.cc                      \
	RMMDIISWavefunctionStepper.cc        \
	SlaterDet.cc                         \
	EnergyFunctional.cc                  \
	SampleStepper.cc                     \
	Basis.cc                             \
	CurrentDensity.cc                    \
	FourierTransform.cc                  \
	Context.cc                           \
	sinft.cc                             \
	UnitCell.cc                          \
	StructureFactor.cc                   \
	ChargeDensity.cc                     \
	SpeciesReader.cc                     \
	SpeciesHandler.cc                    \
	XCPotential.cc                       \
	NonLocalPotential.cc                 \
	SampleReader.cc                      \
	SampleWriter.cc                      \
	StructuredDocumentHandler.cc         \
	SampleHandler.cc                     \
	AtomSetHandler.cc                    \
	WavefunctionHandler.cc               \
	XMLGFPreprocessor.cc                 \
	Base64Transcoder.cc                  \
	CPSampleStepper.cc                   \
	BOSampleStepper.cc                   \
	ParallelOptimizer.cc                 \
	WavefunctionStepper.cc               \
	SDWavefunctionStepper.cc             \
	MDWavefunctionStepper.cc             \
	SDIonicStepper.cc                    \
	MDIonicStepper.cc                    \
	BMDIonicStepper.cc                   \
	PSDWavefunctionStepper.cc            \
	PSDAWavefunctionStepper.cc           \
	JDWavefunctionStepper.cc             \
	SDCellStepper.cc                     \
	ConfinementPotential.cc              \
	Preconditioner.cc                    \
	SimpleConvergenceDetector.cc         \
	MMSpecies.cc                         \
	EmpiricalPotential.cc                \
	EnthalpyFunctional.cc                \
	ExtrapolatorASP.cc                   \
	ExtrapolatorNTC.cc                   \
	ExtrapolatorSimple.cc                \
	Hugoniostat.cc                       \
	release.cc                           \
	qbox_xmlns.cc                        \
	isodate.cc                           \
	BasisMapping.cc                      \
	MLWFTransform.cc                     \
	jade.cc                              \
	ConstraintSet.cc                     \
	Constraint.cc                        \
	DistanceConstraint.cc                \
	PositionConstraint.cc                \
	AngleConstraint.cc                   \
	TorsionConstraint.cc                 \
	AndersonMixer.cc                     \
	SDAIonicStepper.cc                   \
	CGIonicStepper.cc                    \
	jacobi.cc                            \
	HubbardPotential.cc                  \
	SphericalIntegration.cc              \
	PrintMem.cc                          \
	SOTDWavefunctionStepper.cc           \
	SORKTDWavefunctionStepper.cc         \
	FORKTDWavefunctionStepper.cc         \
	TDEULERWavefunctionStepper.cc        \
	EhrenSampleStepper.cc                \
	FCPStepper.cc                        \
	ExponentialWavefunctionStepper.cc    \
	SelfConsistentPotential.cc           \
	clooper.c                           \
	qbLink.cc

qball_SOURCES =                             \
	qb.cc

qball_LDADD =                               \
	$(all_LIBS)

qball_nruns_SOURCES =                       \
	qb-nruns.cc

qball_nruns_LDADD =                         \
	$(all_LIBS)

qball_setupkpts_SOURCES =                   \
	qb-setupkpts.cc                      \
	SymOpSet.cc                          \
	SymOp.cc

qball_setupkpts_LDADD =                     \
	$(all_LIBS)

qballdiff_SOURCES =                         \
	qbdiff.cc

xmlSpecies_SOURCES =                        \
	xmlSpecies.cc                        \
	qbox_xmlns.cc

xmlget_SOURCES =                            \
	xmlget.cc

xmlextract_SOURCES =                        \
	xmlextract.cc
