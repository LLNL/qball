#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.46 2010/08/26 17:44:16 draeger1 Exp $
#------------------------------------------------------------------------------
#
# Get name of architecture-specific include file from hostname

bin_PROGRAMS  =                             \
	qball                               \
	qball-nruns                         \
	qball-setupkpts                     \
	qballdiff                           \
	xmlSpecies                          \
	xmlget                              \
	xmlextract                          \
	testMatrix                          \
	testBasis                           \
	testRowBcast                        \
	testFourierTransform                \
	testBlas                            \
	testGram                            \
	testGramBlock                       \
	testPgemm                           \
	testPgemmBlock                      \
	testPzheev                          \
	testEigenSolvers                    \
	testEigenBlock

lib_LIBRARIES =                             \
	libqbLink.a

include_HEADERS =                           \
	qbLink.h

noinst_HEADERS =                            \
	AndersonMixer.h                     \
	AngleCmd.h                          \
	AngleConstraint.h                   \
	AtomCmd.h                           \
	Atom.h                              \
	AtomsDyn.h                          \
	AtomSet.h                           \
	AtomSetHandler.h                    \
	Base64Transcoder.h                  \
	Basis.h                             \
	BasisMapping.h                      \
	blacs.h                             \
	blas.h                              \
	BLYPFunctional.h                    \
	BMDIonicStepper.h                   \
	BOSampleStepper.h                   \
	CellDyn.h                           \
	Cell.h                              \
	CellLock.h                          \
	CellMass.h                          \
	CellStepFreq.h                      \
	CellStepper.h                       \
	CenterOfMass.h                      \
	CGIonicStepper.h                    \
	ChargeDensity.h                     \
	ChargeMixCoeff.h                    \
	ChargeMixing.h                      \
	ChargeMixNdim.h                     \
	ChargeMixRcut.h                     \
	ClebschGordan.h                     \
	clooper.h                           \
	ComputeMLWFCmd.h                    \
	ConfinementPotential.h              \
	ConstraintCmd.h                     \
	Constraint.h                        \
	ConstraintSet.h                     \
	Context.h                           \
	Control.h                           \
	ConvergenceDetector.h               \
	CPSampleStepper.h                   \
	CurrentDensity.h                    \
	D3vector.h                          \
	Debug.h                             \
	Dimensions.h                        \
	DistanceCmd.h                       \
	DistanceConstraint.h                \
	Dt.h                                \
	Ecutden.h                           \
	Ecut.h                              \
	Ecutprec.h                          \
	Ecuts.h                             \
	EhrenSampleStepper.h                \
	Emass.h                             \
	EmpiricalPotentialCmd.h             \
	EmpiricalPotential.h                \
	EnergyFunctional.h                  \
	EnthalpyFunctional.h                \
	EnthalpyPressure.h                  \
	EnthalpyThreshold.h                 \
	EsmBC.h                             \
	EsmW.h                              \
	ExponentialWavefunctionStepper.h    \
	Extrapolator.h                      \
	ExtrapolatorASP.h                   \
	ExtrapolatorNTC.h                   \
	ExtrapolatorSimple.h                \
	ExtStress.h                         \
	FcpMu.h                             \
	FcpPmass.h                          \
	FCPStepper.h                        \
	FcpThemostat.h                      \
	FcpThermostat.h                     \
	FcpThTemp.h                         \
	FcpThTime.h                         \
	FcpThWidth.h                        \
	FermiTemp.h                         \
	FoldInWsCmd.h                       \
	Force_Complex_WF.h                  \
	FORKTDWavefunctionStepper.h         \
	FOTDWavefunctionStepper.h           \
	FourierTransform.h                  \
	HelpCmd.h                           \
	HubbardPotential.h                  \
	HubbardU.h                          \
	HugDeltaTemp.h                      \
	HugFreq.h                           \
	Hugoniostat.h                       \
	HugoniostatVar.h                    \
	IonicStepper.h                      \
	IPrint.h                            \
	isodate.h                           \
	jacobi.h                            \
	jade.h                              \
	JDWavefunctionStepper.h             \
	KpointCmd.h                         \
	LDAFunctional.h                     \
	LineMinimizer.h                     \
	ListAtomsCmd.h                      \
	ListConstraintsCmd.h                \
	ListSpeciesCmd.h                    \
	LoadCmd.h                           \
	LoaddenCmd.h                        \
	LockCmd.h                           \
	Matrix.h                            \
	MatrixLoc.h                         \
	MDIonicStepper.h                    \
	MDIter.h                            \
	MDSaveCmd.h                         \
	MDWavefunctionStepper.h             \
	Memory.h                            \
	MLWFTransform.h                     \
	MMAtomCmd.h                         \
	MMSpeciesCmd.h                      \
	MMSpecies.h                         \
	NameOf.h                            \
	NA_overlaps.h                       \
	Nempty.h                            \
	NetCharge.h                         \
	Nkpoints.h                          \
	NonLocalPotential.h                 \
	NonLocalPotential.save.h            \
	Non_Selfconsistent_Energy_Output.h  \
	Nparallelkpts.h                     \
	Nrowmax.h                           \
	ParallelOptimizer.h                 \
	ParOptCmd.h                         \
	PBEFunctional.h                     \
	PBERevFunctional.h                  \
	PBESolFunctional.h                  \
	Pblock.h                            \
	PlotCmd.h                           \
	PositionConstraint.h                \
	Preconditioner.h                    \
	PrintCmd.h                          \
	PrintMem.h                          \
	profile.h                           \
	PromoteOccCmd.h                     \
	PSDAWavefunctionStepper.h           \
	PSDWavefunctionStepper.h            \
	qbox_xmlns.h                        \
	QuitCmd.h                           \
	RandomizeRealWfCmd.h                \
	RandomizeVelCmd.h                   \
	RandomizeWfCmd.h                    \
	readTSC.h                           \
	RefCell.h                           \
	release.h                           \
	ResetVcmCmd.h                       \
	RunCmd.h                            \
	RunTimer.h                          \
	Sample.h                            \
	SampleHandler.h                     \
	SampleReader.h                      \
	SampleStepper.h                     \
	SampleWriter.h                      \
	SaveCmd.h                           \
	SavedenCmd.h                        \
	SaveDenFreq.h                       \
	SaveESPCmd.h                        \
	SaveFreq.h                          \
	SavesysCmd.h                        \
	SaveWfFreq.h                        \
	SDAIonicStepper.h                   \
	SDCellStepper.h                     \
	SDIonicStepper.h                    \
	SDWavefunctionStepper.h             \
	SelfConsistentPotential.h           \
	SetCmd.h                            \
	SetVelCmd.h                         \
	SharedFilePtr.h                     \
	ShiftWFCmd.h                        \
	SimpleConvergenceDetector.h         \
	sinft.h                             \
	SlaterDet.h                         \
	Smearing.h                          \
	SmearingWidth.h                     \
	SORKTDWavefunctionStepper.h         \
	SOTDWavefunctionStepper.h           \
	SpeciesCmd.h                        \
	Species.h                           \
	SpeciesHandler.h                    \
	SpeciesReader.h                     \
	SphericalIntegration.h              \
	Spin.h                              \
	Spline.h                            \
	StandardVar.h                       \
	StatusCmd.h                         \
	StrainCmd.h                         \
	Stress.h                            \
	StructuredDocumentHandler.h         \
	StructureFactor.h                   \
	StructureHandler.h                  \
	StrX.h                              \
	SymmetryCmd.h                       \
	Symmetry.h                          \
	SymmetrySet.h                       \
	SymOp.h                             \
	SymOpSet.h                          \
	TDDt.h                              \
	TDEULERWavefunctionStepper.h        \
	TestHandler.h                       \
	Thermostat.h                        \
	ThresholdForce.h                    \
	ThresholdScf.h                      \
	ThresholdStress.h                   \
	ThTemp.h                            \
	ThTime.h                            \
	ThWidth.h                           \
	Timer.h                             \
	TorsionCmd.h                        \
	TorsionConstraint.h                 \
	Unit.h                              \
	UnitCell.h                          \
	UnlockCmd.h                         \
	UserInterface.h                     \
	uuid_str.h                          \
	Wavefunction.h                      \
	WavefunctionHandler.h               \
	WavefunctionStepper.h               \
	WfDiag.h                            \
	WfDyn.h                             \
	WfExtrap.h                          \
	WFPhaseRealCmd.h                    \
	WFPhaseRealVar.h                    \
	XCFunctional.h                      \
	Xc.h                                \
	XCPotential.h                       \
	XMLFile.h                           \
	XMLGFPreprocessor.h

libqbLink_a_SOURCES =                        \
	AtomSet.C                           \
	Atom.C                              \
	SymmetrySet.C                       \
	Symmetry.C                          \
	Species.C                           \
	Wavefunction.C                      \
	SlaterDet.C                         \
	EnergyFunctional.C                  \
	SampleStepper.C                     \
	Basis.C                             \
	CurrentDensity.C                    \
	FourierTransform.C                  \
	Matrix.C                            \
	Context.C                           \
	sinft.C                             \
	Spline.C                            \
	UnitCell.C                          \
	StructureFactor.C                   \
	ChargeDensity.C                     \
	UserInterface.C                     \
	RunCmd.C                            \
	LoadCmd.C                           \
	SaveCmd.C                           \
	SavesysCmd.C                        \
	SavedenCmd.C                        \
	SpeciesCmd.C                        \
	SpeciesReader.C                     \
	SpeciesHandler.C                    \
	XCPotential.C                       \
	LDAFunctional.C                     \
	PBEFunctional.C                     \
	PBESolFunctional.C                  \
	PBERevFunctional.C                  \
	BLYPFunctional.C                    \
	NonLocalPotential.C                 \
	SampleReader.C                      \
	SampleWriter.C                      \
	StructuredDocumentHandler.C         \
	SampleHandler.C                     \
	AtomSetHandler.C                    \
	WavefunctionHandler.C               \
	XMLGFPreprocessor.C                 \
	Base64Transcoder.C                  \
	CPSampleStepper.C                   \
	BOSampleStepper.C                   \
	ParOptCmd.C                         \
	ParallelOptimizer.C                 \
	WavefunctionStepper.C               \
	SDWavefunctionStepper.C             \
	MDWavefunctionStepper.C             \
	SDIonicStepper.C                    \
	MDIonicStepper.C                    \
	BMDIonicStepper.C                   \
	PSDWavefunctionStepper.C            \
	PSDAWavefunctionStepper.C           \
	JDWavefunctionStepper.C             \
	SDCellStepper.C                     \
	ConfinementPotential.C              \
	Preconditioner.C                    \
	SimpleConvergenceDetector.C         \
	MMSpecies.C                         \
	EmpiricalPotential.C                \
	EnthalpyFunctional.C                \
	ExtrapolatorASP.C                   \
	ExtrapolatorNTC.C                   \
	ExtrapolatorSimple.C                \
	Hugoniostat.C                       \
	release.C                           \
	qbox_xmlns.C                        \
	isodate.C                           \
	ComputeMLWFCmd.C                    \
	BasisMapping.C                      \
	MLWFTransform.C                     \
	jade.C                              \
	ConstraintSet.C                     \
	Constraint.C                        \
	DistanceConstraint.C                \
	PositionConstraint.C                \
	AngleConstraint.C                   \
	TorsionConstraint.C                 \
	SaveESPCmd.C                        \
	AndersonMixer.C                     \
	SDAIonicStepper.C                   \
	CGIonicStepper.C                    \
	jacobi.C                            \
	HubbardPotential.C                  \
	SphericalIntegration.C              \
	PrintMem.C                          \
	PlotCmd.C                           \
	SOTDWavefunctionStepper.C           \
	SORKTDWavefunctionStepper.C         \
	FORKTDWavefunctionStepper.C         \
	TDEULERWavefunctionStepper.C        \
	EhrenSampleStepper.C                \
	MDSaveCmd.C                         \
	FCPStepper.C                        \
	ExponentialWavefunctionStepper.C    \
	SelfConsistentPotential.C           \
	clooper.c                           \
	qbLink.C

qball_SOURCES =                             \
	qb.C

qball_LDADD =                               \
	libqbLink.a

qball_nruns_SOURCES =                       \
	qb-nruns.C

qball_nruns_LDADD =                         \
	libqbLink.a

qball_setupkpts_SOURCES =                   \
	qb-setupkpts.C                      \
	SymOpSet.C                          \
	SymOp.C

qball_setupkpts_LDADD =                     \
	libqbLink.a

qballdiff_SOURCES =                         \
	qbdiff.C

testMatrix_SOURCES =                        \
	testMatrix.C                        \
	Matrix.C                            \
	Context.C

testBasis_SOURCES =                         \
	testBasis.C                         \
	Basis.C                             \
	UnitCell.C                          \
	Context.C

testRowBcast_SOURCES =                       \
	testRowBcast.C                      \
	Context.C

testFourierTransform_SOURCES =              \
	testFourierTransform.C              \
	FourierTransform.C                  \
	Basis.C                             \
	UnitCell.C                          \
	Context.C

testBlas_SOURCES =                          \
	testBlas.C

testGram_SOURCES =                          \
	testGram.C                          \
	Matrix.C                            \
	Context.C

testGramBlock_SOURCES =                     \
	testGramBlock.C                     \
	Matrix.C                            \
	Context.C

testPgemm_SOURCES =                         \
	testPgemm.C                         \
	Matrix.C                            \
	Context.C

testPgemmBlock_SOURCES =                    \
	testPgemmBlock.C                    \
	Matrix.C                            \
	Context.C

testPzheev_SOURCES =                        \
	testPzheev.C                        \
	Matrix.C                            \
	Context.C

testEigenSolvers_SOURCES =                  \
	testEigenSolvers.C                  \
	Matrix.C                            \
	Context.C                           \
	jacobi.C

testEigenBlock_SOURCES =                    \
	testEigenBlock.C                    \
	Matrix.C                            \
	Context.C                           \
	jacobi.C

xmlSpecies_SOURCES =                        \
	xmlSpecies.C                        \
	qbox_xmlns.C

xmlget_SOURCES =                            \
	xmlget.C

xmlextract_SOURCES =                        \
	xmlextract.C
