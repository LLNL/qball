#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.46 2010/08/26 17:44:16 draeger1 Exp $
#------------------------------------------------------------------------------

include $(srcdir)/../common-rules.make

noinst_PROGRAMS  =                             \
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

LDADD = $(all_LIBS)

testMatrix_SOURCES = testMatrix.cc
testBasis_SOURCES = testBasis.cc
testRowBcast_SOURCES = testRowBcast.cc 
testFourierTransform_SOURCES = testFourierTransform.cc             
testBlas_SOURCES = testBlas.cc
testGram_SOURCES = testGram.cc
testGramBlock_SOURCES = testGramBlock.cc 
testPgemm_SOURCES = testPgemm.cc 
testPgemmBlock_SOURCES = testPgemmBlock.cc
testPzheev_SOURCES = testPzheev.cc
testEigenSolvers_SOURCES = testEigenSolvers.cc
testEigenBlock_SOURCES = testEigenBlock.cc 
