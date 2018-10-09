#-------------------------------------------------------------------------------
# $Id: Makefile,v 1.46 2010/08/26 17:44:16 draeger1 Exp $
#------------------------------------------------------------------------------

all: $(bin_PROGRAMS)

include $(srcdir)/../common-rules.make

bin_PROGRAMS  =                             \
	qball                               \
	qball-nruns                         \
	qball-setupkpts                     \
	qballdiff

noinst_HEADERS =                            \
	SymOp.h                             \
	SymOpSet.h

qball_SOURCES =                              \
	qb.cc

qball_LDADD =                                \
	$(all_LIBS)

qball_nruns_SOURCES =                        \
	qb-nruns.cc

qball_nruns_LDADD =                          \
	$(all_LIBS)

qball_setupkpts_SOURCES =                    \
	qb-setupkpts.cc                      \
	SymOpSet.cc                          \
	SymOp.cc

qball_setupkpts_LDADD =                      \
	$(all_LIBS)

qballdiff_SOURCES =                          \
	qbdiff.cc

