external_LIBS =                                              \
	$(top_builddir)/external_libs/dftd3/libdftd3.a

internal_LIBS =                                              \
	$(top_builddir)/src/qball/libqbLink.a                \
	$(top_builddir)/src/functionals/libfunctionals.a     \
	$(top_builddir)/src/pseudo/libpseudo.a               \
	$(top_builddir)/src/math/libmath.a

all_LIBS = $(internal_LIBS) $(external_LIBS)

AM_CXXFLAGS = -I$(top_srcdir)/src/
