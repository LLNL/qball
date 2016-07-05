.PHONY: all check check-short check-long

EXTRA_DIST =                                 \
	oct-run_regression_test.pl           \
	oct-run_testsuite.sh

# OCTOPUS_SHARE has to be set to an absolute path because
# oct-run_testsuite.sh changes the working directory for every
# testcase.
# If the testsuite should be skipped, e.g. in make distcheck,
# set the environment variable SKIP_CHECK to some value.
check:
	@CHECKOPTS="" $(MAKE) check-wrapper

check-short:
	@CHECKOPTS="-g short-run" $(MAKE) check-wrapper

check-long:
	@CHECKOPTS="-g long-run" $(MAKE) check-wrapper

check-dryrun:
	@CHECKOPTS="-n" $(MAKE) check-base

check-wrapper:
	@if test "x$(SKIP_CHECK)" == "x"; then \
	  $(MAKE) check-base; \
	else \
	$(MAKE) check-dryrun; \
	fi

check-base:
	@base=`basename "$(top_builddir)/share"` && \
	dir=`dirname "$(top_builddir)/share"` && \
	sharedir="`(cd \"$$dir\" 2> /dev/null && pwd || echo \"$$dir\")`/$$base" && \
	QBALL_SHARE=$$sharedir $(top_srcdir)/testsuite/oct-run_testsuite.sh -c -d $(srcdir) -l $(CHECKOPTS)


