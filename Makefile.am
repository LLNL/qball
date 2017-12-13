ACLOCAL_AMFLAGS = -I m4

SUBDIRS =        \
	external_libs \
	src      \
	testsuite

check-short:
	$(MAKE) -C src check
	$(MAKE) -C testsuite check-short

check-long:
	$(MAKE) -C src check
	$(MAKE) -C testsuite check-long

