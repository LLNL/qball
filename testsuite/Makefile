.PHONY: all check check-short check-long

all:
check:
	local_run=yes ./oct-run_testsuite.sh -g short-run,long-run

check-short:
	local_run=yes ./oct-run_testsuite.sh -g short-run

check-long:
	local_run=yes ./oct-run_testsuite.sh -g long-run
