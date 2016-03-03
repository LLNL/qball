
%.run: %.test
	../oct-run_regression_test.pl -D ../src/ -f $< -p
