#!/usr/bin/env bash
#
# Copyright (C) 2005-2015 Heiko Appel, David Strubbe, Xavier Andrade
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# $Id: oct-run_testsuite.sh.in 14625 2015-10-03 18:33:06Z xavier $

# Paths.
prefix=../
# we do not care about this one, but without it we get a warning from configure:
# config.status: WARNING:  'testsuite/oct-run_testsuite.in' seems to ignore the --datarootdir setting
# this variable is not substituted by older versions of autoconf, e.g. 2.59
# NOTE: datarootdir line must precede pkgdatadir since @datadir@ may be substituted as ${datarootdir}
testsuite=.

# Failure reporting.
failed_tests=0
skipped_tests=0
passed_tests=0
total_tests=0
failure_report=""
NUM=0

# Usage.
function usage() {
    cat <<EOF

 Copyright (C) 2005-2006 by Heiko Appel

Usage: oct-run_testsuite.sh [options]
    
     -h            this message
     -n            dry-run mode (show what would be executed)
     -g LIST       comma-separated list of test groups to run
     -q            query testfiles for the given groups (no tests are run)
     -d DIR        directory where to look for the testsuite
     -l            local run
     -p PREFIX     installation prefix [default: /usr]
     -c            delete all .log files and work directories after the run

Report bugs to <octopus-devel@tddft.org>.
EOF
 exit 0;
}


# Find tests to run. Takes as argument the string passed to the
# -g command line option.
function find_tests() {
    groups="$1"
    if [ -n "$groups" ]; then
	groups_pat=`echo "$groups" | sed 's/ *, */|/g'`
	groups_pat="^TestGroups *:.*($groups_pat)"
    else # No groups given defaults to all tests.
	groups_pat="."
    fi
    testfiles=`find $testsuite -name "*.test" | \
	xargs egrep "$groups_pat" | \
	awk -F: '{print $1}' | sort -u`
    echo "$testfiles"
}

function queue {
	OCT_QUEUE="$OCT_QUEUE $1 "
	OUTLIST="$OUTLIST $2 "
	NUM=$(($NUM+1))
}

function inqueue {
  if [ $OCT_TEST_NJOBS -gt 1 ]
  then
    printf "In execution: $OUTLIST\n\n"
  else
    printf "\n"
  fi
}


function checkqueue {
    inqueue

    still_waiting=1;
    while [ $still_waiting = 1 ]; do
      # find first test in queue which has finished
      for (( ijob=1; ijob <= $NUM; ijob++ )); do
	testpid=`echo $OCT_QUEUE | cut -f $ijob -d " "`
        test_name=`echo $OUTLIST | cut -f $ijob -d " "`
	#echo "checking $test_name $testpid"
        ps -p $testpid > /dev/null;
	# is this pid still running?
        if [ $? -ne 0 ]; then
	  still_waiting=0;
	  # remove it from the queue

          if [ $NUM -eq 1 ]; then
            OCT_QUEUE=
	    OUTLIST=
          elif [ $ijob -eq 1 ]; then
            OCT_QUEUE=`echo $OCT_QUEUE | cut -f 2- -d " "`
	    OUTLIST=`  echo $OUTLIST   | cut -f 2- -d " "`
          else
	    # cut will not accept 0 as an index
            OCT_QUEUE=`echo $OCT_QUEUE | cut -f -$((ijob-1)),$((ijob+1))- -d " "`
	    OUTLIST=`  echo $OUTLIST   | cut -f -$((ijob-1)),$((ijob+1))- -d " "`
          fi
          NUM=$(($NUM-1))
          break;
        fi
      done
      sleep 0.1s;
    done

    logfile=`basename $test_name .test`.log

    failures=`grep "Status:" $logfile| cut -f 2 -d \ `

    if [ x${failures} == x ]; then
	failures="error";
    fi
    
    if [ x${failures} == xskipped ] || [ x${failures} == x0 ]; then
	# Remove Status: line and following blank one.
	sed 'N;$!P;$!D;$d' $logfile
    else
	cat $logfile
    fi

    # Look for failed/passed/skipped testcases and add to failure summary.
    if [ ${failures} == "0" ]; then
	passed_tests=$((passed_tests + 1))
    elif [ ${failures} == "skipped" ]; then
	skipped_tests=$((skipped_tests + 1))
    else
	failed_tests=$((failed_tests + 1))
	name_length=`echo -n $test_name | wc -m`
	space_length=$((60 - name_length))
	spaces=`for((i = 1; i <= space_length; i++)); do echo -n ' '; done`
	failure_report="$failure_report    $test_name$spaces$failures\n"
    fi

    printf "\n ************************\n\n"
}

# Run all tests. Takes as first argument a list of testfile names.
function run_testsuite() {

tests="$1"
if [ -n "$test_groups" ]; then
echo "Test groups: $test_groups"
fi
echo ""

# Check for 'preserve working directories' flag.
if [ "$cleanup" == "yes" ]; then
    preserve_opt=""
else
    preserve_opt="-p"
fi

test_index=0

for y in $tests; do
    total_tests=$((total_tests + 1))
    test_name=`echo $y | sed "s|$testsuite/||"`
    ybase=`basename $y .test`
    if [ "${local_run}" == "yes" ]; then
	runtest_directory=$testsuite
    else
	runtest_directory=$bin_directory
    fi

    test_index=$(($test_index + 1))
    echo Starting test $test_index: $test_name

    if [ "$dry_run" == "yes" ]; then
	echo "	$runtest_directory/oct-run_regression_test.pl -l $opt_n \
	     $preserve_opt -D $bin_directory -f $y > ${ybase}.log 2>&1 &"
    fi
    $runtest_directory/oct-run_regression_test.pl -l $opt_n \
        $preserve_opt -D $bin_directory -f $y > ${ybase}.log 2>&1 &
    ppid=$!
    queue $ppid $test_name
	
    if [ $NUM -ge $OCT_TEST_NJOBS ]
    then
	checkqueue
    fi
done
while [ $NUM -ge 1 ]
do
    checkqueue
done
}


function trap_kill() {
 kill $OCT_QUEUE > /dev/null 2>&1
 echo "Killed"
 exit 1
}

# Show usage info if no args at all.
[ "$#" -eq 0 ] && usage;

# Parse command line.

# Some default settings.
query="no"
cleanup="no"
test_groups=""
dry_run="no"

while getopts "hnlp:d:cg:q:" opt ; do
    case "$opt" in
        h) usage;;
        n) dry_run="yes";; 
        l) local_run="yes";;
        p) prefix="$OPTARG"; testsuite=$prefix/share/octopus/testsuite;;
	d) directory="$OPTARG";;
	c) cleanup="yes";;
	g) test_groups="$OPTARG";;
	q) test_groups="$OPTARG"; query="yes";;
        ?) echo "Error parsing arguments" 1>&2; exit 1;;
    esac
done
shift $[ OPTIND - 1 ]

echo "*****************************"
echo "  Running qball testsuite  "
echo "*****************************"
echo
echo " Note: Running the testsuite can take from several minutes"
echo " to more than two hours."
echo

if [ "${local_run}" == "yes" ]; then
    bin_directory=`pwd`/../src/
else
    bin_directory=$prefix/bin
fi
if [ ! -d "${bin_directory}" ]; then
    echo "Specified binary directory '$bin_directory' does not exist." 1>&2
    exit 1
fi

if [ -z $OCT_TEST_NJOBS ] ; then

    echo
    echo " The environment variable OCT_TEST_NJOBS is not set. qball"
    echo " will determine automatically an appropriate number of jobs:"
    echo
    
    # some Linux
    nproc=`which nproc 2> /dev/null`;
    if [ -n $nproc ] && [ -x $nproc ]; then
	#We need to unset OMP_NUM_THREADS when counting the processors
	OCT_TEST_NJOBS=`OMP_NUM_THREADS= $nproc`
    fi

    # Mac OSX
    if [ -z $OCT_TEST_NJOBS ]; then
	nproc=`which sysctl 2> /dev/null`
	if [ -n $nproc ] && [ -x $nproc ]; then
	    OCT_TEST_NJOBS=`$nproc -n hw.ncpu`
	fi
    fi

    if [ -z $OCT_TEST_NJOBS ] ; then
	OCT_TEST_NJOBS=1
    else
	echo "  - Found $OCT_TEST_NJOBS processors that can be used."
    fi

    if [ $OCT_TEST_NJOBS -gt 1 ] && [ -z $EXEC ]; then
	# execute the commands using nice, to lower the priority and avoid core anxiety
	nice=`which nice 2> /dev/null`
	if [ -n $nice ] && [ -x $nice ] ; then
	    export EXEC=$nice
	fi
    fi
    
    options=" "`$bin_directory/qball -c`" "
    if [[ $options =~ " mpi " ]]; then
	OCT_TEST_NJOBS=$(($OCT_TEST_NJOBS/4))
	echo "  - This is an MPI run. Each job can use up to 4 processors."
    fi

    if [[ $options =~ " openmp " ]]; then
	if [ ! -z $OMP_NUM_THREADS ] ; then
	    OCT_TEST_NJOBS=$(($OCT_TEST_NJOBS/$OMP_NUM_THREADS))
	    echo "  - This is an OpenMP run. OMP_NUM_THREADS is set to $OMP_NUM_THREADS processors per job."
        else
	    echo "  - This is an OpenMP run. Each jobs uses all processors."
	    OCT_TEST_NJOBS=1
	fi
    fi


    if [ $OCT_TEST_NJOBS -eq 0 ] ; then
	OCT_TEST_NJOBS=1
    fi

    echo
    echo " Running $OCT_TEST_NJOBS jobs in parallel. Set environment variable OCT_TEST_NJOBS to adjust this."
    echo
    
    #wait for a bit for the user to read
    sleep 0.75
else
    echo
    echo " Running $OCT_TEST_NJOBS jobs in parallel."
    echo    
fi

# Find testfiles.
if [ -n "$directory" ]; then
    testsuite="$directory"
else
    [ "$local_run" == "yes" ] && testsuite=$(pwd)
fi

# remove any '/' characters so directories (with tab-completion) work as group names
test_groups=`echo $test_groups | tr -d '/'`

testfiles=`find_tests "$test_groups"`

# Query mode? If so, list files and exit.
if [ "$query" == "yes" ]; then
    echo "Testfiles for groups $test_groups:"
    echo ""
    for f in $testfiles; do
	echo ${f##$testsuite/}
    done
    exit 0
fi


# No testfiles found, abort.
if [ -z "$testfiles" ]; then
    echo "No testfiles for group(s) $test_groups found."

# Otherwise, start the whole machinery.
else
    # Get epoch seconds at testsuite start.
    testsuite_start=$(date +%s) 

    # trap signals, so we kill all children before dying
    trap "trap_kill" TERM INT

    if [ "$dry_run" == "yes" ]; then
	opt_n="-n"
    else
	opt_n=""
    fi

    # Run testsuite.
    run_testsuite "$testfiles"

    # Failure reporting to STDOUT.
    echo -e "    Passed:  $passed_tests / $total_tests"
    echo -e "    Skipped: $skipped_tests / $total_tests"
    if [ $failed_tests -gt 0 ]; then
	echo "    Failed:  $failed_tests / $total_tests"
	echo
	echo "    testfile                                                    # failed testcases"
	echo "    ------------------------------------------------------------------------------"
	echo -e "$failure_report"
    else
	if [ $passed_tests -gt 0 ]; then
	    echo -e "\nEverything seems to be OK"
	else
	    echo -e "\nAll tests were skipped."
	    # make sure a failure will be reported by the exit status
	    failed_tests=100
	fi
    fi
    echo
    
    # Clean up.
    [ "$cleanup" == "yes" ] && rm -f *.log

    # Get epoch seconds after we are done and compute time difference.
    testsuite_end=$(date +%s)
    timediff_sec=$[ testsuite_end - testsuite_start ]
    
    RUNTIME="Total run-time of the testsuite: \
    $(printf '%02d:%02d:%02d' $[timediff_sec / 3600] \
    $[(timediff_sec % 3600) / 60 ] $[timediff_sec % 60])"
    
    echo $RUNTIME
    echo ""
fi

exit $failed_tests


# Local Variables:
# mode: shell-script
# coding: utf-8
# End:
