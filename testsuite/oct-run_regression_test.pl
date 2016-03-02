#!/usr/bin/env perl
#
# Copyright (C) 2005-2014 H. Appel, M. Marques, X. Andrade, D. Strubbe
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
# $Id: oct-run_regression_test.pl 15022 2016-01-09 04:19:06Z xavier $

use warnings;
use Getopt::Std;
use File::Basename;
use Fcntl ':mode';
use Time::HiRes qw(gettimeofday tv_interval);
use Scalar::Util qw(looks_like_number);

sub usage {

  print <<EndOfUsage;

 Copyright (C) 2005-2016 H. Appel, M. Marques, X. Andrade, D. Strubbe

Usage: oct-run_regression_test.pl [options]

    -n        dry-run
    -v        verbose
    -h        this usage
    -D        name of the directory where to look for the executables   
    -s        run everything serial
    -f        filename of testsuite [required]
    -p        preserve working directories
    -l        copy output log to current directory
    -m        run matches only (assumes there are work directories)

Exit codes:
    0         all tests passed
    1..253    number of test failures
    254       test skipped
    255       internal error

Report bugs to <octopus-devel\@tddft.org>
EndOfUsage

  exit 0;
}


sub set_precision{
  my $p = $_[0];
  if($p ne "default"){
    $precnum = 1.0*$p;
  } else {
    $precnum = 0.0001
  }
}

# Check whether STDOUT is a terminal. If not, no ANSI sequences are
# emitted.
if(-t STDOUT) {
    $color_start{blue}="\033[34m";
    $color_end{blue}="\033[0m";
    $color_start{red}="\033[31m";
    $color_end{red}="\033[0m";
    $color_start{green}="\033[32m";
    $color_end{green}="\033[0m";
} else {
    $color_start{blue}="";
    $color_end{blue}="";
    $color_start{red}="";
    $color_end{red}="";
    $color_start{green}="";
    $color_end{green}="";
}

if (not @ARGV) { usage; }

$opt_f = "";
getopts("nlvhD:c:f:spm");

# avoid warnings 'used only once: possible typo'
$useless = $opt_h;
$useless = $opt_l;
$useless = $opt_s;

# Default values
use File::Temp qw/tempdir/;

# Handle options
$opt_h && usage;

my $exec_directory;
if($opt_D) {
 $exec_directory = $opt_D;
 if($exec_directory !~ /^\//){
  $exec_directory = get_env("PWD")."/$exec_directory";
 }
} else {
 $exec_directory = "/usr/bin";
}

if(length($opt_f) == 0) {
    die255("You must supply the name of a test file with the -f option.");
}

$aexec = get_env("EXEC");
$global_np = get_env("OCT_TEST_MPI_NPROCS");

# FIXME: all test files should declare Processors
#$np = "serial";
$is_parallel = 0;

# FIXME: could bake in mpiexec at configure time

if(!$opt_s) {
# MPI stuff
    $mpiexec = get_env("MPIEXEC");
    $machinelist = get_env("MACHINELIST");
    if ("$mpiexec" eq "") { $mpiexec = `which mpiexec 2> /dev/null`; }
    chomp($mpiexec);

    if( "$mpiexec" eq "" ) { 
	print "No mpiexec found: running in serial.\n\n";
    } else {
# mpiexec without arguments (to check if it is available)
	$mpiexec_raw = $mpiexec;
	$mpiexec_raw =~ s/\ (.*)//;
	if ( ! -e "$mpiexec_raw" ) {
	    print "mpiexec command ($mpiexec_raw) does not exist: running in serial.\n\n";
	    $mpiexec = "";
	} elsif( ! -x "$mpiexec_raw" ) {
	    print "mpiexec command ($mpiexec_raw) is not executable: running in serial.\n\n";
	    $mpiexec = "";
	} else {
# default number of processors is 1
	    $np = 1;
	    $is_parallel = 1;
	}
    }
} else {
    $mpiexec = "";
}

# default number of processors for MPI runs is 2
$np = 2;
$enabled = ""; # FIXME: should Enabled be optional?
$options_required = "";
$options_required_mpi = "";
$options_are_mpi = 0;

# This variable counts the number of failed testcases.
$failures = 0;

$tempdirpath = get_env("TEMPDIRPATH");
if ("$tempdirpath" eq "") { $tempdirpath = '/tmp'; }
if (! -d $tempdirpath) { mkdir $tempdirpath; }

set_precision("default");
$test_succeeded = 1;

$pwd = get_env("PWD");
if (!$opt_m) {
    my $name = $opt_f;
    $name =~ s/\.\.\///g;
    $name =~ s/\//-/g;
    $workdir = tempdir("$tempdirpath/octopus" . "-" . $name . ".XXXXXX");
    chomp($workdir);

    system ("rm -rf $workdir");
    mkdir $workdir;
      
    $scriptname = "$workdir/matches.sh";
    open(SCRIPT, ">$scriptname") or die255("Could not create '$scriptname'.");
    print SCRIPT "#\!/usr/bin/env bash\n\n";
    print SCRIPT "perl $pwd/$0 -m -D $exec_directory -f $pwd/$opt_f\n";
    close(SCRIPT);
    chmod 0755, $scriptname;
      
    $matchdir = $workdir;
} else {
    $workdir = $pwd;
}

# testsuite
open(TESTSUITE, "<".$opt_f ) or die255("Cannot open testsuite file '$opt_f'.");

while ($_ = <TESTSUITE>) {

    # remove trailing newline 
    chomp; 
    # remove leading whitespace 
    $_ =~ s/^\s+//;
    # remove trailing whitespace
    $_ =~ s/\s+$//;

    # skip blank lines
    next if (length($_) == 0);

    # skip comments
    next if /^#/;

    if ( $_ =~ /^Test\s*:\s*(.*)\s*$/) {
      $test{"name"} = $1;
      if($test{"name"} eq "") {
	  die255("No name was provided with Test tag.");
      }
      print "$color_start{blue} ***** $test{\"name\"} ***** $color_end{blue} \n\n";
      print "Using workdir    : $workdir\n";
      if($opt_p) {
	  print "Workdir will be saved.\n";
      }
      print "Using test file  : $opt_f \n";

    } elsif ( $_ =~ /^Enabled\s*:\s*(.*)\s*$/) {
      %test = ();
      $enabled = $1;
      $enabled =~ s/^\s*//;
      $enabled =~ s/\s*$//;
      $test{"enabled"} = $enabled;

      if ( $enabled eq "No") {
          print STDERR "Test disabled: skipping test\n\n";
	  skip_exit();
      } elsif ( $enabled ne "Yes") {
	  if (!$opt_p && !$opt_m) { system ("rm -rf $workdir"); }
	  die255("Unknown option 'Enabled = $enabled' in testsuite file.");
      }

    } elsif ( $_ =~ /^Options\s*:\s*(.*)\s*$/) {
        $options_required = $1;
	# note: we could implement Options by baking this into the script via configure...
	
    } elsif ( $_ =~ /^Options_MPI\s*:\s*(.*)\s*$/) {
        if ($is_parallel && $np ne "serial") {
	    $options_required_mpi = $1;
	    $options_are_mpi = 1;
	}

    } elsif ( $_ =~ /^Program\s*:\s*(.*)\s*$/) {
	$command = "$exec_directory/$1";

	# FIXME: should we do this for a dry-run?
	if( ! -x "$command") {
	  $command = "$exec_directory/../utils/$1";
	}	
	if( ! -x $command) {
	    die255("Executable '$1' not available.");
        }

	$options_available = `$command -c`;
	chomp($options_available);
	if($is_parallel && $options_available !~ "mpi") {
	    print "Running in serial since executable was not compiled with MPI.\n";
	    $is_parallel = 0;
	    $options_are_mpi = 0;
	}

	if($options_are_mpi) {
	    $options_required = $options_required_mpi;
	}

	if(length($options_required) > 0) {
	    # check if the executable was compiled with the required options
	    foreach my $option (split(/;/, $options_required)){
		if(" $options_available " !~ " $option ") {
		    print "\nSkipping test: executable does not have the required option '$option'";
		    if($options_are_mpi) {
			print " for MPI";
		    }
		    print ".\n";
		    print "Executable: $command\n";
		    print "Available options: $options_available\n\n";
		    skip_exit();
		}
	    }
	}
	# FIXME: import Options to BGW version
    } elsif ( $_ =~ /^TestGroups\s*:/) {
        # handled by oct-run_testsuite.sh
    } else {
      if ( $enabled eq "") {
	die255("Testsuite file must set Enabled tag before another (except Test, Program, Options, TestGroups).");
      }

      if ( $_ =~ /^Util\s*:\s*(.*)\s*$/) {
	$np = "serial";
	$command = "$exec_directory/$1";
	if( ! -x "$command") {
	  $command = "$exec_directory/../utils/$1";
	}

	if( ! -x "$command") {
	    die255("Cannot find utility '$1'.");
	}
      }

      elsif ( $_ =~ /^Processors\s*:\s*(.*)\s*$/) {
	  # FIXME: enforce this is "serial" or numeric
	  $np = $1;
      }

      elsif ( $_ =~ /^Input\s*:\s*(.*)\s*$/) {

        $input_base = $1;
        $input_file = dirname($opt_f) . "/" . $input_base;
      
	if ( $opt_m ) {
	    print "\n\nFor input file : $input_file\n\n";
	    $return_value = 0;
	    # FIXME: this works from outer directory, but not in archived subdirectories.
	    $matchdir = "$workdir/$input_base";
	} else {
          if( -f $input_file ) {
            print "\nUsing input file : $input_file\n";
            $cp_return = system("cp $input_file $workdir/inp");
	    if($cp_return != 0) {
		die255("Copy failed (cp $input_file $workdir/inp)\n");
	    }
            # Ensure that the input file is writable so that it can
            # be overwritten by the next test.
            $mode = (stat "$workdir/inp")[2];
            chmod $mode|S_IWUSR, "$workdir/inp";
          } else {
            die255("Could not find input file '$input_file'.");
          }
      
	  # serial or MPI run?
	  if ( $is_parallel && $np ne "serial") {
            if("$global_np" ne "") {
		$np = $global_np;
            }
	    if ("$mpiexec" =~ /ibrun/) { # used by SGE parallel environment
		$specify_np = "";
		$my_nslots = "MY_NSLOTS=$np";
	    } elsif ("$mpiexec" =~ /runjob/) { # used by BlueGene
		$specify_np = "--np $np --exe";
		$my_nslots = "";
	    } elsif ("$mpiexec" =~ /poe/) { # used by IBM PE 
		$specify_np = ""; 
		$my_nslots = "MP_PROCS=$np"; 
	    } else { # for mpirun and Cray's aprun
		$specify_np = "-n $np";
		$my_nslots = "";
	    }
	    $command_line = "cd $workdir; $my_nslots $mpiexec $specify_np $machinelist $aexec $command inp > out";
	  } else {
	      $command_line = "cd $workdir; $aexec $command inp > out ";
	  }

	  # MPI implementations generally permit using more tasks than actual cores, and running tests this way makes it likely for developers to find race conditions.
	  if($np ne "serial") {
	      if($np > 4) {
		  print "Note: this run calls for more than the standard maximum of 4 MPI tasks.\n";
	      }
	  }

	  print "Executing: " . $command_line . "\n";

	  if ( !$opt_n ) {
	    $test_start = [gettimeofday];
	    $return_value = system("$command_line");
	    $test_end   = [gettimeofday];

	    $elapsed = tv_interval($test_start, $test_end);
	    printf("\tElapsed time: %8.1f s\n\n", $elapsed);

	    if($return_value == 0) {
	      printf "%-40s%s", " Execution", ": \t [ $color_start{green}  OK  $color_end{green} ] \n";
	      
	    } else {
	      print "Test run failed with exit code $return_value.\n";
	      print "These are the last lines of output:\n\n";
	      print "----------------------------------------\n";
	      system("tail -20 $workdir/out");
	      print "----------------------------------------\n\n";

	      printf "%-40s%s", " Execution", ": \t [ $color_start{red} FAIL $color_end{red} ] \n\n";

	      $failures++;
	      $test_succeeded = 0;  
	    }
	    $test{"run"} = 1;
	  }

	  # copy all files of this run to archive directory with the name of the
	  # current input file
	  mkdir "$workdir/$input_base";
	  @wfiles = `ls -d $workdir/* | grep -v inp`;
	  $workfiles = join("",@wfiles);
	  $workfiles =~ s/\n/ /g;
	  $cp_return = system("cp -r $workfiles $workdir/inp $workdir/$input_base");
	  if($cp_return != 0) {
	      die255("Copy failed (cp -r $workfiles $workdir/inp $workdir/$input_base)\n");
	  }
	}
      }

      elsif ( $_ =~ /^Precision\s*:\s*(.*)\s*$/) {
	set_precision($1);
      }

      elsif ( $_ =~ /^ExtraFile\s*:\s*(.*)\s*$/) {
	$file_cp = dirname($opt_f)."/".$1;
	$cp_return = system("cp $file_cp $workdir/");
      }

      elsif ( $_ =~ /^match/ ) {
	  # FIXME: should we do matches even when execution failed?
	  if (!$opt_n && $return_value == 0) {
	      if(run_match_new($_)){
		  printf "%-40s%s", "$name", ":\t [ $color_start{green}  OK  $color_end{green} ] \t (Calculated value = $value) \n";
		  if ($opt_v) { print_hline(); }
	      } else {
		  printf "%-40s%s", "$name", ":\t [ $color_start{red} FAIL $color_end{red} ] \n";
		  print_hline();
		  $test_succeeded = 0;
		  $failures++;
	      }
	  }
      } else {
	  die255("Unknown command '$_'.");
      }
    }

}

if ($opt_l && !$opt_m && !$opt_n)  { system ("cat $workdir/out >> out.log"); }
if (!$opt_p && !$opt_m && $test_succeeded) { system ("rm -rf $workdir"); }

print "\n";
close(TESTSUITE);

print "Status: ".$failures." failures\n";

exit $failures;


sub run_match_new {
  die255("Have to run before matching.") if !$test{"run"} && !opt_m;

  # parse match line
  my ($line, $match, $match_command, $shell_command, $ref_value, $off);
  $line = $_[0];
  $line =~ s/\\;/_COLUMN_/g;
  ($match, $name, $match_command, $ref_value) = split(/;/, $line);
  $match_command =~ s/_COLUMN_/;/g;
  $ref_value =~ s/^\s*//;
  $ref_value =~ s/\s*$//;

  # parse command
  $match_command =~ /\s*(\w+)\s*\((.*)\)/;

  my $func = $1;
  my $params = $2;

  # parse parameters
  $params =~ s/\\,/_COMMA_/g;
  my @par = split(/,/, $params);
  for($params=0; $params <= $#par; $params++){
    $par[$params] =~ s/_COMMA_/,/g;
    $par[$params] =~ s/^\s*//;
    $par[$params] =~ s/\s*$//;
  }

  if($func eq "SHELL"){ # function SHELL(shell code)
    check_num_args(1, 1, $#par, $func);
    $shell_command = $par[0];

  }elsif($func eq "LINE") { # function LINE(filename, line, column)
    check_num_args(3, 3, $#par, $func);
    if($par[1] < 0) { # negative number means from end of file
      $line_num = "`wc -l $par[0] | awk '{print \$1}'`";
      $shell_command = "awk -v n=$line_num '(NR==n+$par[1]+1)' $par[0]";
    } else {
      $shell_command = "awk '(NR==$par[1])' $par[0]";
    }
    $shell_command .= " | cut -b $par[2]-";

  }elsif($func eq "LINEFIELD") { # function LINE(filename, line, field)
    check_num_args(3, 3, $#par, $func);
    if($par[1] < 0) { # negative number means from end of file
      $line_num = "`wc -l $par[0] | awk '{print \$1}'`";
      $shell_command = "awk -v n=$line_num '(NR==n+$par[1]+1) {printf \$$par[2]}' $par[0]";
    } else {
      $shell_command = "awk '(NR==$par[1]) {printf \$$par[2]}' $par[0]";
    }

  }elsif($func eq "GREP") { # function GREP(filename, 're', column <, [offset>])
    check_num_args(3, 4, $#par, $func);
    if($#par == 3) {
	$off = $par[3];
    } else {
	$off = 0;
    }
    # -a means even if the file is considered binary due to a stray funny character, it will work
    $shell_command = "grep -a -A$off $par[1] $par[0] | awk '(NR==$off+1)'";
    $shell_command .= " | cut -b $par[2]-";

  }elsif($func eq "GREPFIELD") { # function GREPFIELD(filename, 're', field <, [offset>])
    check_num_args(3, 4, $#par, $func);
    if($#par == 3) {
	$off = $par[3];
    } else {
	$off = 0;
    }
    # -a means even if the file is considered binary due to a stray funny character, it will work
    $shell_command = "grep -a -A$off $par[1] $par[0]";
    $shell_command .= " | awk '(NR==$off+1) {printf \$$par[2]}'";
    # if there are multiple occurrences found by grep, we will only be taking the first one via awk

  }elsif($func eq "SIZE") { # function SIZE(filename)
    check_num_args(1, 1, $#par, $func);
    $shell_command = "ls -lt $par[0] | awk '{printf \$5}'";

  }else{ # error
    printf STDERR "ERROR: Unknown command '$func'\n";
    return 0;
  }

  # 'set -e; set -o pipefail' (bash 3 only) would make the whole pipe series give an error if any step does;
  # otherwise the error comes only if the last step failed.
  $value = qx(cd $matchdir && $shell_command);
  # Perl gives error code shifted, for some reason.
  $exit_code = $? >> 8;
  if($exit_code) {
      print STDERR "ERROR: Match command failed: $shell_command\n";
      return 0;
  }

  # extract numeric string (including possibility of NaN)
  if($value =~ /([0-9\-+.eEdDnNaA]+)/) {
      $value = $1;
      chomp $value;
  } else {
      $value = "";
  }

  if(length($value) == 0) {
      print STDERR "ERROR: Match command returned nothing: $shell_command\n";
      return 0;
  }

  if(!looks_like_number($value)) {
      print STDERR "ERROR: Match command returned non-numeric value '$value': $shell_command\n";
      return 0;
  }

  if(!looks_like_number($ref_value)) {
      print STDERR "WARNING: Match command has non-numeric reference value '$value'.\n";
      return 0;
  }

  # at this point, we know that the command was successful, and returned a number.
  $success = (abs(($value)-($ref_value)) <= $precnum);

  if(!$success || $opt_v) {
    print_hline();
    print "Match".$name.":\n\n";
    print "   Calculated value : ".$value."\n";
    print "   Reference value  : ".$ref_value."\n";
    print "   Difference       : ".abs($ref_value - $value)."\n";
    print "   Tolerance        : ".$precnum."\n\n";
  }

  return $success;
}

sub print_hline {
  print "\n-----------------------------------------\n\n";
}

# return value of environment variable (specified by string argument), or "" if not set
sub get_env {
    if(exists($ENV{$_[0]})) {
	return $ENV{$_[0]};
    } else {
	return "";
    }
}

# args: min num args, max num args, args given, function name
sub check_num_args {
    my $min_num_args   = $_[0];
    my $max_num_args   = $_[1];
    my $given_num_args = $_[2]+1;
    my $func_name      = $_[3];

    if($given_num_args < $min_num_args) {
	die255("$func_name given $given_num_args argument(s) but needs at least $min_num_args.");
    }
    if($given_num_args > $max_num_args) {
	die255("$func_name given $given_num_args argument(s) but can take no more than $max_num_args.");
    }
}

sub die255 {
    print STDERR "ERROR: " . $_[0] . "\n";
    print "Status: error\n";
    exit 255;
}

sub skip_exit {
    if (!$opt_p && !$opt_m && $test_succeeded) { system ("rm -rf $workdir"); }
    if($failures == 0) {
	print "Status: skipped\n";
	exit 254
    } else {
	print "Status: ".$failures." failures\n";
	exit $failures;
	# if a previous step has failed, mark as failed not skipped
    }
}
