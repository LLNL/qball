#! /usr/bin/perl -w
#
# generate input files for large fcc gold systems for Qbox HPC benchmarking 
#
# written by Erik Draeger, LLNL, 4/12/2013

if ($#ARGV < 1) {
   print "syntax:  qbox_gold_makeinputs [number of atoms] [number of MPI tasks]\n";
  exit;
}

# if you're using the open source version from http://eslab.ucdavis.edu/software/qbox/ set this to one.
# if you're using the qb@LL version, set this to zero.
$useDavisVersion = 0;   


$nAtoms = $ARGV[0];
$nTasks = $ARGV[1];

# print coordinates to .sys file
$species = "gold";
$atomtype = "Au";
$pseudo = "Au_PBE.xml";
$a0 = 7.71;

$nx = 2;
$ny = 2;
if ($nAtoms >= 256)
{
   $nx = 4;
   $ny = 4;
}
$nz = int $nAtoms/(4*$nx*$ny);

if ($nz < 1) { $nz = 1; }

$nFinal = 4*$nx*$ny*$nz;

print "Using natoms = $nFinal\n";

$ax = $a0*$nx;
$ay = $a0*$ny;
$az = $a0*$nz;

$sysfile = join '',"gold.N",$nFinal,".sys";
open SYS, ">$sysfile";
print SYS "set cell $ax 0.0 0.0 0.0 $ay 0.0 0.0 0.0 $az\n";
print SYS "species $species $pseudo\n";

$atomcnt = 1;
for ($i=0; $i<$nx; $i++) { 
  for ($j=0; $j<$ny; $j++) { 
    for ($k=0; $k<$nz; $k++) { 

      $x = 0.0 + $i*$a0;
      $y = 0.0 + $j*$a0;
      $z = 0.0 + $k*$a0;
      printf SYS "atom $atomtype$atomcnt $species %11.7f %11.7f %11.7f\n",$x,$y,$z;
      $atomcnt++;

      $x = 0.5*$a0 + $i*$a0;
      $y = 0.5*$a0 + $j*$a0;
      $z = 0.0 + $k*$a0;
      printf SYS "atom $atomtype$atomcnt $species %11.7f %11.7f %11.7f\n",$x,$y,$z;
      $atomcnt++;

      $x = 0.5*$a0 + $i*$a0;
      $y = 0.0 + $j*$a0;
      $z = 0.5*$a0 + $k*$a0;
      printf SYS "atom $atomtype$atomcnt $species %11.7f %11.7f %11.7f\n",$x,$y,$z;
      $atomcnt++;

      $x = 0.0 + $i*$a0;
      $y = 0.5*$a0 + $j*$a0;
      $z = 0.5*$a0 + $k*$a0;
      printf SYS "atom $atomtype$atomcnt $species %11.7f %11.7f %11.7f\n",$x,$y,$z;
      $atomcnt++;

    }
  }
}
close SYS;

# print input file

$nrowmax = 1024;
if ($nTasks < 128) { $nrowmax = $nTasks; }
elsif ($nTasks <= 4096) { $nrowmax = $nTasks / 4; }
elsif ($nTasks >= 65536) { $nrowmax = 2048; }

$ecut = 130;                     # plane wave basis size, defined by ecut in Rydbergs
$nvalence = 17;                  # number of valence electrons in the pseudopotential
$nstates = $nFinal*$nvalence/2;  # total number of electronic orbitals
$nempty = int 0.2*$nstates;      # number of unoccupied orbitals


$inputfile = join '',"gold.N",$nFinal,".i";
open IN, ">$inputfile";

print IN "set nrowmax $nrowmax\n";
print IN "$sysfile\n";
print IN "\n";
print IN "set xc PBE\n";
print IN "set ecut $ecut\n";
print IN "set nempty $nempty\n";
print IN "set fermi_temp 500.0\n";
print IN "set wf_dyn PSDA\n";
print IN "set ecutprec 8.0\n";
print IN "\n";
print IN "set charge_mix_rcut 3.0\n";
print IN "set charge_mix_coeff 0.8\n";
print IN "\n";
if ($useDavisVersion == 1)
{
   print IN "kpoint del 0 0 0\n";
   print IN "kpoint add 0.0000001 0 0 1.0\n";
}
else
{
   print IN "kpoint 0.0000001 0 0 1.0\n";
}
print IN "\n";
print IN "randomize_wf\n";
print IN "run 0 1 5\n";
print IN "\n";
print IN "quit\n";
close IN;

