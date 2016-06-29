set nrowmax 16
sih4.sys
set ecut 36.0
set wf_dyn PSDA
set xc PBE

set ecutprec 8.0

set atoms_dyn MD
set thermostat SCALING
set th_time 4000
set th_temp 300
set dt 20

set threshold_scf 1.E-5 2

set savefreq 10
set savedenfreq 4

load -states wfsave/sih4
run 31 20
save -states mdtmp/mdchk
savesys mdtmp/mdchk.sys

quit
