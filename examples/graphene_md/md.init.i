set nrowmax 64
graphene.sys
set ecut 45
set XC LDA
set wf_dyn PSDA

set atoms_dyn MD
set center_of_mass fixed
set thermostat SCALING
set th_time 4000
set th_temp 300
set dt 20

set threshold_scf 1.E-5 3
set run_timer 11.5 hours

load -states wfsave/graphene.init
run 20 1000
save -states wfsave/graphene.md
savesys graphene.md.sys
run 20 1000
save -states wfsave/graphene.md
savesys graphene.md.sys
run 20 1000
save -states wfsave/graphene.md
savesys graphene.md.sys
run 20 1000
save -states wfsave/graphene.md
savesys graphene.md.sys
run 20 1000
save -states wfsave/graphene.md
savesys graphene.md.sys


quit
