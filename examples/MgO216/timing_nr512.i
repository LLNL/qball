set nrowmax 512
set force_complex_wf ON
mg108o108h1.sys
set ecut 50
set ecutprec 10
set nempty 11

set wf_dyn PSDA
set ecutprec 8.0

randomize_wf
run 0 5 2


quit
