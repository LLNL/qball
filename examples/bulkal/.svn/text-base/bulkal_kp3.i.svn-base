set nrowmax 2
bulkal_kp3.sys
set ecut 20.0
set wf_dyn PSDA
set nempty 4
set xc LDA

set smearing fermi
set smearing_width 0.05

set ecutprec 5.0

symmetry  1  0  0  0  1  0  0  0  1
symmetry  0  1  -1  1  0  -1  0  0  -1
symmetry  -1  0  0  -1  0  1  -1  1  0
symmetry  0  -1  1  0  -1  0  1  -1  0
symmetry  0  -1  0  -1  0  0  0  0  -1
symmetry  -1  0  1  0  -1  1  0  0  1
symmetry  0  1  0  0  1  -1  -1  1  0
symmetry  1  0  -1  1  0  0  1  -1  0
symmetry  -1  0  0  0  0  -1  0  -1  0
symmetry  1  0  0  1  -1  0  1  0  -1
symmetry  0  1  -1  -1  1  0  0  1  0
symmetry  0  -1  1  0  0  1  -1  0  1
symmetry  -1  1  0  0  1  0  0  1  -1
symmetry  0  0  -1  0  -1  0  -1  0  0
symmetry  0  0  1  -1  0  1  0  -1  1
symmetry  1  -1  0  1  0  -1  1  0  0
symmetry  -1  0  1  -1  1  0  -1  0  0
symmetry  0  1  0  0  0  1  1  0  0
symmetry  1  0  -1  0  0  -1  0  1  -1
symmetry  0  -1  0  1  -1  0  0  -1  1
symmetry  0  0  -1  0  1  -1  1  0  -1
symmetry  -1  1  0  -1  0  0  -1  0  1
symmetry  0  0  1  1  0  0  0  1  0
symmetry  1  -1  0  0  -1  1  0  -1  0
symmetry  -1  0  0  0  -1  0  0  0  -1
symmetry  0  -1  1  -1  0  1  0  0  1
symmetry  1  0  0  1  0  -1  1  -1  0
symmetry  0  1  -1  0  1  0  -1  1  0
symmetry  0  1  0  1  0  0  0  0  1
symmetry  1  0  -1  0  1  -1  0  0  -1
symmetry  0  -1  0  0  -1  1  1  -1  0
symmetry  -1  0  1  -1  0  0  -1  1  0
symmetry  1  0  0  0  0  1  0  1  0
symmetry  -1  0  0  -1  1  0  -1  0  1
symmetry  0  -1  1  1  -1  0  0  -1  0
symmetry  0  1  -1  0  0  -1  1  0  -1
symmetry  1  -1  0  0  -1  0  0  -1  1
symmetry  0  0  1  0  1  0  1  0  0
symmetry  0  0  -1  1  0  -1  0  1  -1
symmetry  -1  1  0  -1  0  1  -1  0  0
symmetry  1  0  -1  1  -1  0  1  0  0
symmetry  0  -1  0  0  0  -1  -1  0  0
symmetry  -1  0  1  0  0  1  0  -1  1
symmetry  0  1  0  -1  1  0  0  1  -1
symmetry  0  0  1  0  -1  1  -1  0  1
symmetry  1  -1  0  1  0  0  1  0  -1
symmetry  0  0  -1  -1  0  0  0  -1  0
symmetry  -1  1  0  0  1  -1  0  1  0

set nkpoints 4
kpoint   0.0000000   0.0000000   0.0000000   0.0740741  crystal
kpoint   0.0000000   0.0000000   0.3333333   0.5925926  crystal
kpoint   0.0000000   0.3333333   0.3333333   0.4444444  crystal
kpoint   0.0000000   0.3333333  -0.3333333   0.8888889  crystal

set threshold_scf 1.E-10 50

randomize_wf
run 0 400 2

quit
