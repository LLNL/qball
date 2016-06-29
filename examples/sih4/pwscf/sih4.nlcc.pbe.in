 &control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='sih4',
    tstress = .true.
    tprnfor = .true.
    verbosity = 'high'
    pseudo_dir = './',
    outdir='WF/'
 /
 &system    
    ibrav=  0, celldm(1) =14.00, nat=  5, ntyp= 2,
    ecutwfc =30.0, nbnd=4, nosym=.true.
 /
 &electrons
    diagonalization='cg'
    mixing_mode = 'plain'
    mixing_beta = 0.7 
    conv_thr =  1.0d-12
 /
CELL_PARAMETERS
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
ATOMIC_SPECIES
 Si 12.0107  Si.pbe-n-van.UPF
 H  1.000    H.pbe-rrkjus.UPF
ATOMIC_POSITIONS bohr
 Si   0.00000000   0.00000000   0.00000000 
 H  1.60000000   1.60000000   1.60000000 
 H  1.60000000  -1.60000000  -1.60000000 
 H -1.60000000   1.60000000  -1.60000000 
 H -1.60000000  -1.60000000   1.60000000 
K_POINTS crystal
1
 0.0 0.0 0.0 1.0
