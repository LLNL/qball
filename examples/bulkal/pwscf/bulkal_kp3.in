 &control
    calculation='scf'
    restart_mode='from_scratch',
    pseudo_dir = '.',
    outdir='WF/',
    prefix='al_lda'
    tprnfor = .true.
    verbosity = 'high'
    tstress = .true.
 /
 &system    
    ibrav=  2, celldm(1) =7.60, nat=  1, ntyp= 1,
    ecutwfc =20.0, occupations='smearing', smearing='fermi-dirac',
    degauss=0.05, nbnd=6
 /
 &electrons
    diagonalization='cg'
    conv_thr = 1d-10
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Al  26.98 TM_CAPW91_Al.UPF
ATOMIC_POSITIONS
 Al 0.00 0.00 0.00 
K_POINTS automatic 
 3 3 3 0 0 0
