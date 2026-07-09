#!/bin/sh
rm -f k.out etot_vs_k.dat
touch etot_vs_k.dat
for k in 2 3 4 5 6 7 8 9 10 11 12 13 ;do
cat > k.in<<EOF
 &CONTROL
  calculation = 'scf'
  outdir = './out/'
  prefix = 'aiida'
  wf_collect = .true.
  pseudo_dir = './pseudo/'
  verbosity = 'high'
/
&SYSTEM
  ecutrho = 550   
  ecutwfc = 55
  ibrav = 2
  a = 5.791172
  nat = 4
  nosym = .false.
  ntyp = 4
  occupations = 'smearing'
  smearing = 'gaussian'
  degauss = 0.001
/ 
&ELECTRONS
  conv_thr = 1.0d-7
  electron_maxstep = 800
  mixing_beta =   4.0000000000d-01
/
ATOMIC_SPECIES
Be     9.012182 Be.pbe-n-rrkjus_psl.1.0.0.UPF
Sc     44.955912 Sc.pbe-spn-rrkjus_psl.1.0.0.UPF
Co     58.933195 Co.pbe-spn-rrkjus_psl.0.3.1.UPF
Si     28.0855 Si.pbe-nl-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
Be           0.7500000000       0.7500000000       0.7500000000
Sc           0.0000000000       0.0000000000       0.0000000000
Co           0.2500000000       0.2500000000       0.2500000000
Si           0.5000000000       0.5000000000       0.5000000000

K_POINTS automatic
$k $k $k 0 0 0

EOF
mpirun -np 2 pw.x <k.in > k.out
# extract Etot from output
etot=`grep -e ! k.out | awk '{print $(NF-1)}'`
echo $k $etot  >> etot_vs_k.dat
done
