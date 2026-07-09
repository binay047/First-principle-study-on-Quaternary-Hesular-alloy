#!/bin/sh
rm -f k.out etot_vs_k.dat
touch etot_vs_k.dat
for k in 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 6.1 6.2 6.3 6.4;do
cat > k.in<<EOF

  
  &CONTROL
  calculation = 'scf'
  outdir = './out/'
  prefix = 'aiida'
  pseudo_dir = './pseudo/'
  verbosity = 'high'
/
&SYSTEM
  ecutrho = 550   
  ecutwfc = 55
  ibrav = 2
  a = $k
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
7 7 7 0 0 0

EOF
mpirun -np 2 pw.x <k.in > k.out
# extract Etot from output
etot=`grep -e ! k.out | awk '{print $(NF-1)}'`
echo $k $etot  >> etot_vs_k.dat
done
