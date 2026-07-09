#bin/sh/     
                                                                                                                                               
NAME="Ecutwfc"
for ecut in 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 ; do
cat > "${NAME}_${ecut}.in" <<EOF 

 &CONTROL
  calculation = 'scf'
  outdir = './out/'
  prefix = 'aiida'
  wf_collect = .true.
  pseudo_dir = './pseudo/'
  verbosity = 'high'
/
&SYSTEM
  ecutrho =   $((10*ecut))
  ecutwfc =   $ecut
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
8 8 8 0 0 0

EOF
  mpirun -np 2 pw.x -inp "${NAME}_${ecut}.in" | tee "${NAME}_${ecut}.out"   
#   pw.x <"${NAME}_${ecut}.in"> "${NAME}_${ecut}.out"
    echo "${NAME}_${ecut}"
    grep "!" "${NAME}_${ecut}.out"

 # Write cut-off and total energies in calcecut.dat.                       
                                                                                  
    awk '/!/ {printf "%d %s\n", ('$ecut'), $5}' "${NAME}_${ecut}.out" >>"ecutwfc.dat"  
done

