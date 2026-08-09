# Thermoelectric-on-BeScCoSi
* In BeScCoSi, you have three phases

#In Dos directory
mpirun -np 8 pw.x <scf.in> scf.out
mpirun -np 8 pw.x <nscf.in> nscf.out
dos.x <dos.in> dos.out
awk 'NR>1 {print $1, $2}' dos.dat > raw_up.dat
awk 'NR>1 {print $1, -$3}' dos.dat > raw_down.dat
projwfc.x < pdos.in > pdos.out
sumpdos.x *\(Be\)* > atom_Be_tot.dat
sumpdos.x *\(Mg\)* > atom_Mg_tot.dat
sumpdos.x *\(Co\)* > atom_Co_tot.dat
sumpdos.x *\(Si\)* > atom_Si_tot.dat


# In band directory
mpirun -np 8 pw.x <scf.in> scf. out
mpirun -np 8 pw.x <band.in> band.out
bands.x <bands.in> bands.out
plot bands_plot.bands.gnu file using xmgrace

#In phonon directory
mpirun -np 8 pw.x <scf.in> scf.out
mpirun -np 8 ph.x <ph.in > ph.out
mpirun -np 8 q2r.x <q2r.in > q2r.out
mpirun -np 8 matdyn.x <matdyn.in > matdyn.out
plotband.x <plotband.in> plotband.out
matdyn.x <phdos.in> phdos.out
awk '{print $1,$2}' phdos.dat> total.dat
awk '{print $1,$3}' phdos.dat> Be.dat
awk '{print $1,$4}' phdos.dat> Sc.dat
awk '{print $1,$5}'  phdos.dat> Co.dat
awk '{print $1,$6}'  phdos.dat> Si.dat

# In the thermo directory using the Slack model,
!we need to create an empty out folder in the working directory, and we need scf. in and thermo_control files, finally run
mpirun -np 8 thermo_pw.x <scf.in> scf. out
! Now, to extract specific heat capacity, free energy and entropy.dat, use the following awk commands inside the therm_files folder generated after running the above code  
awk 'BEGIN{print "#T(K)   Cv(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$5*1312749.8}' output_therm.dat_debye.g1 > Cv.dat
awk 'BEGIN{print "#T(K)   Free_Energy(KJmol^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,($3*1312749.8)/1000}' output_therm.dat_debye.g1 > FreeEnergy.dat
awk 'BEGIN{print "#T(K)   Entropy(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$4*1312749.8}' output_therm.dat_debye.g1 > entropy.dat



