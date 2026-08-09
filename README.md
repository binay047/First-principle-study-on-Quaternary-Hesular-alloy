# First-Principles Study of BeScCoSi 
*This repository contains the input files, scripts, and calculation workflow used to investigate the structural, electronic, magnetic, dynamical, mechanical, piezoelectric, thermodynamic, optical, and thermoelectric properties of BeScCoSi and MgScCoSi using Density Functional Theory (DFT).*
*The calculations were mainly performed using Quantum ESPRESSO, thermo_pw, and BoltzTraP2.*
# 1. Phase Stability
*The total energies of the three phases were calculated and compared.*
*The phase with the lowest total energy was selected as the most stable structure.*
*Phase III, having the LiMgPdSn-type cubic structure with space group F-43m (No. 216), was found to be the most stable phase for both compounds.*

# 2. Convergence Tests
* chmod +x ecut.sh
*./ecut.sh
* **Note:**  Take a converged ecut and place it in k.sh
 * chmod +x k.sh
 * ./k.sh
* **Note:**  Take a converged kpoint and place it in lattice.sh
  * chmod +x lattice.sh
  * ./lattice.sh
* **Note:**  Now, you will get etot_vs_k.dat
  * **Note:** Open terminal in the same lattice directory
* ev.x 
* ang
* noncubic
* 4
* input file name: lattice.dat
* output file name: bin
* **Note:**  Take this a0 from bin into vc_relax.in
# 3. relaxation
*pw.x <vc_relx.in> vc_relax.out
* **Note:** convert cell_parameters into format of a and update "a" and "atomic_position" in new scf.in

# 3. In Dos directory
* mpirun -np 8 pw.x <scf.in> scf.out
* mpirun -np 8 pw.x <nscf.in> nscf.out
* dos.x <dos.in> dos.out
* awk 'NR>1 {print $1, $2}' dos.dat > raw_up.dat
* awk 'NR>1 {print $1, -$3}' dos.dat > raw_down.dat
* projwfc.x < pdos.in > pdos.out
* sumpdos.x *\(Be\)* > atom_Be_tot.dat
* sumpdos.x *\(Mg\)* > atom_Mg_tot.dat
* sumpdos.x *\(Co\)* > atom_Co_tot.dat
* sumpdos.x *\(Si\)* > atom_Si_tot.dat

# 4. In band directory
* mpirun -np 8 pw.x <scf.in> scf. out
* mpirun -np 8 pw.x <band.in> band.out
* **Note:**  note: kpoints in band. in is generated using xcrysden 
* bands.x <bands.in> bands.out
* plot bands_plot.bands.gnu file using xmgrace

# 5. In phonon directory
* mpirun -np 8 pw.x <scf.in> scf.out
* mpirun -np 8 ph.x <ph.in > ph.out
* mpirun -np 8 q2r.x <q2r.in > q2r.out
* mpirun -np 8 matdyn.x <matdyn.in > matdyn.out
* plotband.x <plotband.in> plotband.out
* matdyn.x <phdos.in> phdos.out
* awk '{print $1,$2}' phdos.dat> total.dat
* awk '{print $1,$3}' phdos.dat> Be.dat
* awk '{print $1,$4}' phdos.dat> Sc.dat
* awk '{print $1,$5}'  phdos.dat> Co.dat
* awk '{print $1,$6}'  phdos.dat> Si.dat

# 6. thermo directory
* **Note:** Please create an empty out folder in the working directory, and we need scf. in and thermo_control files, finally run
* mpirun -np 8 thermo_pw.x <scf.in> scf. out
** Now, to extract specific heat capacity, free energy and entropy.dat, use the following awk commands inside the therm_files folder generated after running the above code
* awk 'BEGIN{print "#T(K)   Cv(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$5*1312749.8}' output_therm.dat_debye.g1 > Cv.dat
* awk 'BEGIN{print "#T(K)   Free_Energy(KJmol^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,($3*1312749.8)/1000}' output_therm.dat_debye.g1 > FreeEnergy.dat
* awk 'BEGIN{print "#T(K)   Entropy(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$4*1312749.8}' output_therm.dat_debye.g1 > entropy.dat



