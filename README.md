# First-Principles Study of BeScCoSi 
*This repository contains the input files, scripts, and calculation workflow used to investigate the structural, electronic, magnetic, dynamical, mechanical, piezoelectric, thermodynamic, optical, and thermoelectric properties of BeScCoSi using Density Functional Theory (DFT).*
*The calculations were mainly performed using Quantum ESPRESSO, thermo_pw, and BoltzTraP2.*
# 1. Phase Stability
*The total energies of the three phases were calculated and compared.*
*The phase with the lowest total energy was selected as the most stable structure.*
*Phase III, having the LiMgPdSn-type cubic structure with space group F-43m (No. 216), was found to be the most stable phase.*

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
  * **Note:** Open a terminal in the same lattice directory
* ev.x 
* ang
* noncubic
* 4
* input file name: lattice.dat
* output file name: bin
* **Note:**  Take this a0 from bin into vc_relax.in
*pw.x <vc_relx.in> vc_relax.out
* **Note:** convert cell_parameters into the format of a and update "a" and "atomic_position" in new scf.in

# 2. In Dos directory
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

# 3. In band directory
* mpirun -np 8 pw.x <scf.in> scf. out
* mpirun -np 8 pw.x <band.in> band.out
* **Note:**  note: kpoints in band. in is generated using xcrysden 
* bands.x <bands.in> bands.out
* plot bands_plot.bands.gnu file using xmgrace

# 4. In phonon directory
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

# 5. Thermo directory
* **Note:** Please create an empty out folder in the working directory, and we need scf. in and thermo_control files, finally run
* mpirun -np 8 thermo_pw.x <scf.in> scf. out
* **Note:**  Now, to extract specific heat capacity, free energy and entropy.dat, use the following awk commands inside the therm_files folder generated after running the above code
* awk 'BEGIN{print "#T(K)   Cv(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$5*1312749.8}' output_therm.dat_debye.g1 > Cv.dat
* awk 'BEGIN{print "#T(K)   Free_Energy(KJmol^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,($3*1312749.8)/1000}' output_therm.dat_debye.g1 > FreeEnergy.dat
* awk 'BEGIN{print "#T(K)   Entropy(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$4*1312749.8}' output_therm.dat_debye.g1 > entropy.dat
* **Note:** now update in kl.py C11, C12, C44 from scf. out in GPa by dividing each of them by 10, rho from vc_relax.out by searching near final bfgs, Mavg adding all elements mass from scf.out, Omega_cell from vc_relax.out, at last of kl.py update Tmax(temperature upto which you want to calculate Zt)
* python3 kl.py

# 6. Optical directory
* **Note:** you need non-conserving pseudopotentials for optical properties calculation and add noinv = .true. in the system card in scf. in and nscf.in
* pw.x <scf.in> scf.out
* pw.x <nscf.in> nscf.out
* epsilon.x <epsilon.in> epsilon.out
* awk '{if(FNR<=2){ if(FNR==1) print "# Energy [eV]  Isotropic_Real_Dielectric_Function"; next}eps1=($2+$3+$4)/3;printf " %11.9f%11.9f\n",$1,eps1}' epsr_aiida.dat > dielectric_real_isotropic.dat
* awk '{if(FNR<=2){ if(FNR==1) print "# Energy [eV]  Isotropic_Imaginary_Dielectric_Function"; next} eps2=($2+$3+$4)/3; printf " %11.9f%11.9f\n",$1,eps2}' epsi_aiida.dat >
  dielectric_imaginary_isotropic.dat
* awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Reflectivity [fraction]"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2);  n=sqrt((mod+r[FNR])/2); k=sqrt((mod-r[FNR])/2); R=((n-1)^2 + k^2)/((n+1)^2 + k^2); printf "    %11.9f    %11.9f\n", $1, R }' epsr_aiida.dat epsi_aiida.dat > reflectivity_isotropic.dat
* awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Refractive_Index_n"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2); n=sqrt((mod+r[FNR])/2); printf "    %11.9f    %11.9f\n", $1, n }' epsr_aiida.dat epsi_aiida.dat > refractive_index_isotropic.dat
* awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Extinction_Coefficient_k"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2); k=sqrt((mod-r[FNR])/2); printf "    %11.9f    %11.9f\n", $1, k }' epsr_aiida.dat epsi_aiida.dat > extinction_coefficient_isotropic.dat
* awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Absorption_Coefficient [10^4/cm]"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2); k=sqrt((mod-r[FNR])/2); alpha=(2*$1*k*1.6231012e5)/10000; printf "    %11.9f    %14.6f\n", $1, alpha }' epsr_aiida.dat epsi_aiida.dat > absorption_isotropic_scaled.dat
* awk '{ if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Optical_Conductivity [10^3 Omega^-1 cm^-1]"; next } i_avg=($2+$3+$4)/3; sigma=(1327.21*$1*i_avg)/1000; printf "    %11.9f    %14.6f\n", $1, sigma }' epsi_aiida.dat > optical_conductivity_isotropic_scaled.dat
* awk 'NR==FNR { if(FNR>2) { r_avg=($2+$3+$4)/3; r[FNR]=r_avg } next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_EELS"; next } i_avg=($2+$3+$4)/3; loss = i_avg / (r[FNR]^2 + i_avg^2); printf "    %11.9f    %11.9f\n", $1, loss }' epsr_aiida.dat epsi_aiida.dat > energylossfunction_isotropic.dat

# 7. Raman spectroscopy
* pw.x <scf.in> scf. out
* ph.x <ph_raman.in> ph_raman.out
* **Note:** you need pz-hgh pesudopotentials
* python3 -c "import numpy as np; peak=[(183.62,13.0937), (280.46,251.2447),(401.29,308.4624)]; w=np.linspace(100, 500, 800); fit=sum(I0*np.exp(-((w-w0)/2)**2) for w0, I0 in peak); np.savetxt('raman_curve.dat', np.column_stack((w, fit)), fmt='%.4f')"
 * **Note:** from dynmat.out, look for double frequencies and replace them in the above (183.62,13.0937), (280.46,251.2447),(401.29,308.4624), 183.62, 280.46 and 401.29 are frequencies and 13.0937, 251.2447, 308.4624 are corresponding Raman values

# 8. Piezoelectricity 
* python3 piezo.py
* **Note:** for a different quaternary Heusler alloy, update in `piezo.py`:
  - `COMPOUND_LABEL`, `ALAT_ANGSTROM` — new compound name and relaxed lattice constant
  - `ATOMIC_SPECIES` — masses and pseudopotential filenames for the new elements
  - `ATOMIC_POSITIONS` — fractional coordinates for the DFT-confirmed site ordering (Type I/II/III) of the new compound; do not reuse the old dict as-is
  - `NBND` — recompute from the new total valence electron count
  - `ECUTWFC`, `ECUTRHO` — re-converge for the new pseudopotentials
  - If the new compound is magnetic: add `nspin = 2` and `starting_magnetization(ityp)` to both scf and nscf blocks
  - `GDIR = 1`, `KPTS_BASE`, `NPPSTR`, `STRAIN_VALUES`, and the yz-shear strain function stay unchanged as long as the compound is $F\bar43m$ ($T_d$), since $d_{14}=d_{25}=d_{36}$ still holds
 * pw.x <BeScCoSi_eta00_scf.in> BeScCoSi_eta00_scf.out
 * pw.x <BeScCoSi_eta00_nscf.in> BeScCoSi_eta00_nscf.out
   * **Note:** run all scf. in and nscf.in
  ## Theory: extracting e₁₄ from Berry-phase output (`extract_e14.py`)

### Background

For cubic $T_d$ ($\overline{4}3m$) symmetry, the piezoelectric tensor has a single independent component:

$$
d_{14} = d_{25} = d_{36}
$$

A pure shear strain $\eta_4$ (coupling $y$–$z$) induces a polarization that points **exactly along Cartesian $x$**, with $P_y = P_z = 0$ identically. This is the textbook symmetry argument that justifies why $d_{14}$ is the only independent constant for this point group. The stress constant $e_{14}$ is defined as the linear response:

$$
e_{14} = \left.\frac{dP_x}{d\eta_4}\right|_{\eta_4 = 0}
$$

and the strain constant follows from the elastic constant $C_{44}$:

$$
d_{14} = \frac{e_{14}}{C_{44}}
$$

### The non-Cartesian cell problem

The DFT cell uses `ibrav = 2` primitive FCC lattice vectors, which are **not aligned** with Cartesian $x/y/z$ — the first primitive vector points along

$$
\hat A_1 = \frac{1}{\sqrt2}(-1, 0, 1)
$$

not along $\hat x$. Running Berry-phase with `gdir = 1` therefore does not return $P_x$ directly — QE returns the polarization **projected onto** $\hat A_1$:

$$
P_{\text{measured}} = \vec P \cdot \hat A_1 = P_x \, d_x
$$

(using $P_y = P_z = 0$ from the symmetry argument above, so only the $P_x d_x$ term survives). QE prints $\hat A_1 = (d_x, d_y, d_z)$ explicitly as *"The polarization direction is: ..."*. The true $P_x$ is recovered by:

$$
P_x = \frac{P_{\text{measured}}}{d_x}
$$

computed **per strain point**, not with a fixed constant — straining the cell slightly rotates $\hat A_1$, so $d_x$ drifts point to point.

### What the script does

1. **Parses each `_nscf.out`** — regex-extracts `P = ... (mod ...) C/m^2` and the printed direction vector $(d_x, d_y, d_z)$.
2. **Applies the projection correction** above to recover $P_x$ for each strain point, using that file's own $d_x$.
3. **Sanity-checks each point:**
   - skips files with no polarization block (failed/crashed run)
   - skips points where $d_x \approx 0$ (correction would diverge)
   - flags points where $(d_x,d_y,d_z)$ deviates noticeably from the expected $(-0.707, 0, 0.707)$-type pattern
4. **Linear least-squares fit** (no numpy/scipy — closed-form slope/intercept) of

$$
P_x = e_{14}\,\eta_4 + P_0
$$

   reporting $R^2$; warns if $R^2 < 0.99$, since a poor fit usually signals a **branch jump** — Berry-phase polarization is only defined modulo a quantum $eR/\Omega$, so if QE's chosen branch jumps between adjacent strain points, the data won't lie on a clean line.
5. **Converts to $d_{14}$** using the user-supplied $C_{44}$ (GPa), with SI-to-pm/V unit conversion:

$$
d_{14}\ [\text{pm/V}] = \frac{e_{14}}{C_{44}} \times 10^{12}
$$

### Usage

```bash
python3 extract_e14.py <compound_label> <C44_GPa> [piezo_inputs_dir]
```

Example:

```bash
python3 extract_e14.py BeScCoSi 67.62 ./piezo_inputs
```

### Summary

`piezo.py` generates the strained structures and runs the scf → nscf(lberry) pipeline; `extract_e14.py` is the analysis half — it reads back the Berry-phase results, undoes the coordinate-system artifact from the non-Cartesian cell, fits the slope, and reports the final $d_{14}$.


  

  


