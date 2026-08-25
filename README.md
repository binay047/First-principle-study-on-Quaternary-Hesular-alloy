# First-Principles Study of BeScCoSi

This repository contains the input files, scripts, and calculation workflow used to investigate the structural, electronic, magnetic, dynamical, mechanical, piezoelectric, thermodynamic, optical, and thermoelectric properties of BeScCoSi using Quantum espresso.

The calculations were mainly performed using Quantum ESPRESSO, thermo_pw, and BoltzTraP2.

---

## Table of Contents

1. [Phase Stability](#1-phase-stability)
2. [Convergence Tests](#2-convergence-tests)
3. [Density of States (DOS) and Projected DOS](#3-density-of-states-dos-and-projected-dos)
4. [Band Structure](#4-band-structure)
5. [Phonon Dispersion and Phonon DOS](#5-phonon-dispersion-and-phonon-dos)
6. [Thermodynamic Properties](#6-thermodynamic-properties)
7. [Optical Properties](#7-optical-properties)
8. [Raman Spectroscopy](#8-raman-spectroscopy)
9. [Piezoelectricity](#9-piezoelectricity)
10. [Formation Energy](#10-formation-energy)

---

## 1. Phase Stability

### Theory

Quaternary Heusler alloys of formula $XX'YZ$ can crystallize in several inequivalent Wyckoff-site orderings (commonly labeled Phase I, II, III), all sharing the same space group but differing in which element occupies which of the four sites ($4a$, $4b$, $4c$, $4d$). These configurations are not degenerate — the total energy

$$
E_{\text{tot}} = E_{\text{tot}}(\text{ordering})
$$

differs measurably between them because each ordering places different nearest-neighbor coordination environments around each atom. The thermodynamically stable phase is the one minimizing $E_{\text{tot}}$ at $T=0$.

### Procedure

- The total energies of the three phases were calculated and compared.
- The phase with the lowest total energy was selected as the most stable structure.
- Phase III, having the LiMgPdSn-type cubic structure with space group F-43m (No. 216), was found to be the most stable phase.

---

## 2. Convergence Tests

### Theory

DFT total energies computed with a plane-wave basis are only exact in the limit of infinite basis size and infinite k-point sampling; in practice both are truncated, so every observable must be checked for convergence before being trusted. The plane-wave cutoff $E_{\text{cut}}$ truncates the basis via

$$
\frac{\hbar^2 |\mathbf{G}+\mathbf{k}|^2}{2m_e} \leq E_{\text{cut}}
$$

and convergence is confirmed when $E_{\text{tot}}(E_{\text{cut}})$ changes by less than a chosen threshold (e.g. 1 mRy) upon further increase. Similarly, the Brillouin-zone integral for the charge density is approximated by a finite Monkhorst-Pack k-point mesh, and must be checked until $E_{\text{tot}}$ stabilizes with respect to mesh density. Once $E_{\text{cut}}$ and the k-mesh are fixed, the equilibrium lattice constant $a_0$ is obtained by fitting $E_{\text{tot}}(V)$ to an equation of state (e.g. Birch-Murnaghan), which is what `ev.x` performs on `etot_vs_k.dat` (despite the filename, this is the volume/lattice series, not the k-convergence series).

### Procedure

\`\`\`bash
chmod +x ecut.sh
./ecut.sh
\`\`\`

> **Note:** Take a converged `ecut` and place it in `k.sh`.

\`\`\`bash
chmod +x k.sh
./k.sh
\`\`\`

> **Note:** Take a converged `kpoint` and place it in `lattice.sh`.

\`\`\`bash
chmod +x lattice.sh
./lattice.sh
\`\`\`

> **Note:** You will now have `etot_vs_k.dat`. Open a terminal in the same `lattice` directory and run:

\`\`\`bash
ev.x
\`\`\`

At the prompts, enter:

\`\`\`
ang
noncubic
4
input file name: lattice.dat
output file name: bin
\`\`\`

> **Note:** Take the resulting $a_0$ from `bin` into `vc_relax.in`.

\`\`\`bash
pw.x < vc_relax.in > vc_relax.out
\`\`\`

> **Note:** Convert `CELL_PARAMETERS` into the format of `a`, and update `a` and `ATOMIC_POSITIONS` in the new `scf.in`.

---

## 3. Density of States (DOS) and Projected DOS

### Theory

The total density of states $g(E)$ counts the number of electronic states per unit energy per unit cell,

$$
g(E) = \sum_{n,\mathbf{k}} \delta(E - \varepsilon_{n\mathbf{k}})
$$

and is obtained from a non-self-consistent (`nscf`) calculation on a dense k-mesh, using the density converged in the preceding `scf` step. The projected DOS (PDOS) decomposes $g(E)$ onto atomic orbital character by projecting the Bloch states onto localized atomic basis functions $|\phi_{i}\rangle$,

$$
g_i(E) = \sum_{n,\mathbf{k}} |\langle \phi_i | \psi_{n\mathbf{k}} \rangle|^2 \, \delta(E-\varepsilon_{n\mathbf{k}})
$$

which identifies which atomic species and orbitals dominate the states near the Fermi level or band edges — essential for interpreting bonding character and (for magnetic systems) spin-resolved contributions, hence the separate spin-up/spin-down (`raw_up.dat`/`raw_down.dat`) extraction below.

### Procedure

\`\`\`bash
mpirun -np 8 pw.x < scf.in > scf.out
mpirun -np 8 pw.x < nscf.in > nscf.out
dos.x < dos.in > dos.out
\`\`\`

\`\`\`bash
awk 'NR>1 {print $1, $2}' dos.dat > raw_up.dat
awk 'NR>1 {print $1, -$3}' dos.dat > raw_down.dat
\`\`\`

\`\`\`bash
projwfc.x < pdos.in > pdos.out
\`\`\`

\`\`\`bash
sumpdos.x *\(Be\)* > atom_Be_tot.dat
sumpdos.x *\(Mg\)* > atom_Mg_tot.dat
sumpdos.x *\(Co\)* > atom_Co_tot.dat
sumpdos.x *\(Si\)* > atom_Si_tot.dat
\`\`\`

---

## 4. Band Structure

### Theory

The electronic band structure $\varepsilon_n(\mathbf{k})$ is computed along a path of high-symmetry k-points through the Brillouin zone (chosen here via XCrySDen for the F-43m lattice), starting from the converged charge density of the `scf` run. Plotting $\varepsilon_n(\mathbf{k})$ along this path reveals the fundamental (in)direct band gap, band dispersion/effective masses, and — combined with the PDOS from Section 3 — the orbital origin of the states forming the valence and conduction band edges.

### Procedure

\`\`\`bash
mpirun -np 8 pw.x < scf.in > scf.out
mpirun -np 8 pw.x < band.in > band.out
\`\`\`

> **Note:** k-points in `band.in` are generated using XCrySDen.

\`\`\`bash
bands.x < bands.in > bands.out
\`\`\`

Plot the `bands_plot.bands.gnu` file using xmgrace.

---

## 5. Phonon Dispersion and Phonon DOS

### Theory

Lattice dynamical properties follow from Density Functional Perturbation Theory (DFPT), which computes the dynamical matrix $D(\mathbf{q})$ at a set of q-points directly from the linear response of the electron density to atomic displacements. Diagonalizing $D(\mathbf{q})$,

$$
D(\mathbf{q})\, \mathbf{e}_s(\mathbf{q}) = \omega_s^2(\mathbf{q})\, \mathbf{e}_s(\mathbf{q})
$$

gives the phonon frequencies $\omega_s(\mathbf{q})$ and eigenvectors $\mathbf{e}_s(\mathbf{q})$ for each branch $s$. `q2r.x` Fourier-transforms $D(\mathbf{q})$ from the coarse q-mesh into real-space interatomic force constants, and `matdyn.x` interpolates these back onto a dense q-path (for the dispersion) or dense q-mesh (for the phonon DOS). The absence of imaginary (negative) frequencies throughout the Brillouin zone is the standard confirmation of dynamical stability of the relaxed structure. The atom-resolved phonon DOS (Be/Sc/Co/Si columns) shows which species dominate the low-frequency (acoustic) versus high-frequency (optical) branches, generally reflecting the atomic mass ordering.

### Procedure

\`\`\`bash
mpirun -np 8 pw.x < scf.in > scf.out
mpirun -np 8 ph.x < ph.in > ph.out
mpirun -np 8 q2r.x < q2r.in > q2r.out
mpirun -np 8 matdyn.x < matdyn.in > matdyn.out
plotband.x < plotband.in > plotband.out
\`\`\`

\`\`\`bash
matdyn.x < phdos.in > phdos.out
\`\`\`

\`\`\`bash
awk '{print $1,$2}' phdos.dat > total.dat
awk '{print $1,$3}' phdos.dat > Sc.dat
awk '{print $1,$4}' phdos.dat > Co.dat
awk '{print $1,$5}' phdos.dat > Si.dat
awk '{print $1,$6}' phdos.dat > Be.dat
\`\`\`

---

## 6. Thermodynamic Properties

### Theory

Within the quasi-harmonic Debye model implemented in `thermo_pw`, the vibrational free energy, entropy, and heat capacity are obtained by integrating the phonon (or Debye-approximated) density of states over the Bose-Einstein occupation factor at temperature $T$:

$$
F_{\text{vib}}(T) = k_B T \int g(\omega)\, \ln\!\left[2\sinh\!\left(\frac{\hbar\omega}{2k_BT}\right)\right] d\omega
$$

with $C_v = -T\,\partial^2 F/\partial T^2$ and $S = -\partial F/\partial T$ following directly. These, combined with the elastic constants $C_{11}, C_{12}, C_{44}$, density $\rho$, average atomic mass $M_{\text{avg}}$, and unit-cell volume $\Omega_{\text{cell}}$, feed into the Slack/Debye-based lattice thermal conductivity model used by `kl.py` to estimate the thermoelectric figure of merit $ZT$ up to a chosen maximum temperature.

### Procedure

> **Note:** Create an empty `out` folder in the working directory. You'll need `scf.in` and `thermo_control` files, then run:

\`\`\`bash
mpirun -np 8 thermo_pw.x < scf.in > scf.out
\`\`\`

> **Note:** To extract specific heat capacity, free energy, and entropy, use the following `awk` commands inside the `therm_files` folder generated after running the above:

\`\`\`bash
awk 'BEGIN{print "#T(K)   Cv(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$5*1312749.8}' output_therm.dat_debye.g1 > Cv.dat

awk 'BEGIN{print "#T(K)   Free_Energy(KJmol^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,($3*1312749.8)/1000}' output_therm.dat_debye.g1 > FreeEnergy.dat

awk 'BEGIN{print "#T(K)   Entropy(Jmol^-1K^-1)"}!/^#/{printf "%12.4f  %15.6f\n",$1,$4*1312749.8}' output_therm.dat_debye.g1 > entropy.dat
\`\`\`

> **Note:** Update `kl.py` with C11, C12, C44 from `scf.out` in GPa (dividing each by 10), $\rho$ from `vc_relax.out` (search near "final bfgs"), $M_{\text{avg}}$ (sum of all element masses from `scf.out`), $\Omega_{\text{cell}}$ from `vc_relax.out`, and finally set `Tmax` (the temperature up to which you want to calculate $ZT$).

\`\`\`bash
python3 kl.py
\`\`\`

---

## 7. Optical Properties

### Theory

Linear optical response is obtained from the frequency-dependent complex dielectric function

$$
\varepsilon(\omega) = \varepsilon_1(\omega) + i\,\varepsilon_2(\omega)
$$

computed within the independent-particle random-phase approximation by `epsilon.x`, using interband transitions between occupied and empty Kohn-Sham states obtained on a dense k-mesh (`nscf`, with `noinv = .true.` to retain the full k-point set since inversion symmetry cannot be used to reduce sampling for this response property). From $\varepsilon_1$ and $\varepsilon_2$, all other isotropic optical constants follow via standard relations: the complex refractive index

$$
n(\omega) + i\,k(\omega) = \sqrt{\varepsilon(\omega)}
$$

giving refractive index $n$, extinction coefficient $k$, reflectivity

$$
R(\omega) = \frac{(n-1)^2+k^2}{(n+1)^2+k^2}
$$

absorption coefficient $\alpha(\omega) = \dfrac{2\omega k}{c}$, optical conductivity $\sigma(\omega) \propto \omega\,\varepsilon_2(\omega)$, and the electron energy-loss function $-\mathrm{Im}[1/\varepsilon(\omega)]$, which peaks at the plasma frequency. The isotropic average $(\varepsilon_{xx}+\varepsilon_{yy}+\varepsilon_{zz})/3$ used throughout the `awk` commands is appropriate here since the cubic $F\bar43m$ structure makes the dielectric tensor isotropic.

### Procedure

> **Note:** You need norm-conserving pseudopotentials for the optical properties calculation, and must add `noinv = .true.` to the `SYSTEM` card in both `scf.in` and `nscf.in`.

\`\`\`bash
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out
epsilon.x < epsilon.in > epsilon.out
\`\`\`

\`\`\`bash
awk '{if(FNR<=2){ if(FNR==1) print "# Energy [eV]  Isotropic_Real_Dielectric_Function"; next}eps1=($2+$3+$4)/3;printf " %11.9f%11.9f\n",$1,eps1}' epsr_aiida.dat > dielectric_real_isotropic.dat

awk '{if(FNR<=2){ if(FNR==1) print "# Energy [eV]  Isotropic_Imaginary_Dielectric_Function"; next} eps2=($2+$3+$4)/3; printf " %11.9f%11.9f\n",$1,eps2}' epsi_aiida.dat > dielectric_imaginary_isotropic.dat

awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Reflectivity [fraction]"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2);  n=sqrt((mod+r[FNR])/2); k=sqrt((mod-r[FNR])/2); R=((n-1)^2 + k^2)/((n+1)^2 + k^2); printf "    %11.9f    %11.9f\n", $1, R }' epsr_aiida.dat epsi_aiida.dat > reflectivity_isotropic.dat

awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Refractive_Index_n"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2); n=sqrt((mod+r[FNR])/2); printf "    %11.9f    %11.9f\n", $1, n }' epsr_aiida.dat epsi_aiida.dat > refractive_index_isotropic.dat

awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Extinction_Coefficient_k"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2); k=sqrt((mod-r[FNR])/2); printf "    %11.9f    %11.9f\n", $1, k }' epsr_aiida.dat epsi_aiida.dat > extinction_coefficient_isotropic.dat

awk 'NR==FNR { if(FNR>2) r[FNR]=($2+$3+$4)/3; next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Absorption_Coefficient [10^4/cm]"; next } i_avg=($2+$3+$4)/3; mod=sqrt(r[FNR]^2 + i_avg^2); k=sqrt((mod-r[FNR])/2); alpha=(2*$1*k*1.6231012e5)/10000; printf "    %11.9f    %14.6f\n", $1, alpha }' epsr_aiida.dat epsi_aiida.dat > absorption_isotropic_scaled.dat

awk '{ if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_Optical_Conductivity [10^3 Omega^-1 cm^-1]"; next } i_avg=($2+$3+$4)/3; sigma=(1327.21*$1*i_avg)/1000; printf "    %11.9f    %14.6f\n", $1, sigma }' epsi_aiida.dat > optical_conductivity_isotropic_scaled.dat

awk 'NR==FNR { if(FNR>2) { r_avg=($2+$3+$4)/3; r[FNR]=r_avg } next } { if(FNR<=2) { if(FNR==1) print "# Energy [eV]  Isotropic_EELS"; next } i_avg=($2+$3+$4)/3; loss = i_avg / (r[FNR]^2 + i_avg^2); printf "    %11.9f    %11.9f\n", $1, loss }' epsr_aiida.dat epsi_aiida.dat > energylossfunction_isotropic.dat
\`\`\`

---

## 8. Raman Spectroscopy

### Theory

Raman-active zone-center ($\Gamma$-point, $\mathbf{q}=0$) phonon modes are computed via DFPT in `ph.x`, and their Raman tensor / activity is extracted by `dynmat.x`. Because $F\bar43m$ ($T_d$) is a non-centrosymmetric group, its zone-center optical phonons can be simultaneously Raman- and IR-active. Each identified Raman-active mode has a characteristic frequency $\omega_0$ and Raman intensity $I_0$; the simulated spectrum is constructed as a sum of Gaussian peaks,

$$
I(\omega) = \sum_i I_{0,i} \, \exp\!\left[-\left(\frac{\omega-\omega_{0,i}}{2}\right)^2\right]
$$

which approximates the finite linewidth seen in experimental Raman spectra (arising from anharmonic phonon lifetimes not captured at the harmonic DFPT level) and allows direct visual/qualitative comparison to experimental Raman data.

### Procedure

\`\`\`bash
pw.x < scf.in > scf.out
ph.x < ph_raman.in > ph_raman.out
\`\`\`

> **Note:** You need `pz-hgh` pseudopotentials.

\`\`\`bash
python3 -c "import numpy as np; peak=[(183.62,13.0937), (280.46,251.2447),(401.29,308.4624)]; w=np.linspace(100, 500, 800); fit=sum(I0*np.exp(-((w-w0)/2)**2) for w0, I0 in peak); np.savetxt('raman_curve.dat', np.column_stack((w, fit)), fmt='%.4f')"
\`\`\`

> **Note:** From `dynmat.out`, look for the double frequencies and replace `(183.62,13.0937), (280.46,251.2447), (401.29,308.4624)` above — 183.62, 280.46, and 401.29 are the frequencies, and 13.0937, 251.2447, 308.4624 are the corresponding Raman values.

---

## 9. Piezoelectricity

### Theory: generating the strained inputs (`piezo.py`)

For cubic $T_d$ ($\overline{4}3m$) symmetry, the piezoelectric tensor has a single independent component:

$$
d_{14} = d_{25} = d_{36}
$$

A pure shear strain $\eta_4$ (coupling $y$–$z$) induces polarization exactly along Cartesian $x$ — this symmetry argument is why $d_{14}$ is the only independent constant for this point group. `piezo.py` applies a series of small $\eta_4$ shear strains to the relaxed FCC cell and, for each strain point, writes a two-step calculation: an ordinary `scf` run to converge the charge density, followed by an `nscf` run with `lberry = .true.` (Berry-phase polarization, King-Smith–Vanderbilt formalism) that reads that density. `nosym = .true.` is required on both steps because strain lowers the cell's symmetry below the parent cubic group, and allowing QE to reduce the k-point set via symmetry corrupts the Berry-phase string. Atomic fractional coordinates are held fixed under strain (clamped-ion approximation), so the calculation yields only the electronic piezoelectric response, omitting the internal-strain (relaxed-ion) contribution — a standard, defensible approximation that should be stated explicitly in the thesis methods text.

### Theory: extracting e₁₄ from Berry-phase output (`extract_e14.py`)

**Background.** For cubic $T_d$ ($\overline{4}3m$) symmetry, the piezoelectric tensor has a single independent component:

$$
d_{14} = d_{25} = d_{36}
$$

A pure shear strain $\eta_4$ (coupling $y$–$z$) induces a polarization that points exactly along Cartesian $x$, with $P_y = P_z = 0$ identically. This is the textbook symmetry argument that justifies why $d_{14}$ is the only independent constant for this point group. The stress constant $e_{14}$ is defined as the linear response:

$$
e_{14} = \left.\frac{dP_x}{d\eta_4}\right|_{\eta_4 = 0}
$$

and the strain constant follows from the elastic constant $C_{44}$:

$$
d_{14} = \frac{e_{14}}{C_{44}}
$$

**The non-Cartesian cell problem.** The DFT cell uses `ibrav = 2` primitive FCC lattice vectors, which are not aligned with Cartesian $x/y/z$ — the first primitive vector points along

$$
\hat A_1 = \frac{1}{\sqrt2}(-1, 0, 1)
$$

not along $\hat x$. Running Berry-phase with `gdir = 1` therefore does not return $P_x$ directly — QE returns the polarization projected onto $\hat A_1$:

$$
P_{\text{measured}} = \vec P \cdot \hat A_1 = P_x \, d_x
$$

(using $P_y = P_z = 0$ from the symmetry argument above, so only the $P_x d_x$ term survives). QE prints $\hat A_1 = (d_x, d_y, d_z)$ explicitly as "The polarization direction is: ...". The true $P_x$ is recovered by:

$$
P_x = \frac{P_{\text{measured}}}{d_x}
$$

computed per strain point, not with a fixed constant — straining the cell slightly rotates $\hat A_1$, so $d_x$ drifts point to point.

**What the script does:**

1. Parses each `_nscf.out` — regex-extracts `P = ... (mod ...) C/m^2` and the printed direction vector $(d_x, d_y, d_z)$.
2. Applies the projection correction above to recover $P_x$ for each strain point, using that file's own $d_x$.
3. Sanity-checks each point: skips files with no polarisation block (failed/crashed run), skips points where $d_x \approx 0$ (correction would diverge), and flags points where $(d_x,d_y,d_z)$ deviates noticeably from the expected $(-0.707, 0, 0.707)$-type pattern.
4. Fits a linear least-squares (no numpy/scipy — closed-form slope/intercept) of

$$
P_x = e_{14}\,\eta_4 + P_0
$$

   reporting $R^2$; warns if $R^2 < 0.99$, since a poor fit usually signals a branch jump — Berry-phase polarisation is only defined modulo a quantum $eR/\Omega$, so if QE's chosen branch jumps between adjacent strain points, the data won't lie on a clean line.
5. Converts to $d_{14}$ using the user-supplied $C_{44}$ (GPa), with SI-to-pm/V unit conversion:

$$
d_{14}\ [\text{pm/V}] = \frac{e_{14}}{C_{44}} \times 10^{12}
$$

### Procedure

\`\`\`bash
python3 piezo.py
\`\`\`

> **Note:** For a different quaternary Heusler alloy, update in `piezo.py`:
> - `COMPOUND_LABEL`, `ALAT_ANGSTROM` — new compound name and relaxed lattice constant
> - `ATOMIC_SPECIES` — masses and pseudopotential filenames for the new elements
> - `ATOMIC_POSITIONS` — fractional coordinates for the DFT-confirmed site ordering (Type I/II/III) of the new compound; do not reuse the old dict as-is
> - `NBND` — recompute from the new total valence electron count
> - `ECUTWFC`, `ECUTRHO` — re-converge for the new pseudopotentials
> - If the new compound is magnetic: add `nspin = 2` and `starting_magnetization(ityp)` to both `scf` and `nscf` blocks
> - `GDIR = 1`, `KPTS_BASE`, `NPPSTR`, `STRAIN_VALUES`, and the yz-shear strain function stay unchanged as long as the compound is $F\bar43m$ ($T_d$), since $d_{14}=d_{25}=d_{36}$ still holds

> **Note:** After running `python3 piezo.py`, it will produce a `piezo_inputs` folder — copy the pseudopotentials into it, then run:

\`\`\`bash
pw.x < BeScCoSi_eta00_scf.in > BeScCoSi_eta00_scf.out
pw.x < BeScCoSi_eta00_nscf.in > BeScCoSi_eta00_nscf.out
\`\`\`

> **Note:** Run all `scf.in`/`nscf.in` pairs, then move out of `piezo_inputs`, keep `extract_e14.py` alongside `piezo.py`, and take `C44_GPa` from `thermo`'s `scf.out`.

**Usage:**

\`\`\`bash
python3 extract_e14.py <compound_label> <C44_GPa> [piezo_inputs_dir]
\`\`\`

**Example:**

\`\`\`bash
python3 extract_e14.py BeScCoSi 67.62 ./piezo_inputs
\`\`\`

**Summary:** `piezo.py` generates the strained structures and runs the `scf → nscf(lberry)` pipeline; `extract_e14.py` is the analysis half — it reads back the Berry-phase results, undoes the coordinate-system artefact from the non-Cartesian cell, fits the slope, and reports the final $d_{14}$.

---

## 10. Formation Energy

### Theory

The formation energy quantifies the thermodynamic stability of the compound relative to its constituent elements in their reference (bulk, elemental) phases. It is defined as:

$$
E_f = E_{\text{compound}} - \sum_i n_i\, E_i^{\text{atom}}
$$

where $E_{\text{compound}}$ is the total DFT energy of the relaxed compound per formula unit, $E_i^{\text{atom}} = E_i / N_i$ is the per-atom reference energy of elemental species $i$ (computed separately, from that element's own bulk/ground-state structure, dividing its total energy by its number of atoms in that calculation), and $n_i$ is the number of atoms of species $i$ per formula unit of the compound. A negative $E_f$ indicates the compound is energetically favourable to form from its constituent elements (thermodynamically stable against decomposition into the elemental phases); a positive value indicates instability.

Both $E_{\text{compound}}$ and each $E_i$ are read directly from Quantum ESPRESSO's `scf.out`-type files: the script locates the number of atoms per cell (`"number of atoms/cell"`) and the final converged total energy (the line beginning with `!`, QE's marker for the self-consistent total energy in Ry). The elemental per-atom energies are summed and subtracted from the compound energy to give $E_f$ per formula unit, which is then converted from Ry to eV (1 Ry = 13.605693 eV) and reported both per formula unit and per atom:

$$
E_f\,[\text{eV/f.u.}] = E_f\,[\text{Ry/f.u.}] \times 13.605693, \qquad
E_f\,[\text{eV/atom}] = \frac{E_f\,[\text{eV/f.u.}]}{N_{\text{atoms in compound}}}
$$

> **Caveat worth stating in the thesis methods text:** this definition uses the elements' DFT total energies directly as computed (whatever structure/settings were used for `Mg.out`, `Sc.out`, `Co.out`, `Si.out`), so the reference states must be each element's correct experimental/most-stable bulk allotrope (e.g. hcp-Mg, not an arbitrary or unrelaxed cell) — an inconsistent reference state changes $E_f$ without reflecting a real change in compound stability.

### Procedure

Run separate `scf` calculations for the compound and for each elemental reference:

\`\`\`bash
pw.x < scf.in > scf.out      # compound
pw.x < Be.in > Be.out
pw.x < Sc.in > Sc.out
pw.x < Co.in > Co.out
pw.x < Si.in > Si.out
\`\`\`

Place all output files in the same working directory as the script, then run:

\`\`\`bash
python3 formation_energy.py
\`\`\`

> **Note:** For a different quaternary Heusler alloy, update in the script:
> - `compound_file` — path to the new compound's `scf.out`
> - `elements` dictionary — element labels and their corresponding elemental reference `.out` files, matching the new compound's constituent species
> - Everything else (energy/`nat` parsing, the $E_f$ formula, and unit conversion) stays unchanged

**Output:** Prints each element's per-atom reference energy, followed by the compound energy, summed elemental reference energy, and formation energy in Ry/f.u., eV/f.u., and eV/atom.
