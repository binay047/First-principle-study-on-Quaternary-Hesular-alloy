"""
Piezoelectric (d14) input generator for cubic T_d Heusler XScCoSi
--------------------------------------------------------------------
Generates a TWO-STEP series of pw.x input files (scf + nscf) with small
yz shear strains (eta_4) applied to the relaxed FCC cell, for a
Berry-phase polarization calculation (lberry = .true., gdir = 1).

WHY TWO STEPS:
  QE's Berry-phase (lberry) calculation needs a converged charge density
  first. Running lberry directly inside a single 'scf' calculation is
  unreliable -- the standard, working recipe (confirmed against QE's
  own example10 and multiple verified user reports) is:
    Step 1: ordinary 'scf' run (no lberry) -> converges the density.
    Step 2: 'nscf' run with lberry=.true., reading that density,
            with the k-mesh component along `gdir` set EQUAL to
            `nppstr` (this is required -- a mismatch here is why a
            previous version of this script produced no polarization
            output at all).

  NOTE (patched): nosym = .true. has been added to BOTH steps. Berry-
  phase calculations are unreliable if QE is allowed to use symmetry
  to reduce the k-point set -- the strained cell has lower symmetry
  than the parent cubic cell, and leaving symmetry reduction on is the
  most common reason this type of calculation silently produces no
  usable polarization output.

Physical background:
  For T_d (-43m) symmetry, the only independent piezoelectric strain
  constant is d14 (= d25 = d36 by symmetry). Applying shear strain
  eta_4 (coupling y-z) induces polarization along x, so gdir = 1.

  P_x(eta_4) is linear near eta_4 = 0:
      e14 = dP_x / d(eta_4)  at eta_4 = 0   [Berry-phase piezoelectric
                                              stress constant]
      d14 = e14 / C44                       [piezoelectric strain
                                              constant, using your
                                              already-computed C44]

  CAVEAT: atomic fractional coordinates are held FIXED under strain
  (clamped-ion approximation). This yields the clamped-ion (electronic)
  piezoelectric response only, and omits the internal-strain
  (relaxed-ion) contribution from atoms shifting within the unit cell
  under strain. This should be stated explicitly in the thesis methods
  text -- it is a standard, defensible approximation, not an error.

Workflow:
  1. Edit the INPUT PARAMETERS block below to match your compound.
  2. Run this script -- it writes, per strain value, TWO files into
     ./piezo_inputs/:  <label>_etaNN_scf.in  and  <label>_etaNN_nscf.in
  3. For each strain index NN, run scf FIRST, then nscf (same
     outdir/prefix, so nscf reads the scf density automatically):
       pw.x < X_eta00_scf.in  > X_eta00_scf.out
       pw.x < X_eta00_nscf.in > X_eta00_nscf.out
  4. In the _nscf.out file, look for the Berry-phase polarization
     block (search for "BERRY PHASE", "Polarization", or "P =" --
     exact heading varies slightly by QE version).
  5. Fit P_x vs eta_4 linearly; the slope at eta_4=0 is e14.
  6. d14 = e14 / C44  (use your already-computed C44).
"""

import numpy as np
import os

# ----------------------------------------------------------------------
# 1. INPUT PARAMETERS -- edit to match your compound / scf.in
# ----------------------------------------------------------------------

COMPOUND_LABEL = "BeScCoSi"          # used only for filenames/comments

ALAT_ANGSTROM = 5.81153              # equilibrium lattice constant (A)
                                      # (same as your ibrav=2 "a" value)

ECUTWFC = 55
ECUTRHO = 550
CONV_THR = 1.0e-7
MIXING_BETA = 0.4
ELECTRON_MAXSTEP = 800
NBND = 32

PREFIX = "aiida"
OUTDIR = "./out/"
PSEUDO_DIR = "./pseudo/"

ATOMIC_SPECIES = [
    ("Be", 9.012182,  "Be.pbe-n-rrkjus_psl.1.0.0.UPF"),
    ("Sc", 44.955912, "Sc.pbe-spn-rrkjus_psl.1.0.0.UPF"),
    ("Co", 58.933195, "Co.pbe-spn-rrkjus_psl.0.3.1.UPF"),
    ("Si", 28.0855,   "Si.pbe-nl-rrkjus_psl.1.0.0.UPF"),
]

# Fractional (crystal) coordinates -- phase III, unchanged under
# homogeneous (clamped-ion) strain.
ATOMIC_POSITIONS = {
    "Be": np.array([0.75, 0.75, 0.75]),
    "Sc": np.array([0.00, 0.00, 0.00]),
    "Co": np.array([0.25, 0.25, 0.25]),
    "Si": np.array([0.50, 0.50, 0.50]),
}

# Shear strain values (eta_4, dimensionless engineering strain,
# coupling y-z). Keep these small to stay in the linear regime --
# check linearity of P_x(eta_4) before trusting the fitted slope.
STRAIN_VALUES = [-0.02, -0.015, -0.01, -0.005, 0.0, 0.005, 0.01, 0.015, 0.02]

# k-point mesh for the ORDINARY scf run (step 1), and for the two
# directions PERPENDICULAR to gdir in the nscf/Berry-phase run.
KPTS_BASE = (7, 7, 7)

# Number of k-points along the Berry-phase string (gdir direction).
# In the nscf run, this value REPLACES the base mesh entry in the
# gdir slot (this matching is required -- see header note above).
NPPSTR = 8

# Berry-phase direction: gdir = 1 (x) is correct for eta_4 (yz shear)
# per the T_d symmetry argument above. Do not change unless you are
# computing a different tensor component.
GDIR = 1

OUTPUT_DIR = "./piezo_inputs"

# ----------------------------------------------------------------------
# 2. BUILD FCC PRIMITIVE VECTORS (ibrav = 2 convention) AND APPLY STRAIN
# ----------------------------------------------------------------------

def fcc_primitive_vectors(alat):
    """
    Standard QE ibrav=2 FCC primitive vectors (as columns):
        v1 = a/2 * (-1, 0, 1)
        v2 = a/2 * ( 0, 1, 1)
        v3 = a/2 * (-1, 1, 0)
    """
    a = alat
    v1 = np.array([-a/2, 0.0,  a/2])
    v2 = np.array([ 0.0, a/2,  a/2])
    v3 = np.array([-a/2, a/2,  0.0])
    return np.column_stack([v1, v2, v3])


def apply_yz_shear(V, eta4):
    """
    Applies a symmetric engineering shear strain eta_4 (coupling y-z):
        eps = [[0,      0,       0     ],
               [0,      0,       eta/2 ],
               [0,      eta/2,   0     ]]
        v' = (I + eps) . v
    """
    eps = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, eta4 / 2.0],
        [0.0, eta4 / 2.0, 0.0],
    ])
    return (np.eye(3) + eps) @ V


# ----------------------------------------------------------------------
# 3. SHARED BLOCK BUILDERS
# ----------------------------------------------------------------------

def build_common_blocks(V_strained):
    species_block = "\n".join(
        f"{sym}  {mass:.6f}  {upf}" for sym, mass, upf in ATOMIC_SPECIES
    )
    positions_block = "\n".join(
        f"{sym}  {pos[0]:.10f}  {pos[1]:.10f}  {pos[2]:.10f}"
        for sym, pos in ATOMIC_POSITIONS.items()
    )
    cell_block = "\n".join(
        f"  {V_strained[0, j]:.10f}  {V_strained[1, j]:.10f}  {V_strained[2, j]:.10f}"
        for j in range(3)
    )
    return species_block, positions_block, cell_block


def gdir_mesh(base_mesh, gdir, nppstr):
    """
    Builds the K_POINTS automatic mesh triple with nppstr placed in the
    slot corresponding to gdir (1=x, 2=y, 3=z), matching the pattern
    used in verified QE Berry-phase examples (e.g. mesh = 8 8 24 with
    gdir=3, nppstr=24).
    """
    mesh = list(base_mesh)
    mesh[gdir - 1] = nppstr
    return tuple(mesh)


# ----------------------------------------------------------------------
# 4. WRITE STEP 1 (scf) AND STEP 2 (nscf + lberry) INPUT FILES
# ----------------------------------------------------------------------

def write_scf_input(V_strained, index):
    """Step 1: ordinary SCF run (no lberry) to converge the density for
    this strained structure."""
    fname = os.path.join(
        OUTPUT_DIR, f"{COMPOUND_LABEL}_eta{index:02d}_scf.in"
    )
    species_block, positions_block, cell_block = build_common_blocks(V_strained)

    text = f"""&CONTROL
  calculation = 'scf'
  outdir = '{OUTDIR}'
  prefix = '{PREFIX}_eta{index:02d}'
  pseudo_dir = '{PSEUDO_DIR}'
  verbosity = 'high'
/
&SYSTEM
  ibrav = 0
  nat = 4
  ntyp = 4
  ecutwfc = {ECUTWFC}
  ecutrho = {ECUTRHO}
  nbnd = {NBND}
  occupations = 'fixed'
  nosym = .true.
/
&ELECTRONS
  conv_thr = {CONV_THR}
  electron_maxstep = {ELECTRON_MAXSTEP}
  mixing_beta = {MIXING_BETA}
/
CELL_PARAMETERS angstrom
{cell_block}
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS crystal
{positions_block}
K_POINTS automatic
{KPTS_BASE[0]} {KPTS_BASE[1]} {KPTS_BASE[2]} 0 0 0
"""
    with open(fname, "w") as f:
        f.write(text)
    return fname


def write_nscf_berry_input(V_strained, index):
    """Step 2: nscf run with lberry=.true., reading the density from the
    step-1 scf run (same outdir/prefix). The k-mesh slot along gdir is
    set equal to nppstr -- required for the Berry-phase strings to be
    built correctly."""
    fname = os.path.join(
        OUTPUT_DIR, f"{COMPOUND_LABEL}_eta{index:02d}_nscf.in"
    )
    species_block, positions_block, cell_block = build_common_blocks(V_strained)
    mesh = gdir_mesh(KPTS_BASE, GDIR, NPPSTR)

    text = f"""&CONTROL
  calculation = 'nscf'
  outdir = '{OUTDIR}'
  prefix = '{PREFIX}_eta{index:02d}'
  pseudo_dir = '{PSEUDO_DIR}'
  verbosity = 'high'
  lberry = .true.
  gdir = {GDIR}
  nppstr = {NPPSTR}
/
&SYSTEM
  ibrav = 0
  nat = 4
  ntyp = 4
  ecutwfc = {ECUTWFC}
  ecutrho = {ECUTRHO}
  nbnd = {NBND}
  occupations = 'fixed'
  nosym = .true.
/
&ELECTRONS
  conv_thr = {CONV_THR}
  electron_maxstep = {ELECTRON_MAXSTEP}
  mixing_beta = {MIXING_BETA}
/
CELL_PARAMETERS angstrom
{cell_block}
ATOMIC_SPECIES
{species_block}
ATOMIC_POSITIONS crystal
{positions_block}
K_POINTS automatic
{mesh[0]} {mesh[1]} {mesh[2]} 0 0 0
"""
    with open(fname, "w") as f:
        f.write(text)
    return fname


# ----------------------------------------------------------------------
# 5. MAIN
# ----------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    V0 = fcc_primitive_vectors(ALAT_ANGSTROM)

    print(f"Base FCC primitive vectors (Angstrom) for {COMPOUND_LABEL}:")
    print(np.round(V0, 6))
    print()

    mesh = gdir_mesh(KPTS_BASE, GDIR, NPPSTR)
    print(f"nscf K_POINTS mesh (gdir={GDIR} slot set to nppstr={NPPSTR}): "
          f"{mesh[0]} {mesh[1]} {mesh[2]}")
    print()

    written = []
    for i, eta4 in enumerate(STRAIN_VALUES):
        V_strained = apply_yz_shear(V0, eta4)
        scf_f = write_scf_input(V_strained, i)
        nscf_f = write_nscf_berry_input(V_strained, i)
        written.append((scf_f, nscf_f, eta4))
        print(f"  eta_4 = {eta4:+.4f}:")
        print(f"    {scf_f}")
        print(f"    {nscf_f}")

    print()
    print("Next steps:")
    print("  1. Run EACH PAIR IN ORDER (scf before nscf, same index):")
    print(f"       for i in 00 01 02 03 04 05 06 07 08; do")
    print(f"         pw.x < {OUTPUT_DIR}/{COMPOUND_LABEL}_eta${{i}}_scf.in  > {OUTPUT_DIR}/{COMPOUND_LABEL}_eta${{i}}_scf.out")
    print(f"         pw.x < {OUTPUT_DIR}/{COMPOUND_LABEL}_eta${{i}}_nscf.in > {OUTPUT_DIR}/{COMPOUND_LABEL}_eta${{i}}_nscf.out")
    print(f"       done")
    print("  2. Search each *_nscf.out for the Berry-phase polarization")
    print("     block (heading text varies by QE version -- try:")
    print(f"       grep -i -A5 'berry phase\\|polarization\\|P =' {OUTPUT_DIR}/*_nscf.out")
    print("  3. Fit P_x vs eta_4 linearly; the slope at eta_4=0 is e14.")
    print("  4. d14 = e14 / C44  (use your already-computed C44).")
    print()
    print("Strain series (for your reference when fitting):")
    for scf_f, nscf_f, eta4 in written:
        print(f"    eta_4 = {eta4:+.4f}  ->  {os.path.basename(nscf_f)}")


if __name__ == "__main__":
    main()
