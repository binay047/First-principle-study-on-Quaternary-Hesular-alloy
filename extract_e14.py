"""
Extract P_x vs eta_4 from a series of QE Berry-phase (lberry) nscf
outputs and fit the piezoelectric stress constant e14.

IMPORTANT (read this before trusting results):
  Because the DFT cell uses primitive FCC lattice vectors (ibrav=2),
  which are NOT aligned with Cartesian x/y/z, a gdir=1 Berry-phase run
  reports the polarization projected onto A1 = (-1,0,1)/sqrt(2), not
  the Cartesian P_x directly. QE prints this direction explicitly as
  "The polarization direction is: (dx, dy, dz)" in the output.

  By T_d (-43m) symmetry, a PURE eta_4 (yz) shear strain induces a
  physical polarization vector that points exactly along Cartesian x
  (P_y = P_z = 0 identically -- this is the textbook justification for
  why e14 is the only independent piezoelectric constant for this
  point group). Given that, the measured scalar is simply

      P_measured = P . A1_hat = Px * dx

  so we recover the true Px by dividing by dx (the x-component of the
  printed direction vector) for EACH strain point individually (strain
  slightly rotates A1 away from its unstrained direction, so dx is not
  exactly constant across the series -- always use the value QE prints
  for that specific run, not a hardcoded constant).

Usage:
    python3 extract_e14.py <compound_label> <C44_GPa> [piezo_inputs_dir]

Example:
    python3 extract_e14.py BeScCoSi 67.62 ./piezo_inputs
"""

import sys
import re
import glob
import os

STRAIN_VALUES = [-0.02, -0.015, -0.01, -0.005, 0.0, 0.005, 0.01, 0.015, 0.02]


def parse_polarization_and_direction(filepath):
    """
    Returns (P_measured_C_per_m2, direction_vector) or (None, None) if
    not found. Parses QE's Berry-phase output block:

        P =  -1.0128908  (mod   1.3419264)  C/m^2
        The polarization direction is:  (-0.70709 ,-0.00707 , 0.70709 )
    """
    with open(filepath, "r") as f:
        text = f.read()

    p_matches = re.findall(
        r"P\s*=\s*([-\d.Ee+]+)\s*\(mod\s*[-\d.Ee+]+\)\s*C/m\^2",
        text,
    )
    if not p_matches:
        return None, None
    p_measured = float(p_matches[-1])

    dir_match = re.search(
        r"polarization direction is:\s*\(\s*([-\d.Ee+]+)\s*,\s*"
        r"([-\d.Ee+]+)\s*,\s*([-\d.Ee+]+)\s*\)",
        text,
    )
    if not dir_match:
        return p_measured, None

    direction = (
        float(dir_match.group(1)),
        float(dir_match.group(2)),
        float(dir_match.group(3)),
    )
    return p_measured, direction


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    label = sys.argv[1]
    c44_gpa = float(sys.argv[2])
    piezo_dir = sys.argv[3] if len(sys.argv) > 3 else "./piezo_inputs"

    files = sorted(glob.glob(os.path.join(piezo_dir, f"{label}_eta*_nscf.out")))
    if not files:
        print(f"No files matching {label}_eta*_nscf.out found in {piezo_dir}")
        sys.exit(1)

    data = []
    for fpath in files:
        m = re.search(r"_eta(\d+)_nscf\.out$", fpath)
        if not m:
            continue
        idx = int(m.group(1))
        if idx >= len(STRAIN_VALUES):
            print(f"WARNING: index {idx} from {fpath} has no matching "
                  f"strain value -- skipping.")
            continue
        eta4 = STRAIN_VALUES[idx]

        p_measured, direction = parse_polarization_and_direction(fpath)
        if p_measured is None:
            print(f"WARNING: could not find polarization block in {fpath} "
                  f"-- run may have failed or crashed. Skipping.")
            continue
        if direction is None:
            print(f"WARNING: found P but not the direction vector in "
                  f"{fpath} -- cannot apply projection correction. "
                  f"Skipping this point.")
            continue

        dx = direction[0]
        if abs(dx) < 1e-4:
            print(f"WARNING: direction x-component in {fpath} is "
                  f"~0 ({dx:.5f}) -- projection correction would blow up. "
                  f"Skipping this point (check this run manually).")
            continue

        p_x = p_measured / dx
        data.append((eta4, p_x, p_measured, direction, fpath))

    if len(data) < 2:
        print("Not enough successfully-parsed points to fit a line.")
        sys.exit(1)

    data.sort()
    print(f"{'eta_4':>10} {'P_x (C/m^2)':>14} {'P_meas (C/m^2)':>16} "
          f"{'direction (dx,dy,dz)':>30}   file")
    for eta4, p_x, p_measured, direction, fpath in data:
        dstr = f"({direction[0]:+.4f},{direction[1]:+.4f},{direction[2]:+.4f})"
        print(f"{eta4:>10.4f} {p_x:>14.6f} {p_measured:>16.6f} "
              f"{dstr:>30}   {os.path.basename(fpath)}")
        if abs(direction[1]) > 0.05 or abs(direction[2] - abs(direction[0])) > 0.05:
            print(f"    NOTE: direction vector for this point deviates "
                  f"noticeably from the expected (-0.707,0,0.707)-type "
                  f"pattern -- worth a manual look at {fpath}.")

    n = len(data)
    xs = [d[0] for d in data]
    ys = [d[1] for d in data]
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    sxx = sum((x - xbar) ** 2 for x in xs)
    e14 = sxy / sxx
    p0 = ybar - e14 * xbar

    ss_res = sum((y - (e14 * x + p0)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    print()
    print(f"Linear fit: P_x = e14 * eta_4 + P0")
    print(f"  e14 = {e14:.6e}  C/m^2")
    print(f"  P0  = {p0:.6e}  C/m^2  (should be ~0 at eta_4=0)")
    print(f"  R^2 = {r2:.6f}")
    if r2 < 0.99:
        print("  WARNING: R^2 < 0.99 -- check for a jump/discontinuity "
              "between adjacent points (Berry-phase branch/quantum "
              "ambiguity is the most likely cause), or a bad SCF run.")

    d14_pm_per_V = (e14 / (c44_gpa * 1e9)) * 1e12
    print()
    print(f"Using C44 = {c44_gpa} GPa:")
    print(f"  d14 = e14 / C44 = {d14_pm_per_V:.4f}  pm/V")


if __name__ == "__main__":
    main()
