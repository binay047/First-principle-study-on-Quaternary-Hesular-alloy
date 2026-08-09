#!/usr/bin/env python3

import re

# Files
compound_file = "scf.out"

elements = {
    "Be": "Be.out",
    "Sc": "Sc.out",
    "Co": "Co.out",
    "Si": "Si.out"
}


def get_energy_and_nat(filename):
    energy = None
    nat = None

    with open(filename, "r") as f:
        for line in f:
            if "number of atoms/cell" in line:
                nat = int(line.split("=")[1])

            if line.strip().startswith("!"):
                m = re.search(r'=\s*([-0-9.]+)\s*Ry', line)
                if m:
                    energy = float(m.group(1))

    if energy is None:
        raise ValueError(f"Could not find energy in {filename}")

    if nat is None:
        raise ValueError(f"Could not find nat in {filename}")

    return energy, nat


# Compound energy
E_compound, nat_compound = get_energy_and_nat(compound_file)

# Sum elemental energies per atom
E_elements = 0.0

print("Elemental reference energies:")
for elem, file in elements.items():
    E, nat = get_energy_and_nat(file)
    E_atom = E / nat
    E_elements += E_atom

    print(f"{elem:2s}: {E_atom:12.6f} Ry/atom")

# Formation energy
Ef_Ry = E_compound - E_elements
Ef_eV_fu = Ef_Ry * 13.605693
Ef_eV_atom = Ef_eV_fu / nat_compound

print("\nResults")
print("-" * 40)
print(f"Compound energy      : {E_compound:.6f} Ry")
print(f"Element sum          : {E_elements:.6f} Ry")
print(f"Formation energy     : {Ef_Ry:.6f} Ry/f.u.")
print(f"Formation energy     : {Ef_eV_fu:.6f} eV/f.u.")
print(f"Formation energy     : {Ef_eV_atom:.6f} eV/atom")
