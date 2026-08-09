import numpy as np

# ============================================================
# INPUT PARAMETERS
# ============================================================

# Elastic constants (GPa)
C11 = 198.377761
C12 = 111.980811
C44 = 67.619697

# Density (g/cm^3 -> kg/m^3)
rho =  4.82421* 1000.0

# Average atomic mass (amu)
Mavg = 35.246695
# Atoms per primitive cell
nat = 4

# Primitive-cell volume (bohr^3)
Omega_cell =  331.1373
# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

h = 6.62607015e-34          # J s
kB = 1.380649e-23           # J/K
amu = 1.66053906660e-27     # kg
bohr = 5.29177210903e-11    # m

# ============================================================
# VOIGT-REUSS-HILL ELASTIC MODULI
# ============================================================

B = (C11 + 2*C12)/3.0

GV = (C11 - C12 + 3*C44)/5.0
GR = (5*C44*(C11 - C12))/(4*C44 + 3*(C11 - C12))
G = (GV + GR)/2.0

B_pa = B * 1e9
G_pa = G * 1e9

# ============================================================
# SOUND VELOCITIES
# ============================================================

vt = np.sqrt(G_pa/rho)

vl = np.sqrt((B_pa + 4.0*G_pa/3.0)/rho)

vm = ((1.0/3.0) * (2.0/vt**3 + 1.0/vl**3))**(-1.0/3.0)

# ============================================================
# POISSON'S RATIO
# ============================================================

nu = (3*B - 2*G)/(2*(3*B + G))

# ============================================================
# GRUNEISEN PARAMETER
# ============================================================

gamma = 3.0*(1.0 + nu)/(2.0*(2.0 - 3.0*nu))

if gamma <= 0:
    raise ValueError("Non-physical Gruneisen parameter.")

# ============================================================
# DEBYE TEMPERATURE
# ============================================================

volume_m3 = Omega_cell * (bohr**3)

theta_D = (h/kB) * ((3*nat)/(4*np.pi*volume_m3))**(1/3) * vm

# ============================================================
# ATOMIC VOLUME
# ============================================================

Vatom_ang3 = Omega_cell * (0.529177210903**3) / nat

# ============================================================
# SLACK MODEL
# ============================================================

A = 2.43e-6 / (1.0 - 0.514/gamma + 0.228/gamma**2)

kappa_const = (
    A
    * Mavg
    * theta_D**3
    * Vatom_ang3**(1/3)
    / (gamma**2 * nat**(2/3))
)

# ============================================================
# PRINT RESULTS
# ============================================================

print(f"Bulk modulus B          = {B:.3f} GPa")
print(f"Shear modulus G         = {G:.3f} GPa")
print(f"Poisson ratio           = {nu:.4f}")
print(f"Gruneisen parameter     = {gamma:.4f}")
print(f"Transverse velocity     = {vt:.2f} m/s")
print(f"Longitudinal velocity   = {vl:.2f} m/s")
print(f"Average velocity        = {vm:.2f} m/s")
print(f"Debye temperature       = {theta_D:.2f} K")
print(f"Slack prefactor         = {kappa_const:.4f}")

# ============================================================
# GENERATE kl.dat
# ============================================================

Tmin = 300
Tmax = 1301
step = 100

with open("kl.dat", "w") as f:
    f.write("# Temperature(K)   Lattice_Thermal_Conductivity(W/mK)\n")

    for T in range(Tmin, Tmax + step, step):
        kappa = kappa_const / T
        f.write(f"{T:6d} {kappa:15.6f}\n")

print("\nkl.dat successfully written.")
