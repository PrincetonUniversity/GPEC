#!/usr/bin/env python3
import numpy as np

# Load files
Q_data = np.loadtxt("gpec_recon_Qvec_sol1.out")
Ccov_data = np.loadtxt("gpec_recon_Ccov_sol1.out")
Cphy_data = np.loadtxt("gpec_recon_Ccontra_physics_sol1.out")
Cmet_data = np.loadtxt("gpec_recon_Ccontra_metric_sol1.out")

print("=" * 100)
print("Q VECTOR vs C VECTOR COMPARISON AT FIRST POINT (psi=1e-4, theta=0)")
print("=" * 100)

# First row (theta=0)
Q_row = Q_data[0]
Ccov_row = Ccov_data[0]
Cphy_row = Cphy_data[0]
Cmet_row = Cmet_data[0]

print("\nQ VECTOR (from gpeq_epf):")
print(f"  Q_psi  = {Q_row[4]:.6e} + {Q_row[5]:.6e}i")
print(f"  Q_theta = {Q_row[6]:.6e} + {Q_row[7]:.6e}i")

print("\nC COVARIANT (before adding j×ξ term):")
print(f"  C_psi   = {Ccov_row[4]:.6e} + {Ccov_row[5]:.6e}i")
print(f"  C_theta = {Ccov_row[6]:.6e} + {Ccov_row[7]:.6e}i")

print("\nC CONTRAVARIANT (Physics Method = Q + j×ξ):")
print(f"  C1_phy = {Cphy_row[4]:.6e} + {Cphy_row[5]:.6e}i")
print(f"  C2_phy = {Cphy_row[6]:.6e} + {Cphy_row[7]:.6e}i")
print(f"  C3_phy = {Cphy_row[8]:.6e} + {Cphy_row[9]:.6e}i")

print("\nC CONTRAVARIANT (Metric Method = g_ij * C_j):")
print(f"  C1_met = {Cmet_row[4]:.6e} + {Cmet_row[5]:.6e}i")
print(f"  C2_met = {Cmet_row[6]:.6e} + {Cmet_row[7]:.6e}i")
print(f"  C3_met = {Cmet_row[8]:.6e} + {Cmet_row[9]:.6e}i")

print("\n" + "=" * 100)
print("KEY OBSERVATION:")
print("=" * 100)
Q_theta = Q_row[6] + 1j*Q_row[7]
C2_phy = Cphy_row[6] + 1j*Cphy_row[7]

print(f"Q_theta    = {Q_theta}")
print(f"C2_physics = {C2_phy}")
print(f"C2 == Q_theta? {np.abs(C2_phy - Q_theta) < 1e-10}")
print("\n👉 C2_physics이 Q_theta와 같거나 거의 같은가?")
print("   만약 그렇다면, j×ξ term이 거의 0이라는 뜻!")
