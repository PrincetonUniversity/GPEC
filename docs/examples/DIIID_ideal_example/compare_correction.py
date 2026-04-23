#!/usr/bin/env python3
import numpy as np

# Load physics method
phys_data = np.loadtxt("gpec_recon_Ccontra_physics_sol1.out")
phys_c1 = phys_data[:, 4] + 1j * phys_data[:, 5]
phys_c2 = phys_data[:, 6] + 1j * phys_data[:, 7]
phys_c3 = phys_data[:, 8] + 1j * phys_data[:, 9]

# Load metric method
metric_data = np.loadtxt("gpec_recon_Ccontra_metric_sol1.out")
metric_c1 = metric_data[:, 4] + 1j * metric_data[:, 5]
metric_c2 = metric_data[:, 6] + 1j * metric_data[:, 7]
metric_c3 = metric_data[:, 8] + 1j * metric_data[:, 9]

print("=" * 70)
print("CORRECTION: sq.f1(2) and f1(3) already contain μ₀ and 2π respectively")
print("=" * 70)
print("\nUnderneath first flux surface (ψ = psi_low):")
print("\nPhysics Method Results:")
print(f"  C1_physics = {phys_c1[0]:.6e}")
print(f"  C2_physics = {phys_c2[0]:.6e}")
print(f"  C3_physics = {phys_c3[0]:.6e}")

print("\nMetric Method Results:")
print(f"  C1_metric  = {metric_c1[0]:.6e}")
print(f"  C2_metric  = {metric_c2[0]:.6e}")
print(f"  C3_metric  = {metric_c3[0]:.6e}")

print("\n" + "=" * 70)
print("COMPARISON (Ratio of Physics / Metric):")
print("=" * 70)
c1_ratio = abs(phys_c1[0]) / abs(metric_c1[0])
c2_ratio = abs(phys_c2[0]) / abs(metric_c2[0])
c3_ratio = abs(phys_c3[0]) / abs(metric_c3[0])

print(f"  C1_physics / C1_metric = {c1_ratio:.6f}")
print(f"  C2_physics / C2_metric = {c2_ratio:.6f}")
print(f"  C3_physics / C3_metric = {c3_ratio:.6f}")

print(f"\n  2π ≈ {2 * 3.14159:.6f}")
print(f"\n✓ Physics method should now ≈ Metric method (within numerical error)")
