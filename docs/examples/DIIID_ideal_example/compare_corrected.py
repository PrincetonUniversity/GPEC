#!/usr/bin/env python3
import numpy as np

# Load both methods
phys_data = np.loadtxt("gpec_recon_Ccontra_physics_sol1.out")
metric_data = np.loadtxt("gpec_recon_Ccontra_metric_sol1.out")

print("=" * 80)
print("Sample comparison at different flux surfaces (theta=0)")
print("=" * 80)

# Get unique psi values
psi_values = phys_data[:, 0]
unique_psi = np.unique(psi_values)

for i in [0, len(unique_psi)//4, len(unique_psi)//2, 3*len(unique_psi)//4, -1]:
    psi_val = unique_psi[i]
    idx_phys = np.where((phys_data[:, 0] == psi_val) & (phys_data[:, 1] == 0))[0]
    idx_metric = np.where((metric_data[:, 0] == psi_val) & (metric_data[:, 1] == 0))[0]
    
    if len(idx_phys) > 0 and len(idx_metric) > 0:
        phys_row = phys_data[idx_phys[0]]
        metric_row = metric_data[idx_metric[0]]
        
        c1p = phys_row[4] + 1j*phys_row[5]
        c2p = phys_row[6] + 1j*phys_row[7]
        c3p = phys_row[8] + 1j*phys_row[9]
        
        c1m = metric_row[4] + 1j*metric_row[5]
        c2m = metric_row[6] + 1j*metric_row[7]
        c3m = metric_row[8] + 1j*metric_row[9]
        
        print(f"\nψ = {psi_val:.6e}:")
        print(f"  Physics: C1={abs(c1p):.3e}, C2={abs(c2p):.3e}, C3={abs(c3p):.3e}")
        print(f"  Metric:  C1={abs(c1m):.3e}, C2={abs(c2m):.3e}, C3={abs(c3m):.3e}")
        if abs(c1m) > 1e-15:
            print(f"  Ratio C2_phy/C1_met = {abs(c2p)/abs(c1m):.4f}")
