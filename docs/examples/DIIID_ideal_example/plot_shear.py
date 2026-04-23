#!/usr/bin/env python3
"""
Plot shear diagnostic quantities from GPEC reconstruction output.

Reads spatial domain data:
- gpec_recon_shear_fun_sol*.out: Spatial domain (r-z plane)

Produces single r-z plane plot with shear values colored.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.interpolate import griddata
import sys
import glob

def find_shear_fun_file():
    """Find shear_fun output file with wildcard matching."""
    fun_files = glob.glob('gpec_recon_shear_fun_sol*.out')
    
    if not fun_files:
        print("Error: Could not find gpec_recon_shear_fun_sol*.out file")
        sys.exit(1)
    
    return fun_files[0]

def read_shear_fun_data():
    """Read spatial shear data from output file."""
    fun_file = find_shear_fun_file()
    
    print(f"Reading spatial data from: {fun_file}")
    shear_fun = np.loadtxt(fun_file, skiprows=1)
    
    # Parse spatial data
    psi_fun = shear_fun[:, 0]
    theta_fun = shear_fun[:, 1]
    r_fun = shear_fun[:, 2]
    z_fun = shear_fun[:, 3]
    shear_re_fun = shear_fun[:, 4]
    shear_im_fun = shear_fun[:, 5]
    
    return psi_fun, theta_fun, r_fun, z_fun, shear_re_fun, shear_im_fun

def main():
    """Main plotting routine."""
    print("=" * 70)
    print("GPEC Shear Diagnostic Plotter (r-z plane only)")
    print("=" * 70)
    
    # Read data
    psi_fun, theta_fun, r_fun, z_fun, shear_re_fun, shear_im_fun = read_shear_fun_data()
    
    # Get unique values
    unique_psi_fun = np.unique(psi_fun)
    
    print(f"Spatial grid: {len(unique_psi_fun)} psi levels")
    print(f"Total grid points: {len(r_fun)}")
    
    # Create figure with single r-z plane plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Apply asinh transformation to shear
    shear_re_fun_asinh = np.arcsinh(shear_re_fun)
    
    # Create grid for interpolation
    r_min, r_max = np.min(r_fun), np.max(r_fun)
    z_min, z_max = np.min(z_fun), np.max(z_fun)
    
    grid_r = np.linspace(r_min, r_max, 150)
    grid_z = np.linspace(z_min, z_max, 150)
    grid_r_mesh, grid_z_mesh = np.meshgrid(grid_r, grid_z)
    
    # Interpolate shear values onto grid
    points = np.c_[r_fun, z_fun]
    grid_shear = griddata(points, shear_re_fun_asinh, (grid_r_mesh, grid_z_mesh), 
                          method='cubic', fill_value=np.nan)
    
    # Normalize with 0 at white
    vmin, vmax = np.nanmin(grid_shear), np.nanmax(grid_shear)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    
    # Plot shear as filled contours with scatter points
    contourf = ax.contourf(grid_r_mesh, grid_z_mesh, grid_shear, levels=30, 
                           cmap='RdBu_r', norm=norm)
    
    # Overlay scatter plot for actual data points
    scatter = ax.scatter(r_fun, z_fun, c=shear_re_fun_asinh, cmap='RdBu_r', 
                        norm=norm, s=15, alpha=0.6, edgecolors='none')
    
    # Draw contour at shear=0 in black
    contour_zero = ax.contour(grid_r_mesh, grid_z_mesh, grid_shear, 
                              levels=[0], colors='black', linewidths=2)
    
    ax.set_xlabel('R', fontsize=22, fontweight='bold')
    ax.set_ylabel('Z', fontsize=22, fontweight='bold')
    ax.set_title('Shear (spatial) - r-z plane', fontsize=24, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(labelsize=16)
    
    # Set r-axis ticks
    r_ticks = np.linspace(r_min, r_max, 6)
    ax.set_xticks(r_ticks)
    ax.set_xticklabels([f'{r:.3f}' for r in r_ticks], fontsize=14)
    
    # Set z-axis ticks
    z_ticks = np.linspace(z_min, z_max, 6)
    ax.set_yticks(z_ticks)
    ax.set_yticklabels([f'{z:.3f}' for z in z_ticks], fontsize=14)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('asinh(Shear)', fontsize=18, fontweight='bold')
    cbar.ax.tick_params(labelsize=14)
    
    # Set colorbar ticks to span negative and positive values
    cbar_ticks = np.linspace(vmin, vmax, 7)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f'{val:.1e}' for val in cbar_ticks], fontsize=12)
    
    # Adjust layout and save
    plt.tight_layout()
    
    output_file = 'gpec_recon_shear_rz.png'
    print(f"\nSaving plot to: {output_file}")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print("Plot complete!")
    

if __name__ == '__main__':
    main()
