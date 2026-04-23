#!/usr/bin/env python3
"""
Plot C² reconstruction diagnostics in r-z plane.

For each C² output file, creates a PNG with real and imaginary parts
displayed as scatter plots in the r-z grid.
Color scheme: red (positive) - white (0) - blue (negative)

Usage:
    python plot_c2_reconstruction.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import glob
import os
import sys


def create_rwb_colormap(n_bins=256):
    """
    Create a red-white-blue colormap.
    
    Parameters:
    -----------
    n_bins : int
        Number of bins in the colormap
    
    Returns:
    --------
    cmap : LinearSegmentedColormap
        Red-white-blue colormap (white at center for value 0)
    """
    colors_list = ['red', 'white', 'blue']
    cmap = LinearSegmentedColormap.from_list('rwb', colors_list, N=n_bins)
    return cmap


def plot_c2_reconstruction(filepath):
    """
    Plot C² reconstruction in r-z plane.
    
    Creates a PNG with real and imaginary parts as subplots.
    Each point is colored according to its C² value using red-white-blue
    colormap where white represents zero.
    
    Parameters:
    -----------
    filepath : str
        Path to the C² output file
        
    Returns:
    --------
    output_file : str
        Path to the generated PNG file
    """
    
    try:
        # Read the file - skip header line
        data = np.loadtxt(filepath, skiprows=1)
        
        if data.size == 0:
            print(f"Warning: {filepath} is empty")
            return None
        
        # Handle case where file has only one data point
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        # Extract columns: psi, theta, r, z, C2_re, C2_im
        psi = data[:, 0]
        theta = data[:, 1]
        r = data[:, 2]
        z = data[:, 3]
        c2_re = data[:, 4]
        c2_im = data[:, 5]
        
        # Create red-white-blue colormap
        cmap = create_rwb_colormap(256)
        
        # Calculate max absolute values for symmetric normalization
        vmax_re = np.max(np.abs(c2_re))
        vmax_im = np.max(np.abs(c2_im))
        
        # Use CenteredNorm with white at 0
        if vmax_re > 0:
            norm_re = mcolors.CenteredNorm(vcenter=0, halfrange=vmax_re)
        else:
            norm_re = mcolors.Normalize(vmin=-1, vmax=1)
            
        if vmax_im > 0:
            norm_im = mcolors.CenteredNorm(vcenter=0, halfrange=vmax_im)
        else:
            norm_im = mcolors.Normalize(vmin=-1, vmax=1)
        
        # Create figure with 2 subplots (real and imaginary)
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Extract quantity name from filename
        # e.g., gpec_recon_C2_cw_sol1.out -> C2_cw
        basename = os.path.basename(filepath)
        quantity = basename.replace('gpec_recon_', '').replace('.out', '').rsplit('_sol', 1)[0]
        
        # Plot real part
        scatter_re = axes[0].scatter(r, z, c=c2_re, cmap=cmap, norm=norm_re, 
                                     s=30, alpha=0.8)
        axes[0].set_xlabel('R (m)', fontsize=12, fontweight='bold')
        axes[0].set_ylabel('Z (m)', fontsize=12, fontweight='bold')
        axes[0].set_title(f'{quantity} - Real Part', fontsize=13, fontweight='bold')
        axes[0].grid(True, alpha=0.3, linestyle='--')
        axes[0].set_aspect('equal', adjustable='box')
        
        cbar_re = plt.colorbar(scatter_re, ax=axes[0], label='Value')
        cbar_re.set_label(f'{quantity} (Re)', fontsize=11, fontweight='bold')
        
        # Plot imaginary part
        scatter_im = axes[1].scatter(r, z, c=c2_im, cmap=cmap, norm=norm_im, 
                                     s=30, alpha=0.8)
        axes[1].set_xlabel('R (m)', fontsize=12, fontweight='bold')
        axes[1].set_ylabel('Z (m)', fontsize=12, fontweight='bold')
        axes[1].set_title(f'{quantity} - Imaginary Part', fontsize=13, fontweight='bold')
        axes[1].grid(True, alpha=0.3, linestyle='--')
        axes[1].set_aspect('equal', adjustable='box')
        
        cbar_im = plt.colorbar(scatter_im, ax=axes[1], label='Value')
        cbar_im.set_label(f'{quantity} (Im)', fontsize=11, fontweight='bold')
        
        # Add main title with statistics
        num_points = len(r)
        n_psi = len(np.unique(psi))
        fig.suptitle(f'{quantity} Reconstruction\n(n_psi={n_psi}, n_theta={num_points//n_psi if n_psi > 0 else 0})',
                     fontsize=14, fontweight='bold', y=1.00)
        
        plt.tight_layout()
        
        # Save as PNG
        output_file = filepath.replace('.out', '.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
        plt.close()
        
        return output_file
        
    except Exception as e:
        print(f"✗ Error processing {filepath}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """
    Main function: find and plot all C² output files.
    """
    
    # Find all C² output files
    c2_files = sorted(glob.glob('gpec_recon_C2_*.out'))
    
    if not c2_files:
        print("No C² output files found (gpec_recon_C2_*.out)")
        print(f"Current directory: {os.getcwd()}")
        return
    
    print(f"Found {len(c2_files)} C² output file(s)")
    print("-" * 60)
    
    successful = 0
    for filepath in c2_files:
        print(f"Processing: {filepath}")
        output_file = plot_c2_reconstruction(filepath)
        if output_file:
            successful += 1
    
    print("-" * 60)
    print(f"Completed: {successful}/{len(c2_files)} files processed successfully")


if __name__ == '__main__':
    main()
