#!/usr/bin/env python3
"""
Plot reconstruction diagnostics in r-z plane.

For each reconstruction output file, creates a PNG with real and imaginary parts
displayed as scatter plots in the r-z grid.
Color scheme: red (positive) - white (0) - blue (negative)

Handles different file formats:
- Scalar real-only (5 columns): 1 subplot
- Scalar complex (6 columns): 2 subplots
- 2-component vectors (8 columns): 2 subplots
- 3-component vectors (10 columns): 3 subplots

Usage:
    python plot_reconstruction.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import glob
import os


def create_rwb_colormap(n_bins=256):
    """
    Create a red-white-blue colormap where white is centered at 0.
    
    Parameters:
    -----------
    n_bins : int
        Number of bins in the colormap
    
    Returns:
    --------
    cmap : LinearSegmentedColormap
        Red-white-blue colormap
    """
    colors_list = ['red', 'white', 'blue']
    cmap = LinearSegmentedColormap.from_list('rwb', colors_list, N=n_bins)
    return cmap


def get_component_names(filepath, n_components):
    """
    Extract component names from file based on filename and number of components.
    
    Parameters:
    -----------
    filepath : str
        Path to the output file
    n_components : int
        Number of components (1, 2, or 3)
        
    Returns:
    --------
    names : list
        List of component names
    """
    basename = os.path.basename(filepath)
    quantity = basename.replace('gpec_recon_', '').replace('.out', '').rsplit('_sol', 1)[0]
    
    if n_components == 1:
        return [quantity]
    elif n_components == 2:
        # For 2D vectors (like C_psi, C_theta)
        if 'Ccov' in basename:
            return ['ψ component (C_psi)', 'θ component (C_theta)']
        elif 'Qv' in basename:
            return ['ψ component (Qv_psi)', 'θ component (Qv_theta)']
        elif 'Qw' in basename:
            return ['ψ component (Qw_psi)', 'θ component (Qw_theta)']
        else:
            return ['Component 1', 'Component 2']
    elif n_components == 3:
        # For 3D vectors (like C1, C2, C3 or Psi, Theta, Zeta)
        if 'Ccov' in basename:
            return ['ψ component', 'θ component', 'ζ component']
        elif 'Qv' in basename:
            return ['ψ component', 'θ component', 'ζ component']
        elif 'Qw' in basename:
            return ['ψ component', 'θ component', 'ζ component']
        elif 'Ccontra' in basename:
            return ['C¹ component', 'C² component', 'C³ component']
        else:
            return [f'Component {i+1}' for i in range(n_components)]
    else:
        return [f'Component {i+1}' for i in range(n_components)]


def plot_reconstruction(filepath):
    """
    Plot reconstruction data in r-z plane.
    
    Creates PNG files with real and imaginary parts as subplots.
    
    Parameters:
    -----------
    filepath : str
        Path to the reconstruction output file
        
    Returns:
    --------
    output_files : list
        List of generated PNG file paths
    """
    
    try:
        # Read the file with error handling for inconsistent columns
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Skip header and filter out incomplete lines (blank lines, comment lines)
        data_lines = []
        for line in lines[1:]:  # Skip header
            line = line.strip()
            if not line or line.startswith('#'):  # Skip blank lines
                continue
            # Count columns: split and filter empty strings
            cols = [x for x in line.split() if x]
            # Only keep lines with full column count.
            # Supported formats:
            #   5 cols  = psi theta r z scalar_re
            #   6 cols  = psi theta r z scalar_re scalar_im
            #   8 cols  = psi theta r z vec2_re/im
            #   10 cols = psi theta r z vec3_re/im
            if len(cols) in [5, 6, 8, 10]:
                data_lines.append(line)
        
        if not data_lines:
            print(f"Warning: {filepath} has no valid data lines")
            return []
        
        # Convert to array
        data_str = '\n'.join(data_lines)
        from io import StringIO
        data = np.loadtxt(StringIO(data_str))
        
        if data.size == 0:
            print(f"Warning: {filepath} is empty")
            return []
        
        # Handle case where file has only one data point
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        # Extract columns
        psi = data[:, 0]
        theta = data[:, 1]
        r = data[:, 2]
        z = data[:, 3]
        
        n_cols = data.shape[1]
        scalar_real_only = (n_cols == 5)

        if scalar_real_only:
            n_components = 1
            quantities = [{'re': data[:, 4], 'im': None}]
        else:
            n_data_cols = n_cols - 4  # Subtract psi, theta, r, z
            n_components = n_data_cols // 2  # Each component has (re, im)

            # Extract quantities
            quantities = []
            for i in range(n_components):
                re_col = 4 + i * 2
                im_col = 4 + i * 2 + 1
                quantities.append({
                    're': data[:, re_col],
                    'im': data[:, im_col]
                })
        
        # Get quantity name and component names
        basename = os.path.basename(filepath)
        quantity = basename.replace('gpec_recon_', '').replace('.out', '').rsplit('_sol', 1)[0]
        comp_names = get_component_names(filepath, n_components)
        
        # Create red-white-blue colormap
        cmap = create_rwb_colormap(256)
        
        if scalar_real_only:
            fig, ax = plt.subplots(1, 1, figsize=(8, 7))
            axes = np.array([[ax]])
        else:
            # Create figure with subplots (n_components rows x 2 columns for re/im)
            n_rows = n_components
            fig, axes = plt.subplots(n_rows, 2, figsize=(15, 6*n_rows))

            # Handle single row case (axes not 2D)
            if n_rows == 1:
                axes = axes.reshape(1, -1)
        
        output_files = []
        
        # Plot each component
        for i, (comp_data, comp_name) in enumerate(zip(quantities, comp_names)):
            c2_re = comp_data['re']
            c2_im = comp_data['im']

            # Calculate max absolute values for symmetric normalization
            vmax_re = np.max(np.abs(c2_re))
            if c2_im is not None:
                vmax_im = np.max(np.abs(c2_im))

            # Use CenteredNorm with white at 0
            if vmax_re > 0:
                norm_re = mcolors.CenteredNorm(vcenter=0, halfrange=vmax_re)
            else:
                norm_re = mcolors.Normalize(vmin=-1, vmax=1)

            # Plot real part
            scatter_re = axes[i, 0].scatter(
                r, z, c=c2_re, cmap=cmap, norm=norm_re, s=30, alpha=0.8
            )
            axes[i, 0].set_xlabel('R (m)', fontsize=11, fontweight='bold')
            axes[i, 0].set_ylabel('Z (m)', fontsize=11, fontweight='bold')
            real_title = f'{quantity} - {comp_name}' if scalar_real_only \
                else f'{quantity} - {comp_name} (Real)'
            axes[i, 0].set_title(real_title, fontsize=12, fontweight='bold')
            axes[i, 0].grid(True, alpha=0.3, linestyle='--')
            axes[i, 0].set_aspect('equal', adjustable='box')

            cbar_re = plt.colorbar(scatter_re, ax=axes[i, 0])
            cbar_re.set_label('Value', fontsize=10, fontweight='bold')

            if c2_im is not None:
                if vmax_im > 0:
                    norm_im = mcolors.CenteredNorm(vcenter=0, halfrange=vmax_im)
                else:
                    norm_im = mcolors.Normalize(vmin=-1, vmax=1)

                scatter_im = axes[i, 1].scatter(
                    r, z, c=c2_im, cmap=cmap, norm=norm_im, s=30, alpha=0.8
                )
                axes[i, 1].set_xlabel('R (m)', fontsize=11, fontweight='bold')
                axes[i, 1].set_ylabel('Z (m)', fontsize=11, fontweight='bold')
                axes[i, 1].set_title(
                    f'{quantity} - {comp_name} (Imaginary)',
                    fontsize=12,
                    fontweight='bold'
                )
                axes[i, 1].grid(True, alpha=0.3, linestyle='--')
                axes[i, 1].set_aspect('equal', adjustable='box')

                cbar_im = plt.colorbar(scatter_im, ax=axes[i, 1])
                cbar_im.set_label('Value', fontsize=10, fontweight='bold')
        
        # Add main title with statistics
        num_points = len(r)
        n_psi = len(np.unique(psi))
        n_theta = num_points // n_psi if n_psi > 0 else 0
        
        fig.suptitle(f'{quantity} Reconstruction (n_psi={n_psi}, n_theta={n_theta})',
                     fontsize=14, fontweight='bold', y=0.995)
        
        plt.tight_layout()
        
        # Save as PNG
        output_file = filepath.replace('.out', '.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✓ Saved: {output_file}")
        plt.close()
        
        output_files.append(output_file)
        return output_files
        
    except Exception as e:
        print(f"✗ Error processing {filepath}: {e}")
        import traceback
        traceback.print_exc()
        return []


def main():
    """
    Main function: find and plot all reconstruction output files.
    """
    
    # Find all reconstruction output files (excluding C2 which already has its own script)
    recon_files = sorted(glob.glob('gpec_recon_*.out'))
    
    # Filter out C2 files (they have their own dedicated script)
    recon_files = [f for f in recon_files if 'C2_' not in f]
    
    if not recon_files:
        print("No reconstruction output files found")
        print(f"Current directory: {os.getcwd()}")
        return
    
    print(f"Found {len(recon_files)} reconstruction output file(s)")
    print("-" * 70)
    
    successful = 0
    total_files = 0
    
    for filepath in recon_files:
        print(f"Processing: {filepath}")
        output_files = plot_reconstruction(filepath)
        if output_files:
            successful += 1
            total_files += len(output_files)
    
    print("-" * 70)
    print(f"Completed: {successful}/{len(recon_files)} files processed successfully")
    print(f"Generated {total_files} PNG file(s)")


if __name__ == '__main__':
    main()
