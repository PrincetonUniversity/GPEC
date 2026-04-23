import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def load_gpec_xbrzphi(filepath):
    """
    Load GPEC XBRZPHI output file
    Returns: pandas DataFrame with all columns
    """
    # Skip header lines (first 6 lines contain metadata)
    df = pd.read_csv(filepath, sep=r'\s+', skiprows=6)
    return df

def plot_gpec_displacement_rzplane(filepath, figsize=(18, 12)):
    """
    Plot 6 R-Z plane subplots for displacement components (xr, xz, xp)
    Shows real and imaginary parts - 6 panels total
    """
    df = load_gpec_xbrzphi(filepath)
    
    # Displacement components - real and imaginary
    components = [
        ('real(xr) - Radial Displacement (Real)', df['real(xr)']),
        ('imag(xr) - Radial Displacement (Imag)', df['imag(xr)']),
        ('real(xz) - Vertical Displacement (Real)', df['real(xz)']),
        ('imag(xz) - Vertical Displacement (Imag)', df['imag(xz)']),
        ('real(xp) - Toroidal Displacement (Real)', df['real(xp)']),
        ('imag(xp) - Toroidal Displacement (Imag)', df['imag(xp)'])
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle('GPEC XBRZPHI: Displacement Components in R-Z Plane (n=1)', 
                 fontsize=16, fontweight='bold')
    
    axes_flat = axes.flatten()
    
    for idx, (title, values) in enumerate(components):
        ax = axes_flat[idx]
        
        # Create scatter plot colored by value
        scatter = ax.scatter(df['r'], df['z'], c=values, s=50, 
                           cmap='RdBu_r', alpha=0.8, edgecolors='none')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Value', fontsize=10)
        
        ax.set_xlabel('R (Major Radius)', fontsize=11)
        ax.set_ylabel('Z (Vertical Position)', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    return fig, df

def plot_gpec_magnetic_field_rzplane(filepath, figsize=(18, 12)):
    """
    Plot 6 R-Z plane subplots for magnetic field components (br, bz, bp)
    Shows real and imaginary parts - 6 panels total
    """
    df = load_gpec_xbrzphi(filepath)
    
    # Magnetic field components - real and imaginary
    components = [
        ('real(br) - Radial Magnetic Field (Real)', df['real(br)']),
        ('imag(br) - Radial Magnetic Field (Imag)', df['imag(br)']),
        ('real(bz) - Vertical Magnetic Field (Real)', df['real(bz)']),
        ('imag(bz) - Vertical Magnetic Field (Imag)', df['imag(bz)']),
        ('real(bp) - Toroidal Magnetic Field (Real)', df['real(bp)']),
        ('imag(bp) - Toroidal Magnetic Field (Imag)', df['imag(bp)'])
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle('GPEC XBRZPHI: Magnetic Field Components in R-Z Plane (n=1)', 
                 fontsize=16, fontweight='bold')
    
    axes_flat = axes.flatten()
    
    for idx, (title, values) in enumerate(components):
        ax = axes_flat[idx]
        
        # Create scatter plot colored by value
        scatter = ax.scatter(df['r'], df['z'], c=values, s=50, 
                           cmap='RdBu_r', alpha=0.8, edgecolors='none')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Value', fontsize=10)
        
        ax.set_xlabel('R (Major Radius)', fontsize=11)
        ax.set_ylabel('Z (Vertical Position)', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    return fig, df

# Usage
if __name__ == "__main__":
    filepath = './gpec_xbrzphi_fun_n1.out'
    
    # Plot 1: Displacement components (x) - 6 panels
    fig1, df = plot_gpec_displacement_rzplane(filepath)
    output_path_1 = './gpec_xbrzphi_displacement_6panels.png'
    fig1.savefig(output_path_1, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path_1}")
    
    # Plot 2: Magnetic field components (b) - 6 panels
    fig2, df = plot_gpec_magnetic_field_rzplane(filepath)
    output_path_2 = './gpec_xbrzphi_magnetic_6panels.png'
    fig2.savefig(output_path_2, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path_2}")
    
    # plt.show()