import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from pathlib import Path

class GPECNetCDFViewer:
    """GPEC NetCDF File Viewer and Plotter"""
    
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.ds = None
        self.load_file()
        
    def load_file(self):
        try:
            self.ds = xr.open_dataset(self.filepath)
            print(f"✓ Loaded: {self.filepath.name}")
        except Exception as e:
            print(f"✗ Error: {e}")
            
    def plot_bkappaxi_perp(self):
        """Plot Bkappaxi_perp in Fourier space"""
        if 'Bkappaxi_perp' not in self.ds:
            print("✗ 'Bkappaxi_perp' not found")
            return
        
        data = self.ds['Bkappaxi_perp'].values
        real_part = data[:, :, 0]
        imag_part = data[:, :, 1]
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle('Bkappaxi_perp: Curvature Component (Fourier Space)', 
                     fontsize=14, fontweight='bold')
        
        im1 = axes[0].contourf(real_part, levels=20, cmap='RdBu_r')
        axes[0].set_title('Real Part')
        plt.colorbar(im1, ax=axes[0])
        
        im2 = axes[1].contourf(imag_part, levels=20, cmap='RdBu_r')
        axes[1].set_title('Imaginary Part')
        plt.colorbar(im2, ax=axes[1])
        
        plt.tight_layout()
        return fig
    
    def plot_bkappaxi_perp_fun(self):
        """Plot Bkappaxi_perp_fun in function space"""
        if 'Bkappaxi_perp_fun' not in self.ds:
            print("✗ 'Bkappaxi_perp_fun' not found")
            return
        
        data = self.ds['Bkappaxi_perp_fun'].values
        real_part = data[:, :, 0]
        imag_part = data[:, :, 1]
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle('Bkappaxi_perp_fun: Curvature Component (Function Space)', 
                     fontsize=14, fontweight='bold')
        
        im1 = axes[0].contourf(real_part, levels=20, cmap='RdBu_r')
        axes[0].set_title('Real Part')
        plt.colorbar(im1, ax=axes[0])
        
        im2 = axes[1].contourf(imag_part, levels=20, cmap='RdBu_r')
        axes[1].set_title('Imaginary Part')
        plt.colorbar(im2, ax=axes[1])
        
        plt.tight_layout()
        return fig

# Use the viewer
if __name__ == "__main__":
    filepath = "gpec_output_n1.nc"  # 실제 파일명으로 변경
    viewer = GPECNetCDFViewer(filepath)
    
    if viewer.ds is not None:
        fig1 = viewer.plot_bkappaxi_perp()
        plt.savefig('bkappaxi_perp_fourier.png', dpi=150, bbox_inches='tight')
        
        fig2 = viewer.plot_bkappaxi_perp_fun()
        plt.savefig('bkappaxi_perp_function.png', dpi=150, bbox_inches='tight')
        
        plt.show()
