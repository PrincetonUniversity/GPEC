# GPEC - Generalized Perturbed Equilibrium Code

## Overview

GPEC is a comprehensive suite of nonaxisymmetric stability and perturbed equilibrium codes for tokamak plasma analysis. It calculates stability of tokamak plasmas to nonaxisymmetric modes and performs nonaxisymmetric force balance calculations for stable configurations.

## Primary Languages

- **Fortran 90** (primary) - Core physics codes with modern features (modules, derived types, dynamic allocation)
- **Python** - PyPEC package for postprocessing, running simulations, and visualization

## Key Components/Codes

| Code | Purpose |
|------|---------|
| **DCON** | Direct Criterion of Newcomb - ideal MHD stability analysis |
| **STRIDE** | Optimized DCON with parallel computation support |
| **GPEC** | Generalized Perturbed Equilibrium Code - adds forcing terms from external 3D fields |
| **RDCON** | Resistive DCON - resistive MHD stability |
| **PENTRC** | Kinetic effects analysis (torque matrix calculations) |
| **RMATCH** | Resistive matching for inner-layer solutions |
| **ORBIT** | Full-orbit charged particle tracking |
| **VACUUM** | Vacuum response calculations |
| **SLAYER** | Slab layer model based on linear drift MHD |

## Directory Structure

```
gpec/
├── bin/          # Compiled executables
├── coil/         # Coil geometry and magnetic field calculations
├── dcon/         # Direct Criterion of Newcomb (ideal MHD)
├── equil/        # Equilibrium data processing (EFIT, CHEASE, Miller, etc.)
├── gpec/         # Main GPEC code for 3D perturbations
├── docs/         # Sphinx documentation and examples
├── harvest/      # Git submodule for data client
├── input/        # Template input namelists
├── install/      # Build configuration and Makefile
├── lib/          # Compiled library archives
├── lsode/        # LSODE ODE solver
├── match/        # Resistive matching
├── multi/        # Multiple run processing
├── orbit/        # Particle orbit tracking
├── pentrc/       # Kinetic effects/torque matrix
├── pypec/        # Python analysis package
├── rdcon/        # Resistive DCON
├── rmatch/       # Resistive matching extended
├── slayer/       # Slab layer MHD model
├── stride/       # STRIDE stability code
├── sum/          # Data extraction from MULTI runs
├── vacuum/       # Vacuum response calculations
├── xdraw/        # X11 graphics visualization
├── zlange/       # Complex linear algebra
└── zvode/        # Complex ODE solver
```

## Building

```bash
cd install
make [FFLAGS="flags"] [CC="compiler"] [FC="fortran_compiler"]
```

### Dependencies
- Fortran 90 compiler (ifort, gfortran, pgfortran)
- LAPACK/BLAS (or Intel MKL)
- NetCDF (Fortran and C)
- X11 (for XDRAW graphics)

### Supported Platforms
- PPPL (portal.pppl.gov)
- GA (iris.gat.com)
- NERSC Perlmutter
- macOS (Intel and ARM with gfortran)
- Linux (gfortran, Ubuntu/WSL)

## Key Input Files (Fortran Namelists)

| File | Purpose |
|------|---------|
| `equil.in` | Equilibrium control (format, grid, coordinates) |
| `dcon.in` | DCON/STRIDE stability analysis parameters |
| `vac.in` | Vacuum response calculation settings |
| `gpec.in` | Nonaxisymmetric perturbation configuration |
| `rdcon.in` | Resistive analysis parameters |
| `pentrc.in` | Kinetic torque calculations |
| `coil.in` | Coil geometry specifications |

## Equilibrium Formats Supported

- EFIT (GA code) - g-files
- CHEASE (Lausanne)
- Miller (inverse)
- TRANSP data
- JSOLVER
- Analytical (LAR, Soloviev)

## Coordinate Systems

- Hamada (default)
- PEST
- Boozer
- Equal-arc

## PyPEC (Python Package)

Located in `pypec/` directory:
- `gpec.py` - Main runner wrapper
- `data.py` - Data file reading
- `namelist.py` - Fortran namelist handling
- `modplot.py` - Plotting utilities
- `gui.py` - GUI interface (Enthought Traits-based)
- `synthetics.py` - Synthetic data generation

## Common Workflows

### Running DCON/STRIDE (Ideal MHD Stability)
1. Prepare equilibrium file and `equil.in`
2. Configure `dcon.in` with mode numbers and flags
3. Run `dcon` or `stride` executable
4. View results with `xdraw` or analyze NetCDF outputs

### Running GPEC (3D Perturbations)
1. First run DCON/STRIDE to generate `euler.bin`
2. Configure `gpec.in` with 3D field source (coil, harmonic, or data file)
3. Run `gpec` executable
4. Analyze outputs for torque, stability metrics

### Running Resistive Analysis
1. Run DCON/STRIDE first
2. Configure `rdcon.in`
3. Run `rdcon`, then `rmatch` if needed

## Key Output Files

- `*.bin` - Binary files viewable with XDRAW
- `*.nc` - NetCDF files for analysis
- `euler.bin` - Orthonormal basis from DCON (input to GPEC)
- `globalsol.bin` - Resistive solution from RMATCH

## Module System Usage

```bash
module load gpec           # Default release
module load gpec/0.0       # Development branch
```

Sets `$GPECHOME` environment variable.

## Documentation

Sphinx documentation in `docs/`:
- `docs/index.rst` - Main documentation index
- `docs/outputs.rst` - Output specification guide
- `docs/releases.rst` - Version history
- `docs/examples/` - Example cases and tutorials

Build docs: `cd docs && make html`

## Physics Features

- Mercier criterion (D_I, D_R) for interchange modes
- Ballooning criterion (C_A) for ideal MHD
- Newcomb criterion for low-n modes
- Plasma, vacuum, and total response matrices
- Resistive instability detection
- Kinetic effects (torque matrices)
- Full-orbit and guiding-center particle tracking
