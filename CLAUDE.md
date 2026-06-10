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



## Claude Coding Guidelines for GPEC

### Fortran Fixed-Format (.f files) Requirements

**CRITICAL: Line Width Limit**
- Fixed-format Fortran files (.f, NOT .f90) have a **72-character line limit**
- Lines exceeding 72 characters will cause compilation errors or be silently truncated
- **ALWAYS** check line length before making edits to .f files

#### Line Format Rules
- Columns 1-5: Statement labels (optional)
- Column 6: Continuation character (use `$` or `&` for continuation lines)
- Columns 7-72: Fortran statements
- Columns 73-80: Ignored (historically used for sequence numbers)

#### Continuation Lines
When a statement exceeds 72 characters, split it across multiple lines:
```fortran
      variable_name = very_long_expression_that_would_exceed_limit
     $                + more_of_the_expression
```

#### String Splitting
For long strings, use Fortran string concatenation:
```fortran
      WRITE(*,*) "This is a very long warning message that "//
     $           "needs to be split across multiple lines"
```

### Type Precision

**REAL Variables:**
- Use `0.0_r8` instead of `0.0` for REAL(r8) variables
- Use `1.0_r8` instead of `1.0` for REAL(r8) variables
- Type suffix `_r8` ensures proper precision matching

### Variable Naming Conventions

Use consistent naming patterns within subroutines:
- `hw_*` prefix for half-widths (e.g., `hw_isl`, `hw_v`, `hw_v_crit`, `hw_sat`, `hw_min`)
- Avoid mixing naming styles within the same subroutine

### Safety Checks

**Array Access:**
- Always check if allocatable arrays are ALLOCATED before accessing them:
```fortran
      IF (ALLOCATED(array_name)) THEN
         ! Use array
      ELSE
         ! Handle unallocated case
      ENDIF
```

**Division:**
- Guard against division by zero for potentially small values
- Check bounds for ACOS, ASIN arguments (must be in [-1, 1])
- Check that square root arguments are non-negative

### Before Committing

**Pre-commit Checklist:**
1. Check all modified .f files for lines exceeding 72 characters
2. Verify type suffixes match variable declarations (_r8 for REAL(r8))
3. Ensure continuation lines use proper column 6 markers ($)
4. Test compilation before committing
5. Ask user permission before committing changes

### Verification Commands

Check for lines exceeding 72 characters:
```bash
awk 'length > 72 {print NR": "length" chars"}' filename.f
```

Count characters in a specific line:
```bash
sed -n 'NUMp' filename.f | wc -c
```

## Sign Conventions

GPEC uses right-handed magnetic coordinates (psi, theta, zeta) with Fourier
decomposition exp(im*theta - in*phi). A comprehensive reference is in
`docs/sign_conventions.rst`. Key points:

- **helicity** = ipd * btd (+1 = RH, -1 = LH), computed in `gpec_main` (`gpec/gpec.f`)
- **nn** (toroidal mode number) is always positive; resonant **m** is always positive
- **F = R*Bt** is forced positive via ABS() in `read_eq_efit` (`equil/read_eq.f`)
- **q** is recomputed by field-line integration in `direct_run` (`equil/direct.f`); the g-file q profile is unused, and q is always positive for g-file input
- **omega_E** is positive for rotation in the direction of the toroidal coordinate zeta
- The code does NOT use the COCOS standard
- For SURFMN interface: `m_surfmn = helicity * m_gpec`
- For real-space output: RH configs take complex conjugate (`-helicity*AIMAG(...)` in `gpec/gpout.f`)
