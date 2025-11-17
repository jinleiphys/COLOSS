# COLOSS
**COLOSS: Complex-scaled Optical and couLOmb Scattering Solver** is a program designed to address the scattering problem using a bound-state technique known as complex scaling. In this method, the oscillatory boundary conditions of the wave function are transformed into exponentially decaying ones, accommodating the long-range Coulomb interaction. This program implements the general local optical model potential and the Perey-Buck nonlocal optical potential, with all parameters included in a well-designed input format for ease of use. This design offers users straightforward access to compute **S-matrices** and **cross-sections** of the scattering process.

The program employs **Lagrange-Laguerre functions** as the basis functions and provides two methods to perform numerical integration: one based on the approximation of Lagrange functions, and the other using direct **Gauss-Legendre quadrature**. Additionally, COLOSS incorporates two distinct rotation methods, making it adaptable to potentials without analytical expressions.

For those interested in the theoretical foundations and implementation details of COLOSS, a comprehensive **published paper** is available at [https://www.sciencedirect.com/science/article/abs/pii/S0010465525000712?via%3Dihub](https://www.sciencedirect.com/science/article/abs/pii/S0010465525000712?via%3Dihub). Furthermore, COLOSS now includes a modern **Python GUI application** (see below) as well as a web-based **UI interface input file generator** named `coloss-input-generator.html`, both designed to facilitate the preparation of input files and ensure a smooth start to your computational tasks.

## Table of Contents
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Graphical User Interface (GUI)](#graphical-user-interface-gui)
- [Input Description](#input-description)
- [Output Description](#output-description)
- [Examples](#examples)

---

## Project Structure

```
COLOSS/
├── src/                    # Fortran source files
│   ├── COLOSS.f           # Main program
│   ├── input.f            # Input handling and namelist reading
│   ├── system.f           # System parameters
│   ├── pot_class.f        # Potential parameter definitions
│   ├── rot_potential.f    # Rotated potential calculations
│   ├── scatt.f            # Scattering calculations
│   ├── bound.f            # Bound state calculations
│   ├── matrix_element.f   # Matrix element calculations
│   ├── mesh.f             # Mesh generation
│   ├── generate_laguerre.f # Laguerre basis functions
│   └── ...                # Other modules
├── build/                  # Build artifacts (*.o, *.mod files)
├── adyo_v1_0/             # C++ Coulomb wave function library
│   └── libcwf_cpp.a       # Compiled C++ library
├── coloss_gui/            # Python GUI application
│   ├── main.py            # GUI entry point
│   ├── main_window.py     # Main window
│   ├── input_panel.py     # Input parameter panel
│   ├── plot_widget.py     # Result plotting
│   ├── log_widget.py      # Output logging
│   ├── runner.py          # COLOSS execution wrapper
│   ├── styles.py          # Modern UI styling
│   └── path_utils.py      # Path detection utilities
├── test/                   # Example input files
│   ├── *.in               # Various test cases
│   └── fort.*             # Output files
├── COLOSS                  # Compiled executable (after build)
├── Makefile               # Build configuration
├── compile.sh             # Automated compilation script
├── run_coloss_gui.sh      # GUI launcher script
└── README.md              # This file
```

---

## Getting Started

### Prerequisites

**COLOSS** uses the GNU Compiler Collection (GCC) for Fortran compilation. Make sure you have:
- GCC with gfortran (Fortran compiler)
- LAPACK library for linear algebra
- Conda or Miniconda (for GUI, will be installed automatically by compile.sh if not present)

### Quick Start: Using compile.sh

The easiest way to compile COLOSS and set up the GUI is to use the provided `compile.sh` script:

```bash
./compile.sh
```

#### What compile.sh does:

1. **Checks for required tools**: Verifies gfortran and make are installed
2. **Compiles the Fortran code**: Runs `make` to build the COLOSS executable
3. **Sets up Python environment**:
   - Checks for conda/miniconda (prompts to install if missing)
   - Creates a conda environment named `coloss_gui`
   - Installs PySide6 and other Python dependencies
4. **Creates launcher script**: Generates `run_coloss_gui.sh` for easy GUI access

#### Executable Location

After successful compilation:
- **COLOSS executable**: `./COLOSS` (in the root directory)
- **GUI launcher**: `./run_coloss_gui.sh`

### Manual Compilation (Alternative)

If you prefer to compile manually:

```bash
# Edit Makefile to specify your LAPACK path
# Then run:
make

# To clean build artifacts:
make clean
```

The Makefile is configured to:
- Compile source files from `src/` directory
- Place object files and modules in `build/` directory
- Generate the `COLOSS` executable in the root directory

### Running COLOSS (Command Line)

Transfer to a working directory (e.g., `test/`) and run with standard input:

```bash
cd test
../COLOSS < inputfile.in
```

Or pipe the input directly:

```bash
./COLOSS < test/example.in
```

---

## Graphical User Interface (GUI)

COLOSS includes a modern graphical user interface for easier parameter setup and visualization.

### Launching the GUI

After running `compile.sh`, start the GUI with:

```bash
./run_coloss_gui.sh
```

Or manually activate the conda environment:

```bash
conda activate coloss_gui
python coloss_gui/main.py
```

### GUI Features

- **Tabbed Interface**: Four organized tabs for different parameter groups
  - **General**: Mesh parameters, rotation angle, calculation methods
  - **System**: Projectile/target properties, energy, angular momentum
  - **Potential**: Optical model potential parameters (volume, surface, spin-orbit, Coulomb)
  - **Nonlocal**: Nonlocal potential settings

- **Parameter Tooltips**: Hover over info icons (ⓘ) for detailed parameter descriptions
  - Click the info icon to see full documentation in a popup window

- **Smart Defaults**: Auto-updating values (e.g., Gauss Rmax follows Rmax)

- **Input Validation**: Prevents invalid configurations (e.g., Method 2 with nonlocal)

- **External Potential Support**: File dialog for uploading potential data files

- **Real-time Visualization**:
  - Live output log during calculation
  - Automatic plotting of results (S-matrix, cross sections)

- **Modern Design**: Beautiful gradient-based UI with light/dark theme options

### GUI Workflow

1. **Set Parameters**: Use the tabbed interface to configure your calculation
2. **Save Input**: Save parameters to a `.in` file (optional)
3. **Run Calculation**: Click "Run" to execute COLOSS
4. **View Results**:
   - Monitor progress in the "Output Log" tab
   - View plots in the "Plot" tab after completion

---

## Input Description

COLOSS utilizes the FORTRAN namelist format for user-friendly input. All parameters are optional with sensible defaults.

### General namelist
Parameters for mesh, rotation, and numerical methods:

- **nr** (integer): Number of Lagrange-Laguerre basis functions / order of generalized Laguerre polynomial
- **alpha** (real): Parameter α in generalized Laguerre polynomial L_n^(α)(x)
- **Rmax** (real): Maximum value of scaled Lagrange-Laguerre mesh points (fm)
- **ctheta** (real): Rotation angle for complex scaling (degrees)
- **matgauss** (logical): Use Gauss-Legendre quadrature for Hamiltonian matrix elements
- **bgauss** (logical): Use Gauss-Legendre quadrature for inhomogeneous terms
- **numgauss** (integer): Number of Gauss-Legendre mesh points
- **rmaxgauss** (real): Maximum radius of Gauss-Legendre mesh (fm)
- **method** (integer): Calculation method
  - `1` = Linear Equation Method (works with nonlocal)
  - `2` = Green's Function Method (does NOT work with nonlocal)
- **backrot** (logical): Rotate basis backward (true) or rotate potential (false)
- **cwftype** (integer): Coulomb wave function library
  - `1` = COULCC (Fortran, **RECOMMENDED** - much faster)
  - `2` = cwfcomplex (C++ library)
- **thetah** (real): Angular step size for differential cross section output (degrees)
- **thetamax** (real): Maximum scattering angle for output (degrees)
- **readinpot** (logical): Read external potential from `pot.dat` file
  - **Note**: Must set `backrot=.true.` when using external potential

### System namelist
Reaction system properties:

- **zp, zt** (real): Charge number of projectile and target
- **massp, masst** (real): Mass number of projectile and target
- **namep, namet** (character): Name of projectile and target
- **jmin, jmax** (real): Minimum/maximum total angular momentum J
- **sp** (real): Spin of the projectile
- **elab** (real): Incident kinetic energy in lab frame (MeV)

### Pot namelist
Optical model potential parameters:

- **a1, a2** (real): Mass numbers for radius calculation
  - Default behavior (a1=0, a2=0): Uses a1=0, a2=masst
  - R₀ = r₀ × (a1^(1/3) + a2^(1/3))

- **vv, rv, av** (real): Depth (MeV), radius (fm), diffuseness (fm) of real volume term
- **wv, rw, aw** (real): Depth, radius, diffuseness of imaginary volume term
- **vs, rvs, avs** (real): Depth, radius, diffuseness of real surface term
- **ws, rws, aws** (real): Depth, radius, diffuseness of imaginary surface term
- **vsov, rsov, asov** (real): Depth, radius, diffuseness of real spin-orbit term
- **vsow, rsow, asow** (real): Depth, radius, diffuseness of imaginary spin-orbit term
- **rc** (real): Coulomb radius parameter (fm)

**Note**: The radius parameters use R = r × a13, where a13 is calculated from a1 and a2.

### Nonlocalpot namelist
Nonlocal potential settings:

- **nonlocal** (logical): Use Perey-Buck nonlocal form (Nuclear Physics, 32, 353-380)
- **nlbeta** (real): Nonlocal range parameter β (fm)

**Important**: Nonlocal potential requires `method=1` (Linear Equation Method)

---

## Output Description

COLOSS generates the following output files:

- **fort.1**: Local copy of the input file
- **fort.2**: Table of angular momentum channels considered
- **fort.10**: Scaled Lagrange-Laguerre mesh points and weights
- **fort.60**: S-matrices for different angular momentum channels
- **fort.61**: Nuclear scattering amplitudes for different channels
- **fort.67**: Angular distribution of differential cross section

The GUI automatically reads and visualizes data from fort.60, fort.61, and fort.67.

---

## Examples

### Example 1: ⁴⁰Ca(n,n)⁴⁰Ca at 20 MeV

#### Input file (test/40Ca_n_20MeV.in)

```fortran
&general
    nr=60  alpha=0 Rmax=40 ctheta=7
    matgauss=f bgauss=f method=1
    thetah=1.0 thetamax=180
    cwftype=1 backrot=f /

&system
    zp=0    massp=1   namep='n'
    zt=20   masst=40  namet='40Ca'
    jmin=0 jmax=10  elab=20 sp=0.5 /

&pot
    a1=0.00  a2=0.00
    vv=46.553   rv=1.185    av=0.672
    wv=1.777    rw=1.185    aw=0.672
    vs=0        rvs=0       avs=0
    ws=7.182    rws=1.288   aws=0.538
    vsov=5.343  rsov=0.996  asov=0.590
    vsow=-0.110 rsow=0.996  asow=0.590
    rc=1.698 /

&nonlocalpot
    nonlocal=f nlbeta=0.0 /
```

#### Running the calculation

```bash
cd test
../COLOSS < 40Ca_n_20MeV.in
```

#### Output (excerpt)

```
************************************************************
*                  COLOSS: Complex-scaled                  *
*          Optical and couLOmb Scattering Solver           *
************************************************************

Rotation Angle for Complex Scaling:  7.00 degrees
Rotation Operation: Applied to Potential

-------------------Reaction system-------------------
  Projectile:  n     (A:     1.00 Z:    0.00)
      Target:  40Ca  (A:    40.00 Z:   20.00)
Lab  Frame Energy:   20.00 MeV
C.M. Frame Energy:   19.51 MeV

Generating the channel index
There are   20 different channels
Total J Range:  0.0 <= J <=  10.0

-------------- Lagrange-Laguerre Mesh --------------
Laguerre Polynomial Order:  60
Laguerre Mesh Max Value: Scaled from   219.3181 to    40.0000
Scaling Factor:    0.1824 (fm)

------------- Optical Model Potential -------------
Optical Potential for A =  40.0 Z = 20.0 at  20.000 MeV rc =  1.698
 Vv     rvv    avv    Wv     rw     aw
46.553  1.185  0.672  1.777  1.185  0.672

 Vs     rvs    avs    Ws     rws    aws
 0.000  0.000  0.000  7.182  1.288  0.538

 Vsov   rsov   asov   vsow   rsow   asow
 5.343  0.996  0.590 -0.110  0.996  0.590

------------------------------------------------
Initializing Complex Coulomb Function
Sommerfeld Parameter:    0.0000
Generating Rotated Coulomb Wave Function on Laguerre Mesh
                   (CPU  time =  0.00025700  seconds)

Using Linear Equation Method to Solve
 L     S      J |  (Real(Smat), Imag(Smat))  |  Partial Wave Reac_Xsec (mb)
----------------------------------------------------------------------------
 0   0.5    0.5 |  (  0.349804,   0.223235)  |        28.5534
 1   0.5    0.5 |  (  0.508537,   0.041797)  |        25.5125
 1   0.5    1.5 |  (  0.479616,   0.195699)  |        50.4749
 2   0.5    1.5 |  (  0.321689,  -0.152658)  |        60.2393
 2   0.5    2.5 |  (  0.392012,   0.086108)  |        86.8096
 3   0.5    2.5 |  ( -0.032075,  -0.496351)  |        77.8789
 3   0.5    3.5 |  (  0.317075,  -0.310480)  |       110.8004
 4   0.5    3.5 |  ( -0.134338,  -0.189716)  |       130.5159
 4   0.5    4.5 |  (  0.171986,  -0.508775)  |       122.7206
 5   0.5    4.5 |  (  0.448912,   0.217764)  |       129.5308
 5   0.5    5.5 |  (  0.219593,   0.066785)  |       196.0549
 6   0.5    5.5 |  (  0.878907,   0.126747)  |        43.7629
 6   0.5    6.5 |  (  0.839528,   0.167094)  |        64.5332
 7   0.5    6.5 |  (  0.977509,   0.037757)  |        10.3944
 7   0.5    7.5 |  (  0.975107,   0.047042)  |        12.9563
 8   0.5    7.5 |  (  0.995748,   0.009952)  |         2.3144
 8   0.5    8.5 |  (  0.995614,   0.011765)  |         2.6742
 9   0.5    8.5 |  (  0.999183,   0.002546)  |         0.5053
 9   0.5    9.5 |  (  0.999178,   0.002913)  |         0.5636
10   0.5    9.5 |  (  0.999841,   0.000644)  |         0.1093
----------------------------------------------------------------------------
                   (CPU  time =    0.003733  seconds)

----------------- Files Created -----------------
 1: local copy of input.
 2: list of channel index.
10: scaled Laguerre mesh points and weights.
60: S matrix LSJ distribution.
61: scat amplitude LSJ distribution.
67: cross section angular distribution.
         Total CPU  time =    0.029168  seconds
```

---

## Additional Information

### Performance Tips

1. **Use COULCC** (`cwftype=1`): Much faster than C++ library for most calculations
2. **Optimize mesh size**: Start with moderate `nr` (20-60) and increase if needed
3. **Gauss quadrature**: Enable `matgauss=t` and `bgauss=t` for better accuracy with fewer basis functions

### Troubleshooting

- **Compilation errors**: Check LAPACK path in Makefile
- **GUI won't start**: Ensure conda environment is activated: `conda activate coloss_gui`
- **Missing potential file**: When `readinpot=t`, ensure `pot.dat` exists in working directory

### Citation

If you use COLOSS in your research, please cite our published paper:

**Paper**: Liu, J., Lei, J., & Ren, Z. (2025). COLOSS: Complex-scaled Optical and couLOmb Scattering Solver. *Computer Physics Communications*, 311, 109568.
**DOI**: [10.1016/j.cpc.2025.109568](https://doi.org/10.1016/j.cpc.2025.109568)
**URL**: [https://www.sciencedirect.com/science/article/abs/pii/S0010465525000712?via%3Dihub](https://www.sciencedirect.com/science/article/abs/pii/S0010465525000712?via%3Dihub)

#### BibTeX Entry

```bibtex
@article{LIU2025109568,
    title = {COLOSS: Complex-scaled Optical and couLOmb Scattering Solver},
    journal = {Computer Physics Communications},
    volume = {311},
    pages = {109568},
    year = {2025},
    issn = {0010-4655},
    doi = {https://doi.org/10.1016/j.cpc.2025.109568},
    url = {https://www.sciencedirect.com/science/article/pii/S0010465525000712},
    author = {Junzhe Liu and Jin Lei and Zhongzhou Ren},
    keywords = {Complex scaling, Scattering theory, Nuclear reaction, Optical potential},
    abstract = {We introduce COLOSS, a program designed to address the scattering problem
    using a bound-state technique known as complex scaling. In this method, the oscillatory
    boundary conditions of the wave function are transformed into exponentially decaying ones,
    accommodating the long-range Coulomb interaction. The program implements the general local
    optical potential and the Perey-Buck non-local optical potential, with all potential
    parameters included in a well-designed input format for ease of use. The design offers
    users direct access to compute S-matrices and cross-sections for scattering processes
    involving a projectile of any spin interacting with a spin-0 target. We provide thorough
    discussions on the precision of Lagrange functions and their benefits in evaluating matrix
    elements. Additionally, COLOSS incorporates two distinct rotation methods, making it
    adaptable to potentials without analytical expressions. Comparative results demonstrate
    that COLOSS achieves high accuracy when compared with the direct integration method,
    Numerov, underscoring its utility and effectiveness in scattering calculations.}
}
```

---

## License

GPLv3 - See the LICENSE file for details

## Contact

For questions, bug reports, or contributions, please contact (contact information to be added).
