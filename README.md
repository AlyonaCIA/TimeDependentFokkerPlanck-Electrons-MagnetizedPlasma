<p align="center">
  <h1 align="center">TDFP — Time-Dependent Fokker-Planck Solver<br>for Electrons in Magnetized Plasma</h1>
  <p align="center"><strong>Numerical solution of the relativistic, time-dependent Fokker-Planck kinetic equation for energetic electrons propagating along converging magnetic loops in solar flares.</strong></p>
  <p align="center">Master's thesis in Mathematics — Pontificia Universidad Javeriana, 2021</p>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/fortran-77%20%7C%20F90-734F96?logo=fortran&logoColor=white" alt="Fortran">
  <img src="https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white" alt="Python 3.10+">
  <img src="https://img.shields.io/badge/uv-package%20manager-DE5FE9?logo=uv&logoColor=white" alt="uv">
  <img src="https://img.shields.io/badge/NumPy-1.24+-013243?logo=numpy&logoColor=white" alt="NumPy">
  <img src="https://img.shields.io/badge/Matplotlib-3.7+-11557C?logo=python&logoColor=white" alt="Matplotlib">
  <img src="https://img.shields.io/badge/PyVista-0.43+-4C72B0?logo=python&logoColor=white" alt="PyVista">
  <img src="https://img.shields.io/badge/Ruff-linter-D7FF64?logo=ruff&logoColor=black" alt="Ruff">
  <img src="https://img.shields.io/badge/license-Academic%20Free%20%7C%20Commercial-orange" alt="License">
</p>

---

## Table of Contents

- [Overview](#overview)
- [Astrophysical Context](#astrophysical-context)
- [Numerical Methods](#numerical-methods)
- [Project Architecture](#-project-architecture)
- [Repository Structure](#-repository-structure)
- [Getting Started](#-getting-started)
- [Running the Simulation](#running-the-simulation)
- [Post-Processing & Visualization](#-post-processing--visualization)
- [Output Figures](#-output-figures)
- [Configuration](#configuration)
- [Technology Stack](#-technology-stack)
- [References & Credits](#-references--credits)
- [License](#-license)
- [Authors & Contributors](#-authors--contributors)
- [Changelog](#changelog)

---

## Overview

This project implements a **time-dependent Fokker-Planck (TDFP) solver** that evolves the electron velocity distribution function $f(E, \mu, s, t)$ in a magnetized, fully-ionized hydrogen plasma representative of solar flare coronal loops.

The Fortran numerical core was **originally developed** by **R. Hamilton, E. Lu, and V. Petrosian** (Hamilton, Lu & Petrosian 1990, ApJ 354, 726) and subsequently modified by **J. McTiernan** and **L. Ofman**. This repository contains the codebase associated with a **2021 Master's thesis** in Mathematics at the Pontificia Universidad Javeriana, where the code was used to reproduce and extend the foundational results of **Leach & Petrosian (1981)** (LP81).

> **Note:** This thesis was originally completed in 2021. The repository is being uploaded to GitHub in 2025 for archival, reproducibility, and long-term maintenance purposes.

The solver tracks the phase-space evolution of non-thermal electron beams as they:

1. Are **injected** at the loop top with a power-law energy spectrum.
2. **Stream** along the magnetic field (spatial advection in column depth $\tau$).
3. Lose energy via **Coulomb collisions** with the background plasma (energy diffusion).
4. Undergo **pitch-angle scattering** (diffusion in $\mu$).
5. Experience **magnetic mirroring** in converging flux tubes ($\nabla B$ force).

---

## Astrophysical Context

During a **solar flare**, magnetic reconnection accelerates electrons to non-thermal energies (10 keV – 10 MeV). These electrons propagate along coronal magnetic loops toward the chromosphere, emitting **hard X-ray** (HXR) bremsstrahlung as they interact with the ambient plasma. The shape of the observed HXR spectrum encodes information about the electron transport, acceleration, and energy loss mechanisms.

The **Fokker-Planck kinetic equation** provides the fundamental description of this transport:

$$
\frac{\partial f}{\partial t} + \mu v \frac{\partial f}{\partial s} - \frac{(1-\mu^2)}{2}\frac{d \ln B}{d s} v \frac{\partial f}{\partial \mu} = \left(\frac{\partial f}{\partial t}\right)_{\mathrm{coll}}
$$

where:
| Symbol | Meaning |
|--------|---------|
| $f(E, \mu, s, t)$ | Electron distribution function |
| $E$ | Kinetic energy (keV) |
| $\mu = \cos\theta$ | Pitch-angle cosine |
| $s$ | Distance along the magnetic loop (cm) |
| $B(s)$ | Magnetic field strength along the loop |
| $v$ | Particle velocity ($v = \beta c$) |

The collision term on the right-hand side includes **energy loss**, **energy diffusion**, and **pitch-angle diffusion** due to Coulomb interactions, following the formalism of Rosenbluth, MacDonald & Judd (1957).

---

## Numerical Methods

The code employs **operator splitting** to decompose the full Fokker-Planck equation into four independent sub-operators, each solved sequentially per time step:

| Sub-operator | Physical Process | Numerical Scheme | Source File |
|---|---|---|---|
| **$\tau$-update** | Spatial advection along $B$ | Barton's Monotonic Transport (Upwind) | `tauup.f` |
| **$E$-update** | Coulomb energy losses | Barton's Monotonic Transport (Chang-Cooper-like) | `eup.f` |
| **$B$-update** | Magnetic mirroring ($\nabla B$ drift) | Barton's Monotonic Transport | `bcup.f` |
| **$\mu$-update** | Pitch-angle diffusion | Implicit Crank-Nicolson + Thomas Algorithm (TDMA) | `muup.f` + `cntridag.f` |

### Key Algorithmic Details

- **Barton's Monotonic Transport**: A flux-corrected scheme that prevents numerical oscillations near steep gradients while maintaining conservation. Applied to the advective terms in $\tau$, $E$, and $\mu$ (mirroring).
- **Crank-Nicolson**: A second-order implicit scheme for the diffusive pitch-angle term, unconditionally stable and solved via a tridiagonal matrix (Thomas algorithm / TDMA) implemented in `cntridag.f`.
- **Adaptive Time Stepping**: The time step `dy` is set to satisfy the CFL condition for all sub-operators (Courant, collisional, and mirroring stability limits). Near reporting times, the step is temporarily adjusted for exact output.

---

## 🏗 Project Architecture

The project follows a **dual-language architecture**:

```
┌─────────────────────────────────────────────┐
│              Fortran 77 Core Engine         │
│  (Heavy numerical computation, ~1000 LOC)   │
│                                             │
│  ┌─────────┐   ┌──────┐   ┌──────┐         │
│  │ init.f  │──>│main.f│──>│output│──> .dat  │
│  └─────────┘   └──┬───┘   └──────┘   file   │
│                   │                          │
│    ┌──────────┬───┴───┬──────────┐          │
│    ▼          ▼       ▼          ▼          │
│  tauup.f   eup.f   bcup.f    muup.f        │
│  (space)   (energy) (mirror) (pitch-angle)  │
│                                             │
│  Support: inject.f, bfield.f, dands.f,      │
│           cntridag.f, fkrplk.h              │
└─────────────┬───────────────────────────────┘
              │  writes ASCII output
              ▼
┌─────────────────────────────────────────────┐
│         Python Post-Processing              │
│  (Reads output, integrates, generates plots)│
│                                             │
│  tdfp_reader.py → plot_*.py → figures       │
│  Publication-quality astrophysical plots    │
└─────────────────────────────────────────────┘
```

### Fortran Source Files

| File | Purpose |
|------|---------|
| `main.f` | Program entry point. Controls time-stepping loop, I/O, and reporting. |
| `fkrplk.h` | Header with `PARAMETER` declarations for grid dimensions (`nemax`, `nmumax`, `ntaumax`) and `COMMON` blocks for shared arrays. |
| `init.f` | Initializes all grids ($E$, $\mu$, $\tau$), computes plasma parameters, density profiles, magnetic field, and sets the stable time step. |
| `inject.f` | Injects a non-thermal electron beam: power-law in $E$, isotropic in $\mu$, Gaussian/delta in space, triangular temporal pulse. |
| `tauup.f` | Spatial advection operator ($\partial/\partial\tau$) using monotonic upwind transport with symmetric/asymmetric loop boundary conditions. |
| `eup.f` | Energy update operator ($\partial/\partial E$) for Coulomb losses using Barton's monotonic transport. |
| `bcup.f` | Magnetic mirroring operator ($\nabla B$ force in $\mu$) using monotonic transport. |
| `muup.f` | Pitch-angle diffusion operator using Crank-Nicolson implicit scheme. |
| `cntridag.f` | Thomas Algorithm (TDMA) solver for the tridiagonal system from the Crank-Nicolson discretization. Adapted from *Numerical Recipes* (Press et al. 1986). |
| `bfield.f` | Returns $B(s)$ and $d\ln B/ds$ for a parabolic magnetic field: $B = B_0(1 + s^2/L_B^2)$. |
| `dands.f` | Returns electron number density $\rho(s)$ and spatial coordinate $s$ from column depth $\tau$. Currently assumes constant density $\rho_0 = 10^{11}$ cm$^{-3}$. |
| `main.f90` | Modern Fortran 2008 skeleton reimplementing the solver with modules (`fp_parameters`, `fp_grids`, `fp_physics`, `fp_solver`). |

### Current Simulation Configuration (Configuration 2)

- **Energy range**: 10 keV – 10 MeV (120 grid points, non-uniform spacing)
- **Pitch-angle range**: $\mu \in [-1, 1]$ (60 points per hemisphere, uniform)
- **Spatial range**: Column depth $\tau \in [0, 2 \times 10^{-3}]$ (120 points)
- **Background**: Fully ionized hydrogen, $n_e = 10^{11}$ cm$^{-3}$, $\ln\Lambda = 20$
- **Magnetic field**: $B_0 = 100$ G, mirror ratio $r_m = 2$, coronal length $L_C = 10^9$ cm
- **Injection**: Power-law index $\delta = 5$, triangular pulse (2 s), looptop injection

---

## 📁 Repository Structure

```text
TimeDependentFokkerPlanck-Electrons-MagnetizedPlasma/
│
├── README.md                       # This file
├── CHANGELOG.md                    # Project changelog
├── LICENSE                         # Dual license (Academic Free / Commercial)
├── pyproject.toml                  # Python project config (uv, ruff, mypy)
├── Makefile                        # Build & development workflow commands
├── requirements.txt                # Python dependencies (pip-compatible)
├── .gitignore                      # Fortran + Python + output ignore rules
│
├── src/fortran/                    # Fortran 77 solver source code
│   ├── fkrplk.h                    #   Parameter header (grid dimensions)
│   ├── main.f                      #   Entry point & time-stepping loop
│   ├── init.f                      #   Grid initialization & CFL time step
│   ├── inject.f                    #   Electron injection (power-law source)
│   ├── tauup.f                     #   Spatial advection operator
│   ├── eup.f                       #   Energy loss operator
│   ├── bcup.f                      #   Magnetic mirroring operator
│   ├── muup.f                      #   Pitch-angle diffusion (Crank-Nicolson)
│   ├── cntridag.f                  #   Thomas algorithm (TDMA solver)
│   ├── bfield.f                    #   Magnetic field model B(s)
│   ├── dands.f                     #   Density/distance model ρ(s), s(τ)
│   ├── main.f90                    #   Modern F2008 skeleton (4 modules)
│   └── Makefile                    #   gfortran build system
│
├── scripts/python/                 # Python post-processing & visualization
│   ├── tdfp_reader.py              #   Data reader for Fortran ASCII output
│   ├── tdfp_plots.py               #   Publication-quality plotting + CLI
│   ├── plot_2d_contours.py         #   LP88-style E–μ contour panels
│   ├── plot_3d_coronal_loop.py     #   3D coronal loop (PyVista)
│   ├── plot_density_evolution.py   #   Spatial density n_e(s,t)
│   └── plot_energy_spectrum.py     #   Energy spectrum cooling F(E,t)
│
├── experiments/                    # Simulation experiments (data + notebooks)
│   ├── conf2_E10keV_10MeV/         #   Primary: 120×60×120, 10 keV–10 MeV
│   ├── conf_original/              #   Baseline: 60×30×30, 8 snapshots
│   └── conf4/                      #   High-res: 100×100×120, 120 snapshots
│
├── output/figures/                 # Generated figures (tracked in git)
│   ├── README.md                   #   Figure documentation & physics
│   ├── 01_coronal_loop_3d.png      #   3D coronal loop geometry
│   ├── 02_lp88_contours_mock.png   #   LP88 contours (synthetic data)
│   ├── 03_lp88_contours_real.png   #   LP88 contours (simulation, t=0)
│   ├── 04_lp88_contours_real_t2.png#   LP88 contours (propagated beam)
│   ├── 05_lp88_contours_real_t1.png#   LP88 contours (early propagation)
│   ├── 06_density_evolution.png    #   Density n_e(s) time evolution
│   └── 07_energy_spectrum.png      #   Energy spectrum cooling F(E,t)
│
├── docs/                           # PDF documentation & references
├── archive/                        # Legacy IDL scripts & historical files
│   ├── legacy_idl/                 #   Original IDL plotting scripts
│   ├── scratch/                    #   WIP / scratch notebooks
│   └── old_structure/              #   Pre-reorganization directories
│
└── ci/                             # CI/CD scripts
    ├── run_local.sh                #   Local CI runner
    ├── requirements/               #   CI-specific requirements
    └── scripts/                    #   Formatting, linting, validation
```

---

## 🚀 Getting Started

### Prerequisites

- **Fortran compiler**: `gfortran` (GCC Fortran) or any F77/F90 compatible compiler
- **Python >= 3.10**
- **[uv](https://docs.astral.sh/uv/)** (modern Python package manager)
- **Git**
- **make** (build tool)

### Install uv

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with Homebrew
brew install uv
```

### Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/AlyonaCIA/TimeDependentFokkerPlanck-Electrons-MagnetizedPlasma.git
   cd TimeDependentFokkerPlanck-Electrons-MagnetizedPlasma
   ```

2. **Install Fortran compiler** (if not already available):

   ```bash
   # macOS
   brew install gcc

   # Ubuntu/Debian
   # sudo apt-get install -y gfortran
   ```

3. **Install Python dependencies** (uv creates the virtualenv automatically):

   ```bash
   uv sync
   ```

4. **Install with development tools** (Ruff, mypy, pre-commit):

   ```bash
   make install-dev
   # or: uv sync --extra dev
   ```

5. **Compile the Fortran solver:**

   ```bash
   make build
   # or: cd src/fortran && make clean && make
   ```

---

## Running the Simulation

```bash
make run
# or: cd src/fortran && ./fkrplk.x
```

The program writes an ASCII output file containing:

1. **Header**: Time step `dt`, max time `tmax`, grid dimensions.
2. **Grids**: Energy (keV), pitch-angle cosine $\mu$, column depth $\tau$, spatial position $s$, magnetic field $B$.
3. **Time snapshots**: The full 3D distribution $\phi(E, \mu, \tau)$ at each reporting time $t_r$.

---

## 📊 Post-Processing & Visualization

### Python Modules

Reusable Python post-processing tools are in `scripts/python/`:

| Script | Purpose |
|--------|---------|
| `tdfp_reader.py` | Reads Fortran ASCII output into structured `FPOutput` dataclass |
| `tdfp_plots.py` | Publication-quality plotting functions + CLI |
| `plot_2d_contours.py` | LP88-style $\log_{10}f(E,\mu)$ contour panels |
| `plot_3d_coronal_loop.py` | 3D volumetric coronal loop (PyVista + VTK) |
| `plot_density_evolution.py` | Spatial density $n_e(s,t) = \int f\,d\mu\,dE$ |
| `plot_energy_spectrum.py` | Energy spectrum $F(E,t) = \int f\,d\mu\,ds$ |

### Quick Start

```python
import sys
sys.path.insert(0, "scripts/python")
from tdfp_reader import read_fkrplk_output

data = read_fkrplk_output("experiments/conf_original/fkrplk.test")
# data.phi.shape → (n_time, n_space, n_pitch, n_energy)
```

### Command-Line Plotting

```bash
# LP88 contour panels (real data, propagated beam)
uv run python scripts/python/plot_2d_contours.py \
    --data experiments/conf_original/fkrplk.test --time-index 2 \
    --save output/figures/04_lp88_contours_real_t2.png

# 3D coronal loop
uv run python scripts/python/plot_3d_coronal_loop.py \
    --save output/figures/01_coronal_loop_3d.png --resolution 2560x1440

# Density evolution
uv run python scripts/python/plot_density_evolution.py \
    --data experiments/conf_original/fkrplk.test \
    --save output/figures/06_density_evolution.png

# Energy spectrum cooling
uv run python scripts/python/plot_energy_spectrum.py \
    --data experiments/conf_original/fkrplk.test \
    --save output/figures/07_energy_spectrum.png
```

### Jupyter Notebooks

Analysis notebooks are located within each experiment directory:

| Notebook | Experiment | Description |
|----------|------------|-------------|
| `experiments/conf2_E10keV_10MeV/plots_final_conf2_v0.ipynb` | Config 2 | Full analysis: energy, pitch-angle, spatial, and time plots |
| `experiments/conf_original/test_plot.ipynb` | Original | Basic test plotting |
| `experiments/conf_original/pruebaa1confiorg.ipynb` | Original | Comprehensive analysis |
| `experiments/conf4/test_plot_conf4.ipynb` | Config 4 | High-resolution run analysis |

---

## 🖼 Output Figures

All figures are in `output/figures/` with full physics descriptions in [`output/figures/README.md`](output/figures/README.md).

| # | Figure | Description |
|---|--------|-------------|
| 01 | 3D Coronal Loop | Volumetric semi-circular loop coloured by $\log_{10} n_e$ |
| 02 | LP88 Contours (Mock) | Synthetic $f(E,\mu)$ validating the LP88 contour style |
| 03 | LP88 Contours (Real, t=0) | Initial injection — isotropic power-law at $s=0$ |
| 04 | LP88 Contours (Real, t=2) | Beam propagation — anisotropy develops with depth |
| 05 | LP88 Contours (Real, t=1) | Early propagation — temporal convergence check |
| 06 | Density Evolution | $n_e(s,t)$ — injection → propagation → thermalisation |
| 07 | Energy Spectrum | $F(E,t)$ — power-law injection → collisional cooling |

---

## Configuration

The simulation parameters are set directly in the Fortran source files:

| Parameter | File | Variable | Default |
|-----------|------|----------|---------|
| Energy range | `init.f` | `emin`, `emax` | 10 keV, 10 MeV |
| Energy grid points | `init.f` / `fkrplk.h` | `ne` / `nemax` | 120 |
| Pitch-angle grid points | `init.f` / `fkrplk.h` | `nmu` / `nmumax` | 60 |
| Spatial grid points | `init.f` / `fkrplk.h` | `ntau` / `ntaumax` | 120 |
| Background density | `dands.f` | `rho0` | $10^{11}$ cm$^{-3}$ |
| Coulomb logarithm | `init.f` | `aloglmda` | 20 |
| Magnetic field at loop top | `bfield.f` | `b0` | 100 G |
| Mirror ratio | `bfield.f` | `rm` | 2.0 |
| Coronal loop length | `bfield.f` | `xc` | $10^9$ cm |
| Injection spectral index | `inject.f` | `delta` | 5.0 |
| Injection duration | `inject.f` | (conditional) | 2 s triangular |
| Max simulation time | `main.f` | `tmax` | 10.51 s |

### Available Make Commands

```bash
make help          # Show all available commands
make build         # Compile Fortran solver
make run           # Build and execute simulation
make clean         # Remove build artifacts
make install       # Install Python dependencies (uv sync)
make install-dev   # Install with dev tools (ruff, mypy)
make lint          # Run Ruff linter
make format        # Auto-fix formatting
make check         # lint + type-check
make plot          # Generate standard plots
```

---

## 🛠 Technology Stack

| Category | Technologies |
|---|---|
| **Numerical Core** | Fortran 77 (original solver), Fortran 2008 (modern skeleton) |
| **Compiler** | gfortran (GCC 15.2+) |
| **Python** | Python 3.10+ |
| **Package Manager** | [uv](https://docs.astral.sh/uv/) |
| **Numerical Computing** | NumPy ≥ 1.24 |
| **Visualization** | Matplotlib ≥ 3.7, Seaborn ≥ 0.13 |
| **3D Visualization** | PyVista ≥ 0.43, VTK |
| **Linting & Formatting** | [Ruff](https://docs.astral.sh/ruff/) |
| **Type Checking** | mypy |
| **Build System** | GNU Make |

---

## 📚 References & Credits

### Original Code Authors

The Fortran numerical solver was **originally developed** by the following researchers. This thesis builds upon their work, and all credit for the core algorithm is due to them:

1. **Hamilton, R. J., Lu, E. T. & Petrosian, V.** (1990). *"Numerical Solution of the Time-Dependent Kinetic Equation for Electrons in Magnetized Plasma."* The Astrophysical Journal, 354, 726. [ADS](https://ui.adsabs.harvard.edu/abs/1990ApJ...354..726H)

2. **McTiernan, J. M. & Petrosian, V.** (1990). *"The Behavior of Beams of Relativistic Nonthermal Electrons Under the Influence of Collisions and Synchrotron Losses."* The Astrophysical Journal, 359, 524. [ADS](https://ui.adsabs.harvard.edu/abs/1990ApJ...359..524M)

### Key Physics References

3. **Leach, J. & Petrosian, V.** (1981). *"Impulsive Phase of Solar Flares. I. Characteristics of High Energy Electrons."* The Astrophysical Journal, 251, 781. [ADS](https://ui.adsabs.harvard.edu/abs/1981ApJ...251..781L)

4. **Leach, J. & Petrosian, V.** (1983). *"Impulsive Phase of Solar Flares. II. Characteristics of the Hard X-Rays."* The Astrophysical Journal, 269, 715. [ADS](https://ui.adsabs.harvard.edu/abs/1983ApJ...269..715L)

5. **Petrosian, V.** (1985). *"Directivity of Bremsstrahlung Radiation from Relativistic Beams and the Gamma Rays from Solar Flares."* The Astrophysical Journal, 299, 987.

6. **Rosenbluth, M. N., MacDonald, W. M. & Judd, D. L.** (1957). *"Fokker-Planck Equation for an Inverse-Square Force."* Physical Review, 107, 1.

7. **Press, W. H., Flannery, B. P., Teukolsky, S. A. & Vetterling, W. T.** (1986). *Numerical Recipes.* Cambridge University Press. (Thomas Algorithm implementation)

### IDL Visualization

8. **Holman, G.** (2001). Original IDL visualization script `tdfp_plot.pro` (archived in `archive/legacy_idl/`).

---

## 📄 License

This project is distributed under a **dual license** model:

| Use Case | License | Cost |
|---|---|---|
| Academic, research, educational, non-profit | **Academic & Non-Commercial License** | Free |
| Commercial, enterprise, revenue-generating | **Commercial License** | Paid (contact author) |

See [LICENSE](LICENSE) for full terms.

**Questions?** Contact: alenacivanovaa@gmail.com

---

## Authors & Contributors

<div align="center">

### Principal Investigator

**Alyona Carolina Ivanova Araujo**

Master's Thesis in Mathematics — Pontificia Universidad Javeriana, 2021

**Email:** alenacivanovaa@gmail.com
**GitHub:** [@AlyonaCIA](https://github.com/AlyonaCIA)

---

### Thesis Director

**Jorge Andrés Plazas Vargas**
Profesor Asistente — Departamento de Matemáticas, Pontificia Universidad Javeriana

### Thesis Co-Director

**Juan Carlos Martínez Oliveros**
Assistant Research Physicist — Space Sciences Laboratory, University of California, Berkeley

### Collaborator

**Juan Camilo Guevara Gómez**
PhD — University of Oslo

</div>

---

## Acknowledgments

This thesis would not have been possible without the foundational work of the original code authors — R. Hamilton, E. Lu, V. Petrosian, J. McTiernan, and L. Ofman — whose Fortran solver forms the numerical core of this project.

Special thanks to:
- **V. Petrosian** and collaborators for the Fokker-Planck formalism and original codebase
- **G. Holman** for the IDL visualization tools
- **Pontificia Universidad Javeriana**, Departamento de Matemáticas
- **Space Sciences Laboratory**, UC Berkeley
- The open-source scientific Python community (NumPy, Matplotlib, PyVista)

---

<div align="center">

**If you use this code in your research, please cite the references listed above.**

Thesis completed in 2021 — uploaded to GitHub in 2025 for archival and maintenance.

</div>

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full release history.