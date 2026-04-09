<p align="center">
  <h1 align="center">Time-Dependent Fokker-Planck Solver<br>for Electrons in Magnetized Plasma</h1>
  <p align="center">
    Numerical solution of the relativistic, time-dependent Fokker-Planck kinetic equation<br>
    for energetic electrons propagating along converging magnetic loops in solar flares.
  </p>
</p>

---

## Table of Contents

- [Overview](#overview)
- [Astrophysical Context](#astrophysical-context)
- [Numerical Methods](#numerical-methods)
- [Project Architecture](#project-architecture)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Running the Simulation](#running-the-simulation)
- [Post-Processing & Visualization](#post-processing--visualization)
- [Configuration](#configuration)
- [References](#references)
- [License](#license)

---

## Overview

This project implements a **time-dependent Fokker-Planck (TDFP) solver** that evolves the electron velocity distribution function $f(E, \mu, s, t)$ in a magnetized, fully-ionized hydrogen plasma representative of solar flare coronal loops.

The code was originally developed by **R. Hamilton, E. Lu, and V. Petrosian** (Hamilton, Lu & Petrosian 1990, ApJ 354, 726) and later modified by **J. McTiernan** and **L. Ofman**. This repository is the modernized codebase associated with a 2021 Master's thesis on the numerical solution of this equation, aiming to reproduce and extend the foundational results of **Leach & Petrosian (1981)** (hereafter **LP81**).

The solver tracks the phase-space evolution of non-thermal electron beams as they:

1. Are **injected** at the loop top with a power-law energy spectrum.
2. **Stream** along the magnetic field (spatial advection in column depth $\tau$).
3. Lose energy via **Coulomb collisions** with the background plasma (energy diffusion).
4. Undergo **pitch-angle scattering** (diffusion in $\mu$).
5. Experience **magnetic mirroring** in converging flux tubes ($\nabla B$ force).

---

## Astrophysical Context

During a **solar flare**, magnetic reconnection accelerates electrons to non-thermal energies (10 keV -- 10 MeV). These electrons propagate along coronal magnetic loops toward the chromosphere, emitting **hard X-ray** (HXR) bremsstrahlung as they interact with the ambient plasma. The shape of the observed HXR spectrum encodes information about the electron transport, acceleration, and energy loss mechanisms.

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

## Project Architecture

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
              │  writes binary/ASCII output
              ▼
┌─────────────────────────────────────────────┐
│         Python Post-Processing              │
│  (Reads output, generates plots)            │
│                                             │
│  Jupyter notebooks + matplotlib/numpy       │
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

### Current Simulation Configuration (Configuration 2)

- **Energy range**: 10 keV -- 10 MeV (120 grid points, non-uniform spacing)
- **Pitch-angle range**: $\mu \in [-1, 1]$ (60 points per hemisphere, uniform)
- **Spatial range**: Column depth $\tau \in [0, 2 \times 10^{-3}]$ (120 points)
- **Background**: Fully ionized hydrogen, $n_e = 10^{11}$ cm$^{-3}$, $\ln\Lambda = 20$
- **Magnetic field**: $B_0 = 100$ G, mirror ratio $r_m = 2$, coronal length $L_C = 10^9$ cm
- **Injection**: Power-law index $\delta = 5$, triangular pulse (2 s), looptop injection

---

## Repository Structure

```
TimeDependentFokkerPlanck-Electrons-MagnetizedPlasma/
│
├── README.md                       # This file
├── CHANGELOG.md                    # Project changelog
├── requirements.txt                # Python dependencies
├── .gitignore                      # Git ignore (Fortran + Python + output)
│
├── src/
│   └── fortran/                    # Fortran 77 solver source code
│       ├── fkrplk.h                # Parameter header (grid dimensions)
│       ├── main.f                  # Entry point & time-stepping loop
│       ├── init.f                  # Grid initialization & CFL time step
│       ├── inject.f                # Electron injection (power-law source term)
│       ├── tauup.f                 # Spatial advection operator
│       ├── eup.f                   # Energy loss operator
│       ├── bcup.f                  # Magnetic mirroring operator
│       ├── muup.f                  # Pitch-angle diffusion (Crank-Nicolson)
│       ├── cntridag.f              # Thomas algorithm (TDMA solver)
│       ├── bfield.f                # Magnetic field model B(s)
│       ├── dands.f                 # Density/distance model ρ(s), s(τ)
│       └── Makefile                # gfortran build system
│
├── scripts/
│   └── python/                     # Python post-processing modules
│       ├── tdfp_reader.py          # Data reader for Fortran ASCII output
│       └── tdfp_plots.py           # Publication-quality plotting functions
│
├── experiments/                    # Simulation experiments (data + notebooks)
│   ├── conf2_E10keV_10MeV/        # Primary: 120×60×120 grid, 10 keV–10 MeV
│   ├── conf_original/             # Baseline: 60×30×30 grid, 8 snapshots
│   └── conf4/                     # High-res: 100×100×120 grid, 120 snapshots
│
├── output/                         # Simulation output directory (gitignored)
├── docs/                           # PDF documentation & references
│
└── archive/                        # Legacy and historical files
    ├── legacy_idl/                 # Original IDL scripts & old Makefiles
    ├── scratch/                    # WIP files, scratch notebooks
    └── old_structure/              # Pre-reorganization directories (for reference)
```

---

## Getting Started

### Prerequisites

- **Fortran compiler**: `gfortran` (GCC Fortran) or any F77/F90 compatible compiler
- **Python >= 3.9** with the packages listed in `requirements.txt`
- **make** (build tool)

### Installation

```bash
# Clone the repository
git clone https://github.com/<your-username>/TimeDependentFokkerPlanck-Electrons-MagnetizedPlasma.git
cd TimeDependentFokkerPlanck-Electrons-MagnetizedPlasma

# Install Fortran compiler (macOS)
brew install gcc

# Install Fortran compiler (Ubuntu/Debian)
# sudo apt-get install -y gfortran

# Install Python dependencies
pip install -r requirements.txt
```

### Compiling the Fortran Code

```bash
cd src/fortran/
make clean && make
```

This produces the executable `fkrplk.x` in `src/fortran/`.

---

## Running the Simulation

```bash
cd src/fortran/
./fkrplk.x
```

The program writes the output file `fkrplk.conf2_E10kev10Mev` containing:

1. **Header**: Time step `dt`, max time `tmax`, grid dimensions.
2. **Grids**: Energy (keV), pitch-angle cosine $\mu$, column depth $\tau$, spatial position $s$, magnetic field $B$.
3. **Time snapshots**: The full 3D distribution $\phi(E, \mu, \tau)$ at each reporting time $t_r$.

---

## Post-Processing & Visualization

### Python Modules

Reusable Python post-processing tools are in `scripts/python/`:

- **`tdfp_reader.py`** — Reads the Fortran ASCII output into a structured `FPOutput` dataclass with all grids and the 4D distribution function `φ[time, space, pitch, energy]`.
- **`tdfp_plots.py`** — Publication-quality plotting functions for the four fundamental 1D projections (energy, pitch-angle, spatial, time) plus 2D contour panels and 3D surfaces.

```python
import sys
sys.path.insert(0, "scripts/python")
from tdfp_reader import read_fkrplk_output
from tdfp_plots import plot_energy_spectrum, plot_2d_panels

data = read_fkrplk_output("experiments/conf2_E10keV_10MeV/fkrplk.conf2_E10kev10Mev")
plot_energy_spectrum(data, save="energy.png")
plot_2d_panels(data, save="contours.png")
```

### Command-Line Batch Plotting

Generate all standard figures for any output file:

```bash
python scripts/python/tdfp_plots.py experiments/conf4/fkrplk.conf4 -o output/plots/
```

### Jupyter Notebooks

Analysis notebooks are located within each experiment directory:

| Notebook | Experiment | Description |
|----------|------------|-------------|
| `experiments/conf2_E10keV_10MeV/plots_final_conf2_v0.ipynb` | Config 2 | Full analysis with energy, pitch-angle, spatial, and time plots |
| `experiments/conf_original/test_plot.ipynb` | Original | Basic test plotting |
| `experiments/conf_original/pruebaa1confiorg.ipynb` | Original | Comprehensive analysis |
| `experiments/conf4/test_plot_conf4.ipynb` | Config 4 | Analysis for high-resolution run |

### Legacy IDL Script

The original IDL visualization script `tdfp_plot.pro` (G. Holman, 2001) is archived in `archive/legacy_idl/`.

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

---

## References

1. **Leach, J. & Petrosian, V.** (1981). *"Impulsive Phase of Solar Flares. I. Characteristics of High Energy Electrons."* The Astrophysical Journal, 251, 781. [ADS](https://ui.adsabs.harvard.edu/abs/1981ApJ...251..781L)

2. **Hamilton, R. J., Lu, E. T. & Petrosian, V.** (1990). *"Numerical Solution of the Time-Dependent Kinetic Equation for Electrons in Magnetized Plasma."* The Astrophysical Journal, 354, 726. [ADS](https://ui.adsabs.harvard.edu/abs/1990ApJ...354..726H)

3. **Petrosian, V.** (1985). *"Directivity of Bremsstrahlung Radiation from Relativistic Beams and the Gamma Rays from Solar Flares."* The Astrophysical Journal, 299, 987.

4. **Press, W. H., Flannery, B. P., Teukolsky, S. A. & Vetterling, W. T.** (1986). *Numerical Recipes.* Cambridge University Press. (Thomas Algorithm implementation)

5. **Rosenbluth, M. N., MacDonald, W. M. & Judd, D. L.** (1957). *"Fokker-Planck Equation for an Inverse-Square Force."* Physical Review, 107, 1.

---

## License

This project is provided for academic and research purposes. Please cite the references above if you use this code in publications.