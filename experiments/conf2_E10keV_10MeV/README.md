# Experiment: Configuration 2 — 10 keV to 10 MeV

## Overview

This is the **primary experiment** of the thesis: a high-resolution simulation
of the time-dependent Fokker-Planck equation for electrons accelerated in a
solar flare loop, spanning the energy range **10 keV – 10 MeV**.

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `emin` | 10 keV | Minimum electron energy |
| `emax` | 10 000 keV (10 MeV) | Maximum electron energy |
| `ne` | 120 | Energy grid points |
| `nmu` | 60 | Pitch-angle grid points per hemisphere |
| `ntau` | 120 | Spatial (column depth) grid points |
| `ntr` | 99 | Number of output time snapshots |
| `tmax` | 10.51 s | Total simulation time |
| `dt` | 2.78 × 10⁻⁴ s | Time step |
| `dy` | 1.67 × 10⁻⁵ | Dimensionless time step |
| `rhomax` | 10¹¹ cm⁻³ | Ambient electron density |
| `B₀` | 100 G | Magnetic field at looptop |
| `rm` | 2.0 | Mirror ratio B_footpoint / B_looptop |
| `δ` | 5 | Injection power-law index |
| `E₀` | 10 keV | Injection reference energy |
| Injection | Triangular pulse, 2 s duration | Time profile |

## Grid Dimensions

The full 4D distribution function `φ(E, μ, τ, t)` lives on a grid of size:

```
121 (energy) × 121 (pitch angle) × 121 (space) × 100 (time) ≈ 178 M values
```

## Output File

| File | Size | Description |
|------|------|-------------|
| `fkrplk.conf2_E10kev10Mev` | ~2.2 GB | Full ASCII output with all 100 time snapshots |

> **Note:** The output file is excluded from version control via `.gitignore`
> due to its size. To regenerate it, build and run the solver:
> ```bash
> cd src/fortran && make clean && make && ./fkrplk.x
> ```
> The executable writes to `fkrplk.conf2_E10kev10Mev` in its working directory.

## Analysis

Load and visualize this experiment's output with the Python module:

```python
import sys
sys.path.insert(0, "scripts/python")
from tdfp_reader import read_fkrplk_output
from tdfp_plots import plot_energy_spectrum, plot_pitch_angle

data = read_fkrplk_output("experiments/conf2_E10keV_10MeV/fkrplk.conf2_E10kev10Mev")
plot_energy_spectrum(data, save="energy.png")
```

## Reference Plots

| File | Description |
|------|-------------|
| `nombre.jpg` | Reference plot from original analysis |

## Physical Context

This configuration models a symmetric coronal loop with:
- Constant ambient density along the loop
- Parabolic magnetic field profile B(s) = B₀(1 + (rm − 1)s²/xc²)
- Coulomb collisions, magnetic mirroring, and pitch-angle diffusion
- Power-law electron injection at the looptop (τ = 0)
