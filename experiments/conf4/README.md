# Experiment: Configuration 4

## Overview

Configuration 4 is a **high-resolution parameter study** with a finer grid
in both energy and pitch-angle dimensions (100 × 100) compared to the original
configuration, and 120 output time snapshots for detailed temporal resolution.

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ne` | 100 | Energy grid points |
| `nmu` | 100 | Pitch-angle grid points per hemisphere |
| `ntau` | 120 | Spatial (column depth) grid points |
| `ntr` | 120 | Number of output time snapshots |
| `dt` | 1.77 × 10⁻⁴ s | Time step |
| `dy` | 1.06 × 10⁻⁵ | Dimensionless time step |
| `rhomax` | 10¹¹ cm⁻³ | Ambient electron density |

## Grid Dimensions

```
101 (energy) × 201 (pitch angle) × 121 (space) × 121 (time)
```

## Output Files

| File | Size | Description |
|------|------|-------------|
| `fkrplk.conf4` | ~31 MB | Full output with 121 time snapshots |

## Reference Plots

| File | Description |
|------|-------------|
| `Window0_conf4.jpg` | Energy spectrum f(E) |
| `Window1_conf4.jpg` | Pitch-angle distribution f(μ) |
| `Window2_conf4.jpg` | Spatial profile f(s) |
| `Window3_conf4.jpg` | Time evolution f(t) |

## Analysis

```python
import sys
sys.path.insert(0, "scripts/python")
from tdfp_reader import read_fkrplk_output
from tdfp_plots import plot_energy_spectrum

data = read_fkrplk_output("experiments/conf4/fkrplk.conf4")
plot_energy_spectrum(data)
```

## Notes

- The higher pitch-angle resolution (nmu=100 → 201 points) provides smoother
  angular distributions, important for studying anisotropy in the electron beam.
- This configuration has the finest temporal sampling (120 snapshots), useful
  for animation and detailed time-evolution studies.
