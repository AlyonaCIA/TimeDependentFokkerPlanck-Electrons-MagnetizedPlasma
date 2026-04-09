# Experiment: Original Configuration

## Overview

This is the **original / baseline configuration** used during early development
and testing of the Fokker-Planck solver. It uses a coarser spatial and
pitch-angle grid than the primary experiment, resulting in smaller output files
that are convenient for quick validation and debugging.

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ne` | 60 | Energy grid points |
| `nmu` | 30 | Pitch-angle grid points per hemisphere |
| `ntau` | 30 | Spatial (column depth) grid points |
| `ntr` | 8 | Number of output time snapshots |
| `tmax` | 10.507 s | Total simulation time |
| `dt` | 5.90 × 10⁻⁴ s | Time step |
| `dy` | 3.54 × 10⁻⁵ | Dimensionless time step |
| `rhomax` | 10¹¹ cm⁻³ | Ambient electron density |

## Grid Dimensions

```
61 (energy) × 61 (pitch angle) × 31 (space) × 9 (time)
```

## Output Files

| File | Size | Description |
|------|------|-------------|
| `fkrplk.conforg` | ~13 MB | Original configuration output |
| `fkrplk.test` | ~13 MB | Test run output (identical parameters) |

Both files share the same header parameters and grid dimensions.

## Reference Plots

Plots generated from this configuration (Window*.jpg):

| File | Description |
|------|-------------|
| `Window0.jpg` | Energy spectrum f(E) |
| `Window1.jpg` | Pitch-angle distribution f(μ) |
| `Window2.jpg` | Spatial profile f(s) |
| `Window3.jpg` | Time evolution f(t) |
| `Window0_conf4.jpg` – `Window3_conf4.jpg` | Same four views for Config 4 (cross-reference) |

## Analysis

```python
import sys
sys.path.insert(0, "scripts/python")
from tdfp_reader import read_fkrplk_output
from tdfp_plots import plot_energy_spectrum

data = read_fkrplk_output("experiments/conf_original/fkrplk.conforg")
plot_energy_spectrum(data)
```

## Notes

- This configuration runs much faster than the full-resolution conf2 experiment
  due to the 4× coarser grid in each dimension.
- Useful as a quick smoke-test after modifying the Fortran source.
