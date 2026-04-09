"""I/O module for the Time-Dependent Fokker-Planck solver output.

Reads the ASCII output produced by the Fortran ``fkrplk`` program and returns
structured NumPy arrays of the electron distribution function and all grid
coordinates.

File format (following the original IDL reader ``tdfp_plot.pro``, G. Holman 2001)::

    Line 1:  dt, tmax, dy, rhomax          (4 floats)
    Line 2:  newr, nmu, ntau, ntr          (4 integers)

    Unrolled 1-D block:
        ekev[newr+1]       — Energy grid in keV
        amu[2*nmu+1]       — Pitch-angle cosine grid
        dtau[ntau]         — Column-depth step sizes
        tau[ntau]          — Column-depth grid (0 prepended)
        soft[ntau]         — Spatial position grid (0 prepended)
        bm[ntau+1]         — Magnetic field along loop

    Per time snapshot (ntr+1 total):
        tr, iteration      — Report time and iteration index
        phi[ntau+1][2*nmu+1][newr+1]  — Distribution function block

Example
-------
>>> from tdfp_reader import read_fkrplk_output
>>> data = read_fkrplk_output("experiments/conf2_E10keV_10MeV/fkrplk.conf2")
>>> data.phi.shape  # (time, space, pitch, energy)
(100, 121, 121, 121)
"""

from __future__ import annotations

# -- Standard library --
from dataclasses import dataclass
from pathlib import Path

# -- Third-party --
import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class FPOutput:
    """Container for all arrays read from a ``fkrplk`` output file.

    Attributes
    ----------
    dt : float
        Simulation time-step in seconds.
    tmax : float
        Maximum simulation time in seconds.
    dy : float
        Dimensionless time-step.
    rhomax : float
        Background electron density in cm⁻³.
    newr : int
        Number of energy grid points written to the file.
    nmu : int
        Number of pitch-angle grid points per hemisphere.
    ntau : int
        Number of spatial (column-depth) grid points.
    ntr : int
        Number of output time snapshots minus one.
    ekev : NDArray[np.float64]
        Energy grid in keV, shape ``(newr + 1,)``.
    eta : NDArray[np.float64]
        Dimensionless kinetic energy ``E / mc²``, shape ``(newr + 1,)``.
    beta : NDArray[np.float64]
        Velocity ``v / c``, shape ``(newr + 1,)``.
    amu : NDArray[np.float64]
        Cosine of pitch angle ``μ = cos(θ)``, shape ``(2 * nmu + 1,)``.
    pad : NDArray[np.float64]
        Pitch angle in degrees, shape ``(2 * nmu + 1,)``.
    dtau_arr : NDArray[np.float64]
        Column-depth step sizes, shape ``(ntau,)``.
    tau : NDArray[np.float64]
        Column-depth grid with 0 prepended, shape ``(ntau + 1,)``.
    soft : NDArray[np.float64]
        Spatial position along the loop in cm, shape ``(ntau + 1,)``.
    bm : NDArray[np.float64]
        Magnetic field strength along the loop, shape ``(ntau + 1,)``.
    tr : NDArray[np.float64]
        Report times in seconds, shape ``(ntr + 1,)``.
    iterations : NDArray[np.float64]
        Iteration index at each report time.
    phi : NDArray[np.float64]
        4-D distribution function with axes
        ``[time, space, pitch_angle, energy]``.
    """

    # Simulation metadata
    dt: float
    tmax: float
    dy: float
    rhomax: float
    newr: int
    nmu: int
    ntau: int
    ntr: int

    # Grids
    ekev: NDArray[np.float64]       # Energy (keV), shape (newr+1,)
    eta: NDArray[np.float64]        # Dimensionless energy E/mc^2
    beta: NDArray[np.float64]       # v/c
    amu: NDArray[np.float64]        # cos(pitch angle), shape (2*nmu+1,)
    pad: NDArray[np.float64]        # Pitch angle (degrees)
    dtau_arr: NDArray[np.float64]   # Column depth steps, shape (ntau,)
    tau: NDArray[np.float64]        # Column depth, shape (ntau+1,) with 0 prepended
    soft: NDArray[np.float64]       # Spatial position (cm), shape (ntau+1,)
    bm: NDArray[np.float64]         # Magnetic field, shape (ntau+1,)
    tr: NDArray[np.float64]         # Report times (s), shape (ntr+1,)
    iterations: NDArray[np.float64] # Iteration indices at each report time

    # Distribution function: phi[time, space, pitch_angle, energy]
    phi: NDArray[np.float64]


def read_fkrplk_output(filepath: str | Path) -> FPOutput:
    """Read a ``fkrplk`` ASCII output file into an :class:`FPOutput` dataclass.

    Handles truncated files gracefully — if the file contains fewer time
    snapshots than declared in the header, only the available snapshots
    are returned.

    Parameters
    ----------
    filepath : str | Path
        Path to the output file (e.g. ``fkrplk.conf2``).

    Returns
    -------
    FPOutput
        Structured container with all grids and the 4-D distribution.

    Raises
    ------
    FileNotFoundError
        If *filepath* does not exist.
    """
    filepath = Path(filepath)

    with open(filepath, "r") as fh:
        lines = fh.readlines()

    # --- Line 1: simulation parameters ---
    header = [float(x) for x in lines[0].split() if x.strip()]
    dt, tmax, dy, rhomax = header[0], header[1], header[2], header[3]

    # --- Line 2: grid dimensions ---
    dims = [int(x) for x in lines[1].split() if x.strip()]
    newr, nmu, ntau, ntr = dims[0], dims[1], dims[2], dims[3]

    # --- Unroll all remaining lines into a flat array ---
    values = np.array(
        [float(x) for line in lines[2:] for x in line.split() if x.strip()]
    )

    idx = 0  # running index into the flat array

    # Energy grid
    ekev = values[idx : idx + newr + 1].copy()
    idx += newr + 1
    eta = ekev / 510.99
    beta = np.sqrt(1.0 - 1.0 / (eta + 1.0) ** 2)

    # Pitch-angle grid
    n_amu = 2 * nmu + 1
    amu = values[idx : idx + n_amu].copy()
    idx += n_amu
    pad = 180.0 * np.arccos(amu) / np.pi

    # Column depth steps
    dtau_arr = values[idx : idx + ntau].copy()
    idx += ntau

    # Column depth (prepend 0)
    tau = np.insert(values[idx : idx + ntau], 0, 0.0)
    idx += ntau

    # Spatial position (prepend 0)
    soft = np.insert(values[idx : idx + ntau], 0, 0.0)
    idx += ntau

    # Magnetic field
    bm = values[idx : idx + ntau + 1].copy()
    idx += ntau + 1

    # --- Time snapshots ---
    chunk_size = (newr + 1) * n_amu * (ntau + 1)
    n_times = ntr + 1

    tr = np.zeros(n_times)
    iterations = np.zeros(n_times)
    phi = np.zeros((n_times, ntau + 1, n_amu, newr + 1))

    for itr in range(n_times):
        # Stop gracefully if the file is truncated
        if idx + 2 + chunk_size > len(values):
            n_times = itr
            tr = tr[:n_times]
            iterations = iterations[:n_times]
            phi = phi[:n_times]
            ntr = n_times - 1
            break
        tr[itr] = values[idx]
        iterations[itr] = values[idx + 1]
        idx += 2
        block = values[idx : idx + chunk_size]
        # Fortran writes: do m (space), do l (pitch), write (energy)
        # → fastest index = energy, then pitch, then space
        phi[itr] = block.reshape((ntau + 1, n_amu, newr + 1))
        idx += chunk_size

    return FPOutput(
        dt=dt, tmax=tmax, dy=dy, rhomax=rhomax,
        newr=newr, nmu=nmu, ntau=ntau, ntr=ntr,
        ekev=ekev, eta=eta, beta=beta,
        amu=amu, pad=pad, dtau_arr=dtau_arr,
        tau=tau, soft=soft, bm=bm,
        tr=tr, iterations=iterations,
        phi=phi,
    )
