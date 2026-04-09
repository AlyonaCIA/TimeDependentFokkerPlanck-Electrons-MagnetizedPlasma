"""Time-evolution of spatially-resolved electron density.

Computes the macroscopic electron density by integrating the distribution
function over energy and pitch-angle:

.. math::

    n_e(s, t) = \\int_{E_{\\min}}^{E_{\\max}} \\int_{-1}^{1}
                f(E, \\mu, s, t)\\, d\\mu\\, dE

and plots :math:`n_e(s)` as 1-D lines for several time snapshots, coloured
with a sequential colourmap (viridis) from early (dark) to late (yellow).

Usage
-----
Real Fortran output::

    uv run python scripts/python/plot_density_evolution.py \\
        --data experiments/conf_original/fkrplk.test

Save to file::

    uv run python scripts/python/plot_density_evolution.py \\
        --data experiments/conf_original/fkrplk.test \\
        --save output/figures/06_density_evolution.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure


# ---------------------------------------------------------------------------
# Physics: integrate f(E, mu) -> n_e(s)
# ---------------------------------------------------------------------------

def compute_density(
    phi_t: np.ndarray,
    ekev: np.ndarray,
    amu: np.ndarray,
) -> np.ndarray:
    """Integrate f(E, mu, s) over energy and pitch-angle cosine.

    Parameters
    ----------
    phi_t : ndarray, shape (n_space, n_pitch, n_energy)
        Distribution function at a single time snapshot.
    ekev : ndarray, shape (n_energy,)
        Energy grid in keV.
    amu : ndarray, shape (n_pitch,)
        Pitch-angle cosine grid (μ = cos θ), ordered −1 … +1.

    Returns
    -------
    ne : ndarray, shape (n_space,)
        Electron number density at each spatial point.
    """
    # Replace NaN/Inf with 0 for integration
    phi_clean = np.where(np.isfinite(phi_t), phi_t, 0.0)

    # Integrate over energy (axis=-1), then over pitch-angle (axis=-1)
    integrated_E = np.trapezoid(phi_clean, x=ekev, axis=-1)  # (space, pitch)
    ne = np.trapezoid(integrated_E, x=amu, axis=-1)           # (space,)
    return ne


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_density_evolution(
    soft: np.ndarray,
    ne_all: np.ndarray,
    time_labels: list[str],
    *,
    cmap: str = "viridis",
    save: str | None = None,
) -> Figure:
    """Plot n_e(s) lines for multiple time steps.

    Parameters
    ----------
    soft : ndarray, shape (n_space,)
        Spatial position along the loop in cm.
    ne_all : ndarray, shape (n_times, n_space)
        Electron density at each time snapshot.
    time_labels : list of str
        Legend labels for each time snapshot.
    cmap : str
        Colourmap name.
    save : str or None
        Save figure to this path if given.

    Returns
    -------
    Figure
    """
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

    # Spatial axis in Mm for readability
    s_mm = soft / 1e8

    n_lines = len(time_labels)
    cm = plt.get_cmap(cmap)
    colours = [cm(i / max(n_lines - 1, 1)) for i in range(n_lines)]

    for i, (ne, label) in enumerate(zip(ne_all, time_labels)):
        ax.plot(s_mm, ne, color=colours[i], lw=2.2, label=label)

    ax.set_xlabel("Position along loop  $s$  (Mm)", fontsize=13)
    ax.set_ylabel(r"Electron density  $n_e(s)$  (arb. units)", fontsize=13)
    ax.set_title(
        r"Time-Evolution of Electron Density $n_e(s,\,t)$",
        fontsize=15, fontweight="bold",
    )

    ax.legend(fontsize=11, title="Time step", title_fontsize=12)
    ax.tick_params(labelsize=11)
    ax.grid(True, ls=":", lw=0.5, alpha=0.5)

    # Use log scale for y if range spans > 2 orders of magnitude
    finite = ne_all[np.isfinite(ne_all) & (ne_all > 0)]
    if finite.size > 0 and np.max(finite) / np.clip(np.min(finite), 1e-300, None) > 100:
        ax.set_yscale("log")
        ax.set_ylabel(r"Electron density  $n_e(s)$  (arb. units, log scale)", fontsize=13)

    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Saved: {save}")
        plt.close(fig)
    else:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot time-evolution of spatially-resolved electron density.",
    )
    parser.add_argument(
        "--data", type=str, required=True,
        help="Path to a fkrplk output file.",
    )
    parser.add_argument(
        "--save", type=str, default=None,
        help="Save figure to this path instead of displaying.",
    )
    parser.add_argument(
        "--cmap", type=str, default="viridis",
        help="Colourmap name (default: viridis).",
    )
    parser.add_argument(
        "--time-indices", type=str, default=None,
        help="Comma-separated time indices to plot (e.g. '0,2,4,6,8'). "
             "Defaults to ~5 evenly-spaced snapshots.",
    )
    return parser


def main() -> None:
    """Entry point for the CLI."""
    args = _build_parser().parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tdfp_reader import read_fkrplk_output

    data = read_fkrplk_output(args.data)

    n_times = data.phi.shape[0]

    # Select time indices
    if args.time_indices is not None:
        t_indices = [int(x) for x in args.time_indices.split(",")]
    else:
        # Pick ~5 evenly-spaced snapshots including first and last meaningful
        # Detect last meaningful time step (finite phi values with reasonable range)
        last_good = n_times - 1
        for ti in range(n_times - 1, -1, -1):
            p = data.phi[ti]
            finite = p[np.isfinite(p)]
            if finite.size > 0 and np.max(np.abs(finite)) > 1e-20:
                last_good = ti
                break
        n_sel = min(5, last_good + 1)
        t_indices = [int(i) for i in np.linspace(0, last_good, n_sel)]

    # Compute density for each selected time step
    ne_all = np.zeros((len(t_indices), data.phi.shape[1]))
    time_labels: list[str] = []

    for i, ti in enumerate(t_indices):
        ne_all[i] = compute_density(data.phi[ti], data.ekev, data.amu)
        # Label with iteration number
        iter_val = data.iterations[ti]
        frac = ti / max(n_times - 1, 1)
        time_labels.append(
            f"t[{ti}]  (iter {int(iter_val):,}, {frac:.0%})"
        )

    plot_density_evolution(
        data.soft, ne_all, time_labels,
        cmap=args.cmap, save=args.save,
    )


if __name__ == "__main__":
    main()
