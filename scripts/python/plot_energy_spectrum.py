"""Energy spectrum cooling — total flux vs. energy at multiple times.

Computes the total energy spectrum by integrating the distribution function
over pitch-angle and space:

.. math::

    F(E, t) = \\int_0^L \\int_{-1}^{1} f(E, \\mu, s, t)\\, d\\mu\\, ds

and plots :math:`F(E)` on a log–log scale for several time snapshots,
coloured with a sequential colourmap (viridis) from early (dark) to late
(yellow).  The progressive steepening of the spectrum reveals **collisional
cooling**: Coulomb collisions deplete low-energy electrons faster than
high-energy ones.

Usage
-----
Real Fortran output::

    uv run python scripts/python/plot_energy_spectrum.py \\
        --data experiments/conf_original/fkrplk.test

Save to file::

    uv run python scripts/python/plot_energy_spectrum.py \\
        --data experiments/conf_original/fkrplk.test \\
        --save output/figures/07_energy_spectrum.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure


# ---------------------------------------------------------------------------
# Physics: integrate f over pitch-angle and space -> F(E)
# ---------------------------------------------------------------------------

def compute_energy_spectrum(
    phi_t: np.ndarray,
    amu: np.ndarray,
    soft: np.ndarray,
) -> np.ndarray:
    """Integrate f(E, mu, s) over pitch-angle cosine and space.

    Parameters
    ----------
    phi_t : ndarray, shape (n_space, n_pitch, n_energy)
        Distribution function at a single time snapshot.
    amu : ndarray, shape (n_pitch,)
        Pitch-angle cosine grid (μ = cos θ), ordered −1 … +1.
    soft : ndarray, shape (n_space,)
        Spatial position along the loop in cm.

    Returns
    -------
    F_E : ndarray, shape (n_energy,)
        Total energy spectrum (integrated flux).
    """
    # Replace NaN/Inf with 0 for integration
    phi_clean = np.where(np.isfinite(phi_t), phi_t, 0.0)

    # Integrate over pitch-angle (axis=1), then over space (axis=0)
    integrated_mu = np.trapezoid(phi_clean, x=amu, axis=1)  # (space, energy)
    F_E = np.trapezoid(integrated_mu, x=soft, axis=0)        # (energy,)
    return F_E


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_energy_spectrum(
    ekev: np.ndarray,
    spectra: np.ndarray,
    time_labels: list[str],
    *,
    cmap: str = "viridis",
    save: str | None = None,
) -> Figure:
    """Log-log plot of F(E) for multiple time steps.

    Parameters
    ----------
    ekev : ndarray, shape (n_energy,)
        Energy grid in keV.
    spectra : ndarray, shape (n_times, n_energy)
        Total energy spectrum at each time snapshot.
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

    n_lines = len(time_labels)
    cm = plt.get_cmap(cmap)
    colours = [cm(i / max(n_lines - 1, 1)) for i in range(n_lines)]

    for i, (spec, label) in enumerate(zip(spectra, time_labels)):
        # Mask non-positive values for log-log
        pos = spec > 0
        if np.any(pos):
            ax.plot(ekev[pos], spec[pos], color=colours[i], lw=2.2, label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy  $E$  (keV)", fontsize=13)
    ax.set_ylabel(r"Total spectrum  $F(E)$  (arb. units)", fontsize=13)
    ax.set_title(
        r"Energy Spectrum Cooling — $F(E,\,t) = \int f\, d\mu\, ds$",
        fontsize=15, fontweight="bold",
    )

    ax.legend(fontsize=11, title="Time step", title_fontsize=12)
    ax.tick_params(labelsize=11)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

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
        description="Plot energy spectrum cooling F(E, t) on a log-log scale.",
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
        last_good = n_times - 1
        for ti in range(n_times - 1, -1, -1):
            p = data.phi[ti]
            finite = p[np.isfinite(p)]
            if finite.size > 0 and np.max(np.abs(finite)) > 1e-20:
                last_good = ti
                break
        n_sel = min(5, last_good + 1)
        t_indices = [int(i) for i in np.linspace(0, last_good, n_sel)]

    # Compute spectrum for each selected time step
    spectra = np.zeros((len(t_indices), data.phi.shape[3]))
    time_labels: list[str] = []

    for i, ti in enumerate(t_indices):
        spectra[i] = compute_energy_spectrum(data.phi[ti], data.amu, data.soft)
        iter_val = data.iterations[ti]
        frac = ti / max(n_times - 1, 1)
        time_labels.append(
            f"t[{ti}]  (iter {int(iter_val):,}, {frac:.0%})"
        )

    plot_energy_spectrum(
        data.ekev, spectra, time_labels,
        cmap=args.cmap, save=args.save,
    )


if __name__ == "__main__":
    main()
