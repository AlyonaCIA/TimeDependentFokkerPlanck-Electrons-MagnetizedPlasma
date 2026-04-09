"""Recreate the classic LP88 contour plots of f(E, mu) from Fokker-Planck output.

Plots filled + line contours of log10(f) in the Energy (keV) vs. Pitch-Angle
(degrees) plane, reproducing the style from Leach & Petrosian (1988).  Multiple
spatial positions and time snapshots can be displayed as a panel grid.

Usage
-----
Mock data (no simulation file needed)::

    uv run python scripts/python/plot_2d_contours.py

Real Fortran output::

    uv run python scripts/python/plot_2d_contours.py \\
        --data experiments/conf_original/fkrplk.test --time-index 0

Save to file::

    uv run python scripts/python/plot_2d_contours.py --save contours.png
"""

from __future__ import annotations

# -- Standard library --------------------------------------------------------
import argparse
import sys
from pathlib import Path

# -- Third-party -------------------------------------------------------------
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure


# ---------------------------------------------------------------------------
# Mock data generation
# ---------------------------------------------------------------------------

def _mock_fp_data(
    n_energy: int = 80,
    n_pitch: int = 60,
    n_space: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate a synthetic f(E, mu, s) mimicking LP88 Figure 2.

    The mock distribution has:
    - Power-law decay in energy.
    - Anisotropy peaking at small pitch angles (beam along the field).
    - Spatial evolution: beam becomes more isotropic at larger column depths.

    Returns
    -------
    ekev : ndarray, shape (n_energy,)
        Energy grid in keV.
    pad : ndarray, shape (n_pitch,)
        Pitch angle in degrees (0–180).
    space_labels : ndarray of str, shape (n_space,)
        Labels for each spatial position.
    phi : ndarray, shape (n_space, n_pitch, n_energy)
        Distribution function f(s, mu, E).
    """
    ekev = np.logspace(np.log10(10.0), np.log10(1e4), n_energy)
    pad = np.linspace(0.0, 180.0, n_pitch)
    mu = np.cos(np.radians(pad))

    # Column-depth labels (in units of 10^19 cm^-2)
    tau_vals = np.array([0.0, 0.5, 1.0, 2.0, 5.0])[:n_space]
    space_labels = np.array([rf"$\tau = {t:.1f}\times10^{{19}}$" for t in tau_vals])

    phi = np.zeros((n_space, n_pitch, n_energy))

    for i_s, tau in enumerate(tau_vals):
        # Isotropisation factor: beam narrows at small tau, broadens at large
        sigma_mu = 0.3 + 0.5 * (tau / 5.0)

        for i_mu, mu_val in enumerate(mu):
            # Angular part: Gaussian beam centred on mu=1 (forward)
            angular = np.exp(-((mu_val - 1.0) ** 2) / (2.0 * sigma_mu**2))
            # Add a small isotropic floor
            angular = angular + 0.02

            # Energy part: broken power-law
            e_break = 100.0  # keV
            spectral = np.where(
                ekev < e_break,
                (ekev / e_break) ** (-1.5),
                (ekev / e_break) ** (-3.5 - 0.3 * tau),
            )
            phi[i_s, i_mu, :] = 1e8 * angular * spectral

    return ekev, pad, space_labels, phi


def _load_simulation_data(
    filepath: str | Path,
    time_index: int,
    space_indices: list[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load f(E, mu) slices from a Fortran output file.

    Returns the same tuple format as ``_mock_fp_data``.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tdfp_reader import read_fkrplk_output

    data = read_fkrplk_output(filepath)

    if space_indices is None:
        # Pick 5 evenly-spaced spatial positions
        n_s = data.ntau + 1
        space_indices = [int(i) for i in np.linspace(0, n_s - 1, 5)]

    phi_t = data.phi[time_index]  # (space, pitch, energy)
    phi_sel = phi_t[space_indices, :, :]  # (n_sel, pitch, energy)

    # Spatial labels from soft (position in cm → Mm)
    space_labels = np.array(
        [rf"$s = {data.soft[j] / 1e8:.1f}$ Mm" for j in space_indices]
    )

    return data.ekev, data.pad, space_labels, phi_sel


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _style_ax(ax: matplotlib.axes.Axes) -> None:
    """Apply consistent tick styling to *ax*."""
    ax.tick_params(axis="both", which="major", labelsize=11)
    ax.tick_params(axis="both", which="minor", labelsize=9)


def plot_lp88_contours(
    ekev: np.ndarray,
    pad: np.ndarray,
    space_labels: np.ndarray,
    phi: np.ndarray,
    *,
    n_levels: int = 15,
    cmap: str = "inferno",
    save: str | None = None,
) -> Figure:
    """Plot LP88-style contour panels of log10(f) in E–pitch-angle space.

    Parameters
    ----------
    ekev : ndarray, shape (n_energy,)
        Energy grid in keV.
    pad : ndarray, shape (n_pitch,)
        Pitch angle in degrees.
    space_labels : ndarray of str
        Label for each spatial panel.
    phi : ndarray, shape (n_panels, n_pitch, n_energy)
        Distribution function slices.
    n_levels : int
        Number of contour levels.
    cmap : str
        Colourmap name.
    save : str | None
        Save the figure to this path if given.

    Returns
    -------
    Figure
    """
    n_panels = len(space_labels)
    ncols = min(n_panels, 3)
    nrows = (n_panels + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(6.0 * ncols, 5.0 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    fig.suptitle(
        r"LP88 Contours — $\log_{10}\,f(E,\,\mu)$",
        fontsize=16, fontweight="bold", y=1.02,
    )
    axes_flat = axes.ravel()

    # Shared contour levels across all panels
    with np.errstate(divide="ignore", invalid="ignore"):
        log_phi = np.log10(np.clip(phi, 1e-30, None))

    vmin = float(np.nanmin(log_phi[np.isfinite(log_phi)]))
    vmax = float(np.nanmax(log_phi[np.isfinite(log_phi)]))
    levels = np.linspace(vmin, vmax, n_levels)

    # Meshgrid for contour: x = pitch angle, y = energy (log-scale)
    PA, E = np.meshgrid(pad, ekev, indexing="ij")

    for idx in range(len(axes_flat)):
        ax = axes_flat[idx]
        if idx >= n_panels:
            ax.set_visible(False)
            continue

        Z = log_phi[idx]

        # Filled contours
        cf = ax.contourf(PA, E, Z, levels=levels, cmap=cmap, extend="both")
        # Line contours (white for visibility)
        ax.contour(
            PA, E, Z,
            levels=levels[::2],
            colors="white",
            linewidths=0.5,
            alpha=0.6,
        )

        ax.set_yscale("log")
        ax.set_xlim(pad.min(), pad.max())
        ax.set_ylim(ekev.min(), ekev.max())

        ax.set_xlabel("Pitch Angle  (degrees)", fontsize=12)
        ax.set_ylabel("Energy  (keV)", fontsize=12)

        # Subtle grid for readability
        ax.grid(True, which="major", ls=":", lw=0.4, color="white", alpha=0.25)

        # Panel label
        ax.text(
            0.03, 0.95,
            space_labels[idx],
            transform=ax.transAxes,
            fontsize=11,
            color="white",
            va="top",
            ha="left",
            bbox=dict(facecolor="black", alpha=0.6, edgecolor="none", pad=3),
        )
        _style_ax(ax)

    # Shared colourbar
    cbar = fig.colorbar(
        cf, ax=axes_flat[:n_panels].tolist(),
        label=r"$\log_{10}\, f$",
        pad=0.02,
        shrink=0.85,
    )
    cbar.ax.tick_params(labelsize=11)

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
        description="LP88-style 2D contour plots of f(E, pitch-angle).",
    )
    parser.add_argument(
        "--data", type=str, default=None,
        help="Path to a fkrplk output file. Uses mock data if omitted.",
    )
    parser.add_argument(
        "--time-index", type=int, default=0,
        help="Time-snapshot index to plot (default: 0).",
    )
    parser.add_argument(
        "--space-indices", type=str, default=None,
        help="Comma-separated spatial indices (e.g. '0,10,30,60,100').",
    )
    parser.add_argument(
        "--save", type=str, default=None,
        help="Save figure to this path instead of displaying.",
    )
    parser.add_argument(
        "--cmap", type=str, default="inferno",
        help="Colourmap name (default: inferno).",
    )
    parser.add_argument(
        "--levels", type=int, default=15,
        help="Number of contour levels (default: 15).",
    )
    return parser


def main() -> None:
    """Entry point for the CLI."""
    args = _build_parser().parse_args()

    space_idx = None
    if args.space_indices:
        space_idx = [int(x) for x in args.space_indices.split(",")]

    if args.data is not None:
        ekev, pad, labels, phi = _load_simulation_data(
            args.data,
            time_index=args.time_index,
            space_indices=space_idx,
        )
    else:
        ekev, pad, labels, phi = _mock_fp_data()

    plot_lp88_contours(
        ekev, pad, labels, phi,
        n_levels=args.levels,
        cmap=args.cmap,
        save=args.save,
    )


if __name__ == "__main__":
    main()
