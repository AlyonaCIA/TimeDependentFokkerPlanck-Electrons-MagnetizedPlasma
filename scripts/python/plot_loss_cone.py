"""Loss-cone phase-space diagram — trapped vs. precipitating electrons.

Plots the electron distribution function :math:`f(E, \\mu)` at selected
spatial positions as filled contours in the :math:`(\\mu, E)` plane,
with the loss-cone boundary overlaid.

The **loss-cone angle** at position :math:`s` is:

.. math::

    \\mu_{\\mathrm{loss}}(s)
    = \\pm\\sqrt{1 - \\frac{B(s)}{B_{\\max}}}

Electrons with :math:`|\\mu| > \\mu_{\\mathrm{loss}}` have their
velocity nearly parallel to **B** and can stream to the footpoints
without mirroring — they **precipitate** into the dense chromosphere
and are lost.  Electrons with :math:`|\\mu| < \\mu_{\\mathrm{loss}}`
are reflected by the converging field and remain **trapped** in the
coronal portion of the loop.

For a mirror ratio :math:`r_m = B_{\\max}/B_{\\min} = 2`
(the default in this simulation), the loss-cone boundary at the loop
apex (:math:`B = B_{\\min}`) is
:math:`\\mu_{\\mathrm{loss}} = \\pm\\sqrt{1 - 1/r_m} \\approx \\pm 0.707`.

Usage
-----
Real data::

    uv run python scripts/python/plot_loss_cone.py \\
        --data experiments/conf_original/fkrplk.test --time-index 2

Save a figure::

    uv run python scripts/python/plot_loss_cone.py \\
        --data experiments/conf_original/fkrplk.test --time-index 2 \\
        --save output/figures/10_loss_cone.png

Mock data (no simulation file needed)::

    uv run python scripts/python/plot_loss_cone.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.figure import Figure


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

def loss_cone_boundary(B_local: float, B_max: float) -> float:
    r"""Return the loss-cone cosine at a given position.

    .. math::

        \mu_{\mathrm{loss}} = \sqrt{1 - B(s) / B_{\max}}

    Parameters
    ----------
    B_local : float
        Magnetic field strength at position *s*.
    B_max : float
        Maximum field strength along the loop (at footpoints).

    Returns
    -------
    float
        Absolute value of the loss-cone cosine :math:`|\mu_c|`.
        Returns 0 if :math:`B(s) \ge B_{\max}`.
    """
    ratio = B_local / B_max
    if ratio >= 1.0:
        return 0.0
    return np.sqrt(1.0 - ratio)


# ---------------------------------------------------------------------------
# Mock data
# ---------------------------------------------------------------------------

def _mock_loss_cone(
    n_energy: int = 80,
    n_pitch: int = 100,
    mirror_ratio: float = 2.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Generate a synthetic f(E, μ) with a clear loss-cone depletion.

    Returns
    -------
    ekev : ndarray, shape (n_energy,)
    amu : ndarray, shape (n_pitch,)
    f_E_mu : ndarray, shape (n_pitch, n_energy)
    mu_loss : float
        Loss-cone cosine boundary.
    """
    ekev = np.logspace(np.log10(10), np.log10(5000), n_energy)
    amu = np.linspace(-1, 1, n_pitch)

    mu_loss = np.sqrt(1.0 - 1.0 / mirror_ratio)

    MU, E = np.meshgrid(amu, ekev, indexing="ij")

    # Trapped population: strong in |μ| < μ_loss
    trapped = 1e8 * (E / 10.0) ** (-2.5) * np.exp(
        -(MU**2) / (2 * 0.5**2)
    )

    # Loss-cone depletion: suppress particles with |μ| > μ_loss
    depletion = np.where(
        np.abs(MU) <= mu_loss,
        1.0,
        np.exp(-((np.abs(MU) - mu_loss) ** 2) / (2 * 0.05**2)) * 0.1,
    )

    # Small isotropic floor
    floor = 1e2 * (E / 10.0) ** (-4.0)

    f_E_mu = trapped * depletion + floor

    return ekev, amu, f_E_mu, mu_loss


# ---------------------------------------------------------------------------
# Load simulation data
# ---------------------------------------------------------------------------

def _load_simulation(
    filepath: str | Path,
    time_index: int = -1,
    space_index: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Load f(E, μ) from simulation and compute the loss-cone boundary.

    Parameters
    ----------
    filepath : path
        Path to fkrplk output file.
    time_index : int
        Time snapshot to use.
    space_index : int
        Spatial position index (0 = loop top / injection).

    Returns
    -------
    ekev, amu, f_E_mu, mu_loss
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tdfp_reader import read_fkrplk_output

    data = read_fkrplk_output(filepath)
    phi_t = data.phi[time_index]       # (space, pitch, energy)
    f_E_mu = phi_t[space_index, :, :]  # (pitch, energy)

    # Compute loss-cone boundary at this spatial position
    B_local = data.bm[space_index]
    B_max = data.bm.max()
    mu_loss = loss_cone_boundary(B_local, B_max)

    return data.ekev, data.amu, f_E_mu, mu_loss


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_loss_cone(
    ekev: np.ndarray,
    amu: np.ndarray,
    f_E_mu: np.ndarray,
    mu_loss: float,
    *,
    title: str | None = None,
    cmap: str = "magma",
    n_levels: int = 20,
    save: str | None = None,
    dpi: int = 300,
    figsize: tuple[float, float] = (10, 7),
) -> Figure:
    r"""Plot loss-cone phase-space diagram.

    Parameters
    ----------
    ekev : ndarray, shape (n_energy,)
        Energy grid in keV.
    amu : ndarray, shape (n_pitch,)
        Pitch-angle cosine grid (:math:`\mu = \cos\theta`).
    f_E_mu : ndarray, shape (n_pitch, n_energy)
        Electron distribution function.
    mu_loss : float
        Absolute value of the loss-cone cosine boundary.
    title : str | None
        Figure title override.
    cmap : str
        Colourmap (default: ``'magma'``).
    n_levels : int
        Number of contour levels.
    save : str | None
        Save figure to this path.
    dpi : int
        Figure resolution.
    figsize : tuple
        Figure size (width, height) in inches.

    Returns
    -------
    Figure
    """
    # -- Prepare log10(f) ----------------------------------------------------
    f_clean = np.where(np.isfinite(f_E_mu) & (f_E_mu > 0), f_E_mu, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_f = np.log10(f_clean)

    finite = log_f[np.isfinite(log_f)]
    if finite.size == 0:
        print("Warning: f(E,μ) is all zeros / NaN. Nothing to plot.")
        return plt.figure()

    vmin = float(np.nanmin(finite))
    vmax = float(np.nanmax(finite))
    if vmax - vmin < 1e-10:
        vmax = vmin + 1.0
    levels = np.linspace(vmin, vmax, n_levels)

    # Meshgrid: x = μ, y = E
    MU, E = np.meshgrid(amu, ekev, indexing="ij")

    # -- Figure --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # Filled contours
    cf = ax.contourf(MU, E, log_f, levels=levels, cmap=cmap, extend="both")

    # Thin white line contours for readability
    ax.contour(
        MU, E, log_f,
        levels=levels[::3],
        colors="white",
        linewidths=0.4,
        alpha=0.5,
    )

    ax.set_yscale("log")
    ax.set_xlim(amu.min(), amu.max())
    ax.set_ylim(ekev[ekev > 0].min(), ekev.max())

    # -- Loss-cone boundary lines --------------------------------------------
    erange = np.array([ekev[ekev > 0].min(), ekev.max()])

    for sign in (+1, -1):
        mu_c = sign * mu_loss
        ax.plot(
            [mu_c, mu_c], erange,
            color="cyan", lw=2.5, ls="--",
            zorder=5,
        )

    # -- Shade the loss-cone regions (|μ| > μ_loss → precipitating) ----------
    # Right loss cone: μ > +μ_loss
    ax.fill_betweenx(
        erange,
        mu_loss, amu.max(),
        color="gray", alpha=0.30, zorder=4,
        hatch="///", edgecolor="white", linewidth=0.0,
    )
    # Left loss cone: μ < -μ_loss
    ax.fill_betweenx(
        erange,
        amu.min(), -mu_loss,
        color="gray", alpha=0.30, zorder=4,
        hatch="///", edgecolor="white", linewidth=0.0,
    )

    # -- Annotation labels on regions ----------------------------------------
    mid_trapped_y = np.sqrt(erange[0] * erange[1])  # geometric mean in log

    ax.text(
        0.0, mid_trapped_y * 2.5,
        "TRAPPED",
        fontsize=14, fontweight="bold", color="white",
        ha="center", va="center",
        bbox=dict(facecolor="black", alpha=0.5, edgecolor="none", pad=4),
        zorder=6,
    )

    # "LOST" labels in both loss-cone wings
    for mu_pos in [(mu_loss + amu.max()) / 2, (amu.min() - mu_loss) / 2]:
        ax.text(
            mu_pos, mid_trapped_y * 2.5,
            "LOST",
            fontsize=12, fontweight="bold", color="#FF6666",
            ha="center", va="center",
            bbox=dict(facecolor="black", alpha=0.5, edgecolor="none", pad=3),
            zorder=6,
        )

    # -- μ = 0 reference line (perpendicular to B) ---------------------------
    ax.axvline(0, color="white", lw=0.8, ls=":", alpha=0.4, zorder=3)

    # -- Colour bar -----------------------------------------------------------
    cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=30, shrink=0.92)
    cbar.set_label(
        r"$\log_{10}\; f(E,\,\mu)$",
        fontsize=14,
    )
    cbar.ax.tick_params(labelsize=11)

    # -- Axes labels ----------------------------------------------------------
    ax.set_xlabel(
        r"Pitch-Angle Cosine  $\mu = \cos\theta$",
        fontsize=14,
    )
    ax.set_ylabel(
        r"Kinetic Energy  $E$  (keV)",
        fontsize=14,
    )

    if title is None:
        title = (
            r"Loss-Cone Phase-Space Diagram"
            rf" — $\mu_{{\mathrm{{loss}}}} = \pm{mu_loss:.3f}$"
        )
    ax.set_title(title, fontsize=16, fontweight="bold", pad=12)

    ax.tick_params(labelsize=12)
    ax.grid(True, which="major", ls=":", lw=0.3, alpha=0.3, color="white")

    # -- Legend ---------------------------------------------------------------
    trapped_patch = mpatches.Patch(
        facecolor="none", edgecolor="cyan", linestyle="--", linewidth=2,
        label=r"Loss-cone boundary $|\mu_c|$",
    )
    lost_patch = mpatches.Patch(
        facecolor="gray", alpha=0.3, edgecolor="white",
        hatch="///",
        label=r"Loss cone (precipitating)",
    )
    ax.legend(
        handles=[trapped_patch, lost_patch],
        loc="lower right",
        fontsize=11,
        framealpha=0.7,
        edgecolor="gray",
    )

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Saved: {save}")
        plt.close(fig)

    return fig


# ---------------------------------------------------------------------------
# Multi-panel: loss cone at several spatial positions
# ---------------------------------------------------------------------------

def plot_loss_cone_panels(
    ekev: np.ndarray,
    amu: np.ndarray,
    phi_slices: list[np.ndarray],
    mu_loss_vals: list[float],
    panel_labels: list[str],
    *,
    cmap: str = "magma",
    n_levels: int = 20,
    save: str | None = None,
    dpi: int = 300,
) -> Figure:
    """Plot loss-cone diagrams at multiple spatial positions as panels.

    Parameters
    ----------
    ekev : ndarray, shape (n_energy,)
    amu : ndarray, shape (n_pitch,)
    phi_slices : list of ndarray, each shape (n_pitch, n_energy)
    mu_loss_vals : list of float
        Loss-cone cosine at each position.
    panel_labels : list of str
        Label for each panel.
    """
    n = len(phi_slices)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(6.5 * ncols, 5.5 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    fig.suptitle(
        r"Loss-Cone Evolution Along the Coronal Loop",
        fontsize=17, fontweight="bold",
    )
    axes_flat = axes.ravel()

    # Global contour levels
    all_log = []
    for phi_s in phi_slices:
        f_pos = np.where(np.isfinite(phi_s) & (phi_s > 0), phi_s, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            lf = np.log10(f_pos)
        fin = lf[np.isfinite(lf)]
        if fin.size > 0:
            all_log.append(fin)
    if all_log:
        combined = np.concatenate(all_log)
        vmin, vmax = float(combined.min()), float(combined.max())
    else:
        vmin, vmax = -30.0, 0.0
    if vmax - vmin < 1e-10:
        vmax = vmin + 1.0
    levels = np.linspace(vmin, vmax, n_levels)

    MU, E = np.meshgrid(amu, ekev, indexing="ij")

    last_cf = None
    for idx in range(len(axes_flat)):
        ax = axes_flat[idx]
        if idx >= n:
            ax.set_visible(False)
            continue

        f_clean = np.where(
            np.isfinite(phi_slices[idx]) & (phi_slices[idx] > 0),
            phi_slices[idx], np.nan,
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            log_f = np.log10(f_clean)

        cf = ax.contourf(MU, E, log_f, levels=levels, cmap=cmap, extend="both")
        ax.contour(MU, E, log_f, levels=levels[::3],
                   colors="white", linewidths=0.3, alpha=0.4)
        last_cf = cf

        ax.set_yscale("log")
        ax.set_xlim(amu.min(), amu.max())
        ax.set_ylim(ekev[ekev > 0].min(), ekev.max())

        mulc = mu_loss_vals[idx]
        erange = np.array([ekev[ekev > 0].min(), ekev.max()])

        for sign in (+1, -1):
            ax.plot([sign * mulc, sign * mulc], erange,
                    color="cyan", lw=2, ls="--", zorder=5)

        ax.fill_betweenx(erange, mulc, amu.max(),
                         color="gray", alpha=0.25, hatch="///",
                         edgecolor="white", linewidth=0.0, zorder=4)
        ax.fill_betweenx(erange, amu.min(), -mulc,
                         color="gray", alpha=0.25, hatch="///",
                         edgecolor="white", linewidth=0.0, zorder=4)

        ax.set_xlabel(r"$\mu = \cos\theta$", fontsize=12)
        ax.set_ylabel(r"$E$  (keV)", fontsize=12)
        ax.set_title(
            panel_labels[idx] + rf"  ($\mu_c = {mulc:.3f}$)",
            fontsize=12,
        )
        ax.tick_params(labelsize=10)
        ax.grid(True, which="major", ls=":", lw=0.3, alpha=0.3, color="white")

    if last_cf is not None:
        cbar = fig.colorbar(
            last_cf, ax=axes_flat[:n].tolist(),
            label=r"$\log_{10}\, f(E,\,\mu)$",
            pad=0.02, shrink=0.85,
        )
        cbar.ax.tick_params(labelsize=11)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Saved: {save}")
        plt.close(fig)

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Loss-cone phase-space diagram f(E, μ) with mirror boundary.",
    )
    p.add_argument(
        "--data", type=str, default=None,
        help="Path to a fkrplk output file. Uses mock data if omitted.",
    )
    p.add_argument(
        "--time-index", type=int, default=-1,
        help="Time-snapshot index (default: -1 = last meaningful).",
    )
    p.add_argument(
        "--space-index", type=int, default=0,
        help="Spatial position index for single-panel mode (default: 0 = apex).",
    )
    p.add_argument(
        "--panels", action="store_true",
        help="Multi-panel mode: show loss cone at 5 spatial positions.",
    )
    p.add_argument(
        "--save", type=str, default=None,
        help="Save figure to this path.",
    )
    p.add_argument(
        "--cmap", type=str, default="magma",
        help="Colourmap (default: magma).",
    )
    p.add_argument(
        "--dpi", type=int, default=300,
        help="Figure DPI (default: 300).",
    )
    p.add_argument(
        "--mirror-ratio", type=float, default=None,
        help="Override mirror ratio for mock data (default: 2.0).",
    )
    return p


def main() -> None:
    """Entry point for the CLI."""
    args = _build_parser().parse_args()

    if args.data is not None:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from tdfp_reader import read_fkrplk_output

        data = read_fkrplk_output(args.data)

        if args.panels:
            # Multi-panel: 5 positions along the loop
            n_s = data.ntau + 1
            space_indices = [int(i) for i in np.linspace(0, n_s - 1, 5)]
            B_max = data.bm.max()

            phi_slices = []
            mu_loss_vals = []
            panel_labels = []
            for si in space_indices:
                phi_slices.append(data.phi[args.time_index, si, :, :])
                mu_loss_vals.append(loss_cone_boundary(data.bm[si], B_max))
                panel_labels.append(rf"$s = {data.soft[si]/1e8:.1f}$ Mm")

            plot_loss_cone_panels(
                data.ekev, data.amu,
                phi_slices, mu_loss_vals, panel_labels,
                cmap=args.cmap, save=args.save, dpi=args.dpi,
            )
        else:
            ekev, amu, f_E_mu, mu_loss = _load_simulation(
                args.data, args.time_index, args.space_index,
            )
            s_mm = data.soft[args.space_index] / 1e8
            title = (
                rf"Loss-Cone Diagram at $s = {s_mm:.1f}$ Mm"
                rf" — $\mu_{{\mathrm{{loss}}}} = \pm{mu_loss:.3f}$"
            )
            plot_loss_cone(
                ekev, amu, f_E_mu, mu_loss,
                title=title, cmap=args.cmap,
                save=args.save, dpi=args.dpi,
            )
    else:
        rm = args.mirror_ratio if args.mirror_ratio is not None else 2.0
        ekev, amu, f_E_mu, mu_loss = _mock_loss_cone(mirror_ratio=rm)
        plot_loss_cone(
            ekev, amu, f_E_mu, mu_loss,
            cmap=args.cmap, save=args.save, dpi=args.dpi,
        )

    if args.save is None:
        plt.show()


if __name__ == "__main__":
    main()
