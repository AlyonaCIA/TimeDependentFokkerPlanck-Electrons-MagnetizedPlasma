"""3-D coronal loop coloured by Hard X-ray emissivity (matplotlib).

Maps the 1-D spatial coordinate *s* onto a semi-circular coronal-loop
geometry and renders it as a glowing scatter-plot tube whose colour and
size encode :math:`\\log_{10} I_{\\mathrm{HXR}}(s)`, simulating what the
RHESSI / Solar Orbiter STIX satellite would observe.

Hard X-ray emissivity is approximated via thin-target bremsstrahlung:

.. math::

    I_{\\mathrm{HXR}}(s) \\;\\propto\\; n_e(s)
    \\int_{E_{\\min}}^{E_{\\max}}
    \\frac{1}{E}
    \\int_{-1}^{1} f(E, \\mu, s)\\, d\\mu\\, dE

where the :math:`1/E` weighting reflects the Kramers bremsstrahlung
cross-section.

Usage
-----
Real data::

    uv run python scripts/python/plot_3d_xray_loop.py \\
        --data experiments/conf_original/fkrplk.test

Save a figure::

    uv run python scripts/python/plot_3d_xray_loop.py \\
        --data experiments/conf_original/fkrplk.test \\
        --save output/figures/08_xray_coronal_loop.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Path3DCollection  # noqa: F401


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

def compute_xray_emissivity(
    phi_t: np.ndarray,
    ekev: np.ndarray,
    amu: np.ndarray,
    n_bg: float = 1.0e11,
) -> np.ndarray:
    r"""Compute thin-target HXR emissivity as a function of position.

    .. math::

        I_{\mathrm{HXR}}(s) \propto n_{\mathrm{bg}}
        \int \frac{1}{E}
        \int_{-1}^{1} f(E,\mu,s)\, d\mu\, dE

    Parameters
    ----------
    phi_t : ndarray, shape (n_space, n_pitch, n_energy)
        Electron distribution function at one time snapshot.
    ekev : ndarray, shape (n_energy,)
        Energy grid in keV.
    amu : ndarray, shape (n_pitch,)
        Pitch-angle cosine grid.
    n_bg : float
        Background electron density in cm⁻³ (for proportionality).

    Returns
    -------
    I_xray : ndarray, shape (n_space,)
        Hard X-ray emissivity at each spatial point (arb. units).
    """
    phi_clean = np.where(np.isfinite(phi_t), phi_t, 0.0)

    # Integrate over pitch-angle: shape → (n_space, n_energy)
    f_E_s = np.trapezoid(phi_clean, x=amu, axis=1)

    # Weight by 1/E (Kramers cross-section) and integrate over energy
    # Guard against division by zero at E=0
    inv_E = np.where(ekev > 0, 1.0 / ekev, 0.0)
    integrand = f_E_s * inv_E[np.newaxis, :]  # (n_space, n_energy)

    I_xray = n_bg * np.trapezoid(integrand, x=ekev, axis=-1)  # (n_space,)
    return np.abs(I_xray)


def _mock_emissivity(n_points: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """Synthetic HXR emissivity peaking at the two footpoints.

    This is what the LP81 / LP88 results predict: the hard X-ray
    emission is concentrated at the dense chromospheric footpoints.
    """
    loop_length = 50.0  # Mm
    s = np.linspace(0.0, loop_length, n_points)
    t = s / loop_length

    I_xray = (
        1.0e4 * np.exp(-((t - 0.0) ** 2) / (2 * 0.05**2))
        + 1.0e4 * np.exp(-((t - 1.0) ** 2) / (2 * 0.05**2))
        + 5.0e1  # tenuous coronal emission
        + 2.0e2 * np.exp(-((t - 0.5) ** 2) / (2 * 0.12**2))  # weak looptop
    )
    return s, I_xray


def _emissivity_from_simulation(
    filepath: str | Path,
    time_index: int = -1,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute HXR emissivity from Fortran output data.

    Returns
    -------
    s : ndarray
        Spatial coordinate in Mm.
    I_xray : ndarray
        Hard X-ray emissivity (arb. units).
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tdfp_reader import read_fkrplk_output

    data = read_fkrplk_output(filepath)
    phi_t = data.phi[time_index]  # (space, pitch, energy)

    I_xray = compute_xray_emissivity(
        phi_t, data.ekev, data.amu, n_bg=data.rhomax,
    )
    # Convert spatial coordinate cm → Mm
    s = data.soft / 1.0e8
    return s, I_xray


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def _loop_xyz(s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map arc-length *s* onto a 3-D semi-circular loop (x–z plane).

    Returns x, y, z arrays.  Footpoints sit on the z = 0 plane.
    """
    R = s[-1] / np.pi
    theta = s / R  # [0, π]
    x = R * np.cos(theta)
    y = np.zeros_like(s)
    z = R * np.sin(theta)
    return x, y, z


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def plot_xray_loop(
    s: np.ndarray,
    I_xray: np.ndarray,
    *,
    save: str | None = None,
    cmap: str = "hot",
    dpi: int = 200,
    figsize: tuple[float, float] = (12, 9),
) -> Figure:
    r"""Render a 3-D coronal loop coloured by HXR emissivity.

    Parameters
    ----------
    s : ndarray, shape (N,)
        Arc-length coordinate (Mm).
    I_xray : ndarray, shape (N,)
        Hard X-ray emissivity at each spatial point.
    save : str | None
        If given, save the figure to this path.
    cmap : str
        Colourmap (default: ``'hot'``).
    dpi : int
        Figure resolution.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    Figure
        The matplotlib Figure object.
    """
    x, y, z = _loop_xyz(s)
    R = s[-1] / np.pi

    # Log-scale emissivity for colour mapping
    I_floor = I_xray[I_xray > 0].min() if np.any(I_xray > 0) else 1.0
    I_safe = np.clip(I_xray, I_floor, None)
    log_I = np.log10(I_safe)

    # Normalise to [0, 1] for size/alpha scaling
    log_min, log_max = log_I.min(), log_I.max()
    if log_max - log_min < 1e-12:
        log_norm = np.ones_like(log_I)
    else:
        log_norm = (log_I - log_min) / (log_max - log_min)

    # Point sizes: bigger where emission is bright
    sizes = 8 + 120 * log_norm**2

    # Build RGBA colours with per-point alpha baked in
    colormap = plt.colormaps[cmap]
    rgba_core = colormap(log_norm).copy()
    rgba_core[:, 3] = 0.15 + 0.85 * log_norm**1.5   # bright at footpoints

    rgba_mid = colormap(log_norm).copy()
    rgba_mid[:, 3] = 0.15

    rgba_glow = colormap(log_norm).copy()
    rgba_glow[:, 3] = 0.06

    # ── Figure setup ────────────────────────────────────────────────────
    fig = plt.figure(figsize=figsize, facecolor="black")
    ax = fig.add_subplot(111, projection="3d", facecolor="black")

    # ── Outer glow layer (large, very transparent) ──────────────────────
    ax.scatter(
        x, y, z,
        color=rgba_glow,
        s=sizes * 4.5,
        edgecolors="none",
        depthshade=True,
    )

    # ── Mid glow layer ──────────────────────────────────────────────────
    ax.scatter(
        x, y, z,
        color=rgba_mid,
        s=sizes * 2.0,
        edgecolors="none",
        depthshade=True,
    )

    # ── Core emission layer ─────────────────────────────────────────────
    sc = ax.scatter(
        x, y, z,
        c=log_I,
        cmap=cmap,
        s=sizes,
        edgecolors="none",
        depthshade=True,
    )

    # ── Footpoint markers (bright spheres) ──────────────────────────────
    for fx, fy, fz in [(x[0], y[0], z[0]), (x[-1], y[-1], z[-1])]:
        ax.scatter(
            [fx], [fy], [fz],
            c="#FF4500",
            s=600,
            alpha=0.95,
            edgecolors="#FFD700",
            linewidths=1.5,
            zorder=10,
            depthshade=False,
        )
        # Extra glow ring around footpoints
        ax.scatter(
            [fx], [fy], [fz],
            c="#FF6000",
            s=2000,
            alpha=0.12,
            edgecolors="none",
            depthshade=False,
        )

    # ── Chromosphere surface (z = 0 plane) ──────────────────────────────
    extent = R * 1.5
    xx_plane, yy_plane = np.meshgrid(
        np.linspace(-extent, extent, 2),
        np.linspace(-extent * 0.7, extent * 0.7, 2),
    )
    zz_plane = np.zeros_like(xx_plane)
    ax.plot_surface(
        xx_plane,
        yy_plane,
        zz_plane,
        color="#1A0800",
        alpha=0.45,
        shade=False,
        zorder=0,
    )

    # ── Colour-bar ──────────────────────────────────────────────────────
    cbar = fig.colorbar(sc, ax=ax, shrink=0.55, aspect=20, pad=0.02)
    cbar.set_label(
        r"$\log_{10}\; I_{\mathrm{HXR}}$  [arb. units]",
        fontsize=13,
        color="white",
        labelpad=10,
    )
    cbar.ax.yaxis.set_tick_params(color="white", labelcolor="white", labelsize=10)
    cbar.outline.set_edgecolor("white")
    cbar.outline.set_linewidth(0.5)

    # ── Axes styling (space theme) ──────────────────────────────────────
    ax.set_xlabel(r"$x$  [Mm]", fontsize=12, color="white", labelpad=10)
    ax.set_ylabel(r"$y$  [Mm]", fontsize=12, color="white", labelpad=10)
    ax.set_zlabel(r"$z$  [Mm]", fontsize=12, color="white", labelpad=10)

    ax.set_title(
        r"Coronal Loop — Hard X-ray Emissivity  $I_{\mathrm{HXR}}(s)$",
        fontsize=16,
        color="white",
        pad=20,
        fontweight="bold",
    )

    # Remove grid lines and set tick colours
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("gray")
    ax.yaxis.pane.set_edgecolor("gray")
    ax.zaxis.pane.set_edgecolor("gray")
    ax.xaxis.pane.set_alpha(0.05)
    ax.yaxis.pane.set_alpha(0.05)
    ax.zaxis.pane.set_alpha(0.05)

    ax.tick_params(axis="x", colors="white", labelsize=9)
    ax.tick_params(axis="y", colors="white", labelsize=9)
    ax.tick_params(axis="z", colors="white", labelsize=9)

    # Camera angle: slightly elevated, looking at the loop from the side
    ax.view_init(elev=25, azim=-60)
    ax.set_box_aspect([1.6, 0.8, 1.0])

    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight", facecolor="black")
        print(f"Saved: {save}")

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "3-D coronal loop coloured by Hard X-ray emissivity "
            "(matplotlib scatter, RHESSI-style)."
        ),
    )
    p.add_argument(
        "--data",
        type=str,
        default=None,
        help="Path to a fkrplk output file.  If omitted, uses mock data.",
    )
    p.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time-snapshot index (default: -1 = last meaningful).",
    )
    p.add_argument(
        "--save",
        type=str,
        default=None,
        help="Save figure to this path.",
    )
    p.add_argument(
        "--cmap",
        type=str,
        default="hot",
        help="Colourmap name (default: hot).",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Figure DPI (default: 200).",
    )
    return p


def main() -> None:
    """Entry point for the CLI."""
    args = _build_parser().parse_args()

    if args.data is not None:
        s, I_xray = _emissivity_from_simulation(args.data, args.time_index)
    else:
        s, I_xray = _mock_emissivity()

    fig = plot_xray_loop(s, I_xray, save=args.save, cmap=args.cmap, dpi=args.dpi)

    if args.save is None:
        plt.show()


if __name__ == "__main__":
    main()
