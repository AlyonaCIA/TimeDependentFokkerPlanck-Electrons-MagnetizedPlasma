"""Publication-quality plotting for TDFP solver output.

Provides functions for the four fundamental 1-D projections of the electron
distribution function *f(E, μ, s, t)* and 2-D / 3-D visualisations.

Example
-------
>>> from tdfp_reader import read_fkrplk_output
>>> from tdfp_plots import plot_energy_spectrum, plot_2d_panels
>>> data = read_fkrplk_output("experiments/conf2_E10keV_10MeV/fkrplk.conf2")
>>> plot_energy_spectrum(data, save="energy.png")
"""

from __future__ import annotations

# -- Standard library --
from math import floor, log10
from typing import Sequence

# -- Third-party --
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.figure import Figure
from matplotlib.ticker import LinearLocator

# -- Local --
from tdfp_reader import FPOutput


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _sci(num: float, decimal_digits: int = 1) -> str:
    """Format *num* as LaTeX scientific notation (e.g. ``$1.2\\times 10^{3}$``)."""
    if num == 0:
        return r"$0$"
    exp = int(floor(log10(abs(num))))
    coeff = round(num / 10.0**exp, decimal_digits)
    return rf"${coeff:.{decimal_digits}f}\times 10^{{{exp:d}}}$"


def _style_ax(ax: matplotlib.axes.Axes) -> None:
    """Apply consistent tick styling to *ax*."""
    ax.tick_params(axis="both", which="major", labelsize=14)


# ---------------------------------------------------------------------------
# 1D projections
# ---------------------------------------------------------------------------
def plot_energy_spectrum(
    data: FPOutput,
    time_indices: Sequence[int] = (0, 1, 2, 3, 4),
    pitch_index: int = 0,
    space_index: int = 10,
    ax: matplotlib.axes.Axes | None = None,
    save: str | None = None,
) -> Figure:
    """Plot *f(E)* at fixed pitch angle and spatial position.

    Parameters
    ----------
    data : FPOutput
        Loaded simulation data.
    time_indices : Sequence[int]
        Indices into the time axis to overlay.
    pitch_index : int
        Index into the pitch-angle (``amu``/``pad``) array.
    space_index : int
        Index into the spatial (``tau``/``soft``) array.
    ax : matplotlib.axes.Axes | None
        Pre-existing axes; a new figure is created when *None*.
    save : str | None
        If given, save the figure to this file path.

    Returns
    -------
    Figure
        The matplotlib figure containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    for it in time_indices:
        label = f"t = {data.tr[it]:.2f} s"
        ax.plot(data.ekev, data.phi[it, space_index, pitch_index, :], label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy [keV]", fontsize=16)
    ax.set_ylabel("Normalized density", fontsize=16)
    ax.set_title("Electron Distribution — Energy Dependence", fontsize=18)

    ax.text(
        0.97, 0.97,
        rf"$\alpha$ = {data.pad[pitch_index]:.1f}°"
        + "\n"
        + f"s = {_sci(data.soft[space_index])} cm",
        transform=ax.transAxes, fontsize=13,
        va="top", ha="right",
        bbox=dict(boxstyle="round", fc="white", alpha=0.8),
    )

    locmin = mticker.LogLocator(
        base=10.0,
        subs=np.arange(0.1, 1.0, 0.1),
        numticks=12,
    )
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.legend(fontsize=12)
    _style_ax(ax)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


def plot_pitch_angle(
    data: FPOutput,
    time_indices: Sequence[int] = (0, 1, 2, 3, 4),
    energy_index: int = 10,
    space_index: int = 10,
    ax: matplotlib.axes.Axes | None = None,
    save: str | None = None,
) -> Figure:
    """Plot *f(μ)* at fixed energy and spatial position.

    Parameters
    ----------
    data : FPOutput
        Loaded simulation data.
    time_indices : Sequence[int]
        Time-axis indices to overlay.
    energy_index : int
        Index into the energy (``ekev``) array.
    space_index : int
        Index into the spatial (``soft``) array.
    ax : matplotlib.axes.Axes | None
        Pre-existing axes; creates a new figure when *None*.
    save : str | None
        If given, save the figure to this path.

    Returns
    -------
    Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    for it in time_indices:
        label = f"t = {data.tr[it]:.2f} s"
        ax.plot(data.pad, data.phi[it, space_index, :, energy_index], label=label)

    ax.set_yscale("log")
    ax.set_xlabel("Pitch Angle [degrees]", fontsize=16)
    ax.set_ylabel("Normalized density", fontsize=16)
    ax.set_title("Electron Distribution — Pitch-Angle Dependence", fontsize=18)

    ax.text(
        0.97, 0.97,
        f"E = {data.ekev[energy_index]:.2f} keV"
        + "\n"
        + f"s = {_sci(data.soft[space_index])} cm",
        transform=ax.transAxes, fontsize=13,
        va="top", ha="right",
        bbox=dict(boxstyle="round", fc="white", alpha=0.8),
    )

    ax.legend(fontsize=12)
    _style_ax(ax)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


def plot_spatial_profile(
    data: FPOutput,
    time_indices: Sequence[int] = (0, 1, 2, 3, 4),
    energy_index: int = 10,
    pitch_index: int = 0,
    ax: matplotlib.axes.Axes | None = None,
    save: str | None = None,
) -> Figure:
    """Plot *f(s)* at fixed energy and pitch angle.

    Parameters
    ----------
    data : FPOutput
        Loaded simulation data.
    time_indices : Sequence[int]
        Time-axis indices to overlay.
    energy_index : int
        Index into the energy array.
    pitch_index : int
        Index into the pitch-angle array.
    ax : matplotlib.axes.Axes | None
        Pre-existing axes; creates a new figure when *None*.
    save : str | None
        If given, save the figure to this path.

    Returns
    -------
    Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    for it in time_indices:
        label = f"t = {data.tr[it]:.2f} s"
        ax.plot(data.soft, data.phi[it, :, pitch_index, energy_index], label=label)

    ax.set_yscale("log")
    ax.set_xlabel("Distance [cm]", fontsize=16)
    ax.set_ylabel("Normalized density", fontsize=16)
    ax.set_title("Electron Distribution — Spatial Dependence", fontsize=18)

    ax.text(
        0.97, 0.97,
        f"E = {data.ekev[energy_index]:.2f} keV"
        + "\n"
        + rf"$\alpha$ = {data.pad[pitch_index]:.1f}°",
        transform=ax.transAxes, fontsize=13,
        va="top", ha="right",
        bbox=dict(boxstyle="round", fc="white", alpha=0.8),
    )

    ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax.legend(fontsize=12)
    _style_ax(ax)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


def plot_time_evolution(
    data: FPOutput,
    space_indices: Sequence[int] = (0, 5, 10, 20, 30),
    energy_index: int = 10,
    pitch_index: int = 0,
    ax: matplotlib.axes.Axes | None = None,
    save: str | None = None,
) -> Figure:
    """Plot *f(t)* at fixed energy and pitch angle for several positions.

    Parameters
    ----------
    data : FPOutput
        Loaded simulation data.
    space_indices : Sequence[int]
        Spatial-axis indices to overlay.
    energy_index : int
        Index into the energy array.
    pitch_index : int
        Index into the pitch-angle array.
    ax : matplotlib.axes.Axes | None
        Pre-existing axes; creates a new figure when *None*.
    save : str | None
        If given, save the figure to this path.

    Returns
    -------
    Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    for si in space_indices:
        s_km = data.soft[si] / 1e5
        label = f"s = {_sci(s_km)} km" if s_km != 0 else "s = 0 km"
        ax.plot(data.tr, data.phi[:, si, pitch_index, energy_index], label=label)

    ax.set_yscale("log")
    ax.set_xlabel("Time [s]", fontsize=16)
    ax.set_ylabel("Normalized density", fontsize=16)
    ax.set_title("Electron Distribution — Time Dependence", fontsize=18)

    ax.text(
        0.97, 0.97,
        f"E = {data.ekev[energy_index]:.2f} keV"
        + "\n"
        + rf"$\alpha$ = {data.pad[pitch_index]:.1f}°",
        transform=ax.transAxes, fontsize=13,
        va="top", ha="right",
        bbox=dict(boxstyle="round", fc="white", alpha=0.8),
    )

    ax.legend(fontsize=12, loc="lower right")
    _style_ax(ax)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 2D contour panels (pcolormesh)
# ---------------------------------------------------------------------------
def plot_2d_panels(
    data: FPOutput,
    energy_index: int = 10,
    pitch_index: int = 50,
    space_index: int = 55,
    cmap: str = "YlGnBu",
    save: str | None = None,
) -> Figure:
    """Three-panel 2-D colour map: pitch vs time, energy vs time, distance vs time.

    All panels use ``log10(phi)`` for the colour scale.

    Parameters
    ----------
    data : FPOutput
        Loaded simulation data.
    energy_index, pitch_index, space_index : int
        Fixed-axis indices for the cross-sections shown.
    cmap : str
        Matplotlib colour-map name.
    save : str | None
        If given, save the figure to this path.

    Returns
    -------
    Figure
    """
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    phi_t = data.phi.transpose((3, 2, 1, 0))  # → (energy, pitch, space, time)

    with np.errstate(divide="ignore", invalid="ignore"):
        # -- Panel 1: Pitch Angle vs Time (fixed energy & distance) --
        z1 = np.log10(phi_t[energy_index, :, space_index, :])
        axes[0].pcolormesh(data.tr, data.pad, z1, cmap=cmap, shading="auto")
        axes[0].set_xlabel("Time [s]", fontsize=14)
        axes[0].set_ylabel("Pitch angle [deg]", fontsize=14)
        axes[0].set_title(
            f"E = {data.ekev[energy_index]:.0f} keV — "
            f"s = {data.soft[space_index]/1e5:.0f} km",
            fontsize=14,
        )
        axes[0].axhline(y=data.pad[pitch_index], color="b", ls=":", lw=2)

        # -- Panel 2: Energy vs Time (fixed pitch & distance) --
        z2 = np.log10(phi_t[:, pitch_index, space_index, :])
        axes[1].pcolormesh(data.tr, data.ekev, z2, cmap=cmap, shading="auto")
        axes[1].set_xlabel("Time [s]", fontsize=14)
        axes[1].set_ylabel("Energy [keV]", fontsize=14)
        axes[1].set_title(
            f"s = {data.soft[space_index]/1e5:.0f} km — "
            rf"$\alpha$ = {data.pad[pitch_index]:.0f}°",
            fontsize=14,
        )
        axes[1].axhline(y=data.ekev[energy_index], color="b", ls=":", lw=2)

        # -- Panel 3: Distance vs Time (fixed energy & pitch) --
        z3 = np.log10(phi_t[energy_index, pitch_index, :, :])
        axes[2].pcolormesh(data.tr, data.soft / 1e5, z3, cmap=cmap, shading="auto")
        axes[2].set_xlabel("Time [s]", fontsize=14)
        axes[2].set_ylabel("Distance [km]", fontsize=14)
        axes[2].set_title(
            f"E = {data.ekev[energy_index]:.0f} keV — "
            rf"$\alpha$ = {data.pad[pitch_index]:.0f}°",
            fontsize=14,
        )
        axes[2].axhline(y=data.soft[space_index] / 1e5, color="b", ls=":", lw=2)

    for a in axes:
        _style_ax(a)
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 3D surface panels
# ---------------------------------------------------------------------------
def plot_3d_surfaces(
    data: FPOutput,
    energy_index: int = 12,
    pitch_index: int = 30,
    space_index: int = 30,
    cmap: str = "YlGnBu",
    save: str | None = None,
) -> Figure:
    """Three-panel 3-D surface: energy vs time, pitch vs time, distance vs time.

    Parameters
    ----------
    data : FPOutput
        Loaded simulation data.
    energy_index, pitch_index, space_index : int
        Fixed-axis indices for the cross-sections shown.
    cmap : str
        Matplotlib colour-map name.
    save : str | None
        If given, save the figure to this path.

    Returns
    -------
    Figure
    """
    fig, axes = plt.subplots(
        1, 3, subplot_kw={"projection": "3d"}, figsize=(22, 6),
    )

    phi_t = data.phi.transpose((3, 2, 1, 0))

    # Panel 1: Energy vs Time
    T, E = np.meshgrid(data.tr, data.ekev)
    axes[0].plot_surface(
        T, E, phi_t[:, pitch_index, space_index, :],
        cmap=cmap, linewidth=0, antialiased=False,
    )
    axes[0].set_xlabel("Time [s]", fontsize=12, labelpad=8)
    axes[0].set_ylabel("Energy [keV]", fontsize=12, labelpad=10)
    axes[0].set_zlabel("Density", fontsize=12, labelpad=12)
    axes[0].set_title(
        rf"$\alpha$={data.pad[pitch_index]:.0f}° — s={data.soft[space_index]/1e5:.0f} km",
        fontsize=13,
    )
    axes[0].zaxis.set_major_locator(LinearLocator(8))

    # Panel 2: Pitch vs Time
    T, PA = np.meshgrid(data.tr, data.pad)
    axes[1].plot_surface(
        T, PA, phi_t[energy_index, :, space_index, :],
        cmap=cmap, linewidth=0, antialiased=False,
    )
    axes[1].set_xlabel("Time [s]", fontsize=12, labelpad=8)
    axes[1].set_ylabel("Pitch Angle [deg]", fontsize=12, labelpad=10)
    axes[1].set_zlabel("Density", fontsize=12, labelpad=12)
    axes[1].set_title(
        f"E={data.ekev[energy_index]:.0f} keV — s={data.soft[space_index]/1e5:.0f} km",
        fontsize=13,
    )
    axes[1].zaxis.set_major_locator(LinearLocator(8))

    # Panel 3: Distance vs Time
    T, S = np.meshgrid(data.tr, data.soft / 1e5)
    axes[2].plot_surface(
        T, S, phi_t[energy_index, pitch_index, :, :],
        cmap=cmap, linewidth=0, antialiased=False,
    )
    axes[2].set_xlabel("Time [s]", fontsize=12, labelpad=8)
    axes[2].set_ylabel("Distance [km]", fontsize=12, labelpad=10)
    axes[2].set_zlabel("Density", fontsize=12, labelpad=12)
    axes[2].set_title(
        rf"E={data.ekev[energy_index]:.0f} keV — $\alpha$={data.pad[pitch_index]:.0f}°",
        fontsize=13,
    )
    axes[2].zaxis.set_major_locator(LinearLocator(8))

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=200, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# CLI entry point — generate all standard plots for a given output file
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    from pathlib import Path

    from tdfp_reader import read_fkrplk_output

    parser = argparse.ArgumentParser(
        description="Generate standard TDFP distribution plots.",
    )
    parser.add_argument("datafile", help="Path to fkrplk output file")
    parser.add_argument(
        "-o", "--outdir", default=".", help="Directory for saved figures",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {args.datafile} ...")
    data = read_fkrplk_output(args.datafile)
    print(
        f"  Grid: {data.newr+1} energy × {2*data.nmu+1} pitch × "
        f"{data.ntau+1} space × {data.ntr+1} time"
    )

    print("Plotting energy spectrum ...")
    plot_energy_spectrum(data, save=str(outdir / "energy_spectrum.png"))

    print("Plotting pitch-angle distribution ...")
    plot_pitch_angle(data, save=str(outdir / "pitch_angle.png"))

    print("Plotting spatial profile ...")
    plot_spatial_profile(data, save=str(outdir / "spatial_profile.png"))

    print("Plotting time evolution ...")
    plot_time_evolution(data, save=str(outdir / "time_evolution.png"))

    print("Plotting 2D panels ...")
    plot_2d_panels(data, save=str(outdir / "contour_panels.png"))

    print("Plotting 3D surfaces ...")
    plot_3d_surfaces(data, save=str(outdir / "surface_3d.png"))

    print(f"All figures saved to {outdir}/")
