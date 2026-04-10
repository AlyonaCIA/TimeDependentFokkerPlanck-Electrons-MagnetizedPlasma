"""Dynamic HXR spectrogram — Soft-Hard-Soft spectral evolution.

Computes the time-resolved Hard X-ray photon spectrum via thin-target
bremsstrahlung and displays it as a 2-D dynamic spectrogram
:math:`J(\\varepsilon, t)`:

.. math::

    J(\\varepsilon, t)
    \\;\\propto\\;
    \\frac{n_{\\mathrm{bg}}}{\\varepsilon}
    \\int_{\\varepsilon}^{E_{\\max}}
    \\frac{1}{E}
    \\int_0^L \\int_{-1}^{1} f(E, \\mu, s, t)\\, d\\mu\\, ds\\, dE

The :math:`1/(\\varepsilon E)` weighting is the non-relativistic Bethe-Heitler
(Kramers) bremsstrahlung cross-section.  The integral runs over all electron
energies *above* the photon energy :math:`\\varepsilon`, because only electrons
with :math:`E \\ge \\varepsilon` can radiate at that photon energy.

The classic **Soft-Hard-Soft (SHS)** pattern appears as the spectral peak
(hardest emission) coinciding with the impulsive injection phase, then
softening as collisional cooling steepens the electron distribution.

Usage
-----
Real Fortran output::

    uv run python scripts/python/plot_dynamic_spectrogram.py \\
        --data experiments/conf_original/fkrplk.test

Save to file::

    uv run python scripts/python/plot_dynamic_spectrogram.py \\
        --data experiments/conf_original/fkrplk.test \\
        --save output/figures/09_dynamic_spectrogram.png

With the high-resolution conf4 data (120 time snapshots)::

    uv run python scripts/python/plot_dynamic_spectrogram.py \\
        --data experiments/conf4/fkrplk.conf4 \\
        --save output/figures/09_dynamic_spectrogram.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure


# ---------------------------------------------------------------------------
# Physics: compute photon spectrum J(ε, t)
# ---------------------------------------------------------------------------

def compute_photon_spectrum(
    phi_t: np.ndarray,
    ekev: np.ndarray,
    amu: np.ndarray,
    soft: np.ndarray,
    E_photon: np.ndarray,
    n_bg: float = 1.0e11,
) -> np.ndarray:
    r"""Compute thin-target HXR photon spectrum for one time snapshot.

    .. math::

        J(\varepsilon) \propto \frac{n_{\rm bg}}{\varepsilon}
        \int_{\varepsilon}^{E_{\max}} \frac{1}{E}
        \int_0^L \int_{-1}^{1} f(E,\mu,s)\, d\mu\, ds\, dE

    Parameters
    ----------
    phi_t : ndarray, shape (n_space, n_pitch, n_energy)
        Distribution function at one time snapshot.
    ekev : ndarray, shape (n_energy,)
        Electron energy grid in keV.
    amu : ndarray, shape (n_pitch,)
        Pitch-angle cosine grid.
    soft : ndarray, shape (n_space,)
        Spatial position along the loop in cm.
    E_photon : ndarray, shape (n_photon,)
        Photon energy grid in keV at which to evaluate the spectrum.
    n_bg : float
        Background electron density in cm⁻³.

    Returns
    -------
    J : ndarray, shape (n_photon,)
        Photon flux in arbitrary units at each photon energy.
    """
    phi_clean = np.where(np.isfinite(phi_t), phi_t, 0.0)

    # F(E) = ∫∫ f(E,μ,s) dμ ds  — the spatially integrated electron spectrum
    integrated_mu = np.trapezoid(phi_clean, x=amu, axis=1)   # (space, energy)
    F_E = np.trapezoid(integrated_mu, x=soft, axis=0)        # (energy,)
    F_E = np.where(np.isfinite(F_E), F_E, 0.0)

    # Integrand g(E) = F(E) / E  (Kramers cross-section)
    inv_E = np.where(ekev > 0, 1.0 / ekev, 0.0)
    g_E = F_E * inv_E  # (n_energy,)

    # For each photon energy ε, integrate g(E) from ε to E_max
    J = np.zeros(E_photon.shape[0])
    for i, eps in enumerate(E_photon):
        mask = ekev >= eps
        if np.count_nonzero(mask) < 2:
            J[i] = 0.0
            continue
        J[i] = np.trapezoid(g_E[mask], x=ekev[mask])

    # Multiply by n_bg / ε
    inv_eps = np.where(E_photon > 0, 1.0 / E_photon, 0.0)
    J = n_bg * J * inv_eps

    return np.abs(J)


def compute_spectrogram(
    phi: np.ndarray,
    ekev: np.ndarray,
    amu: np.ndarray,
    soft: np.ndarray,
    E_photon: np.ndarray,
    n_bg: float = 1.0e11,
    t_indices: np.ndarray | None = None,
) -> np.ndarray:
    """Compute the 2-D spectrogram J(ε, t) for selected time snapshots.

    Parameters
    ----------
    phi : ndarray, shape (n_time, n_space, n_pitch, n_energy)
        Full 4-D distribution function.
    ekev, amu, soft : ndarray
        Grid arrays.
    E_photon : ndarray, shape (n_photon,)
        Photon energy grid in keV.
    n_bg : float
        Background density.
    t_indices : ndarray or None
        Time indices to include.  Defaults to all.

    Returns
    -------
    spectrogram : ndarray, shape (len(t_indices), n_photon)
        Photon flux at each (time, photon energy).
    """
    if t_indices is None:
        t_indices = np.arange(phi.shape[0])

    spectrogram = np.zeros((len(t_indices), E_photon.shape[0]))
    for row, ti in enumerate(t_indices):
        spectrogram[row] = compute_photon_spectrum(
            phi[ti], ekev, amu, soft, E_photon, n_bg,
        )
    return spectrogram


def _mock_spectrogram() -> (
    tuple[np.ndarray, np.ndarray, np.ndarray]
):
    """Synthetic Soft-Hard-Soft spectrogram for testing.

    Simulates a triangular injection pulse where the spectral index
    hardens during the rise and softens during the decay.
    """
    n_time = 100
    t_array = np.linspace(0, 10, n_time)  # seconds

    n_photon = 80
    E_photon = np.logspace(np.log10(5), np.log10(500), n_photon)  # keV

    spectrogram = np.zeros((n_time, n_photon))
    t_peak = 3.0  # peak of impulsive phase (seconds)
    dur = 4.0     # total injection duration

    for it, t in enumerate(t_array):
        # Triangular intensity envelope
        if t < t_peak:
            amplitude = t / t_peak
        elif t < t_peak + dur:
            amplitude = max(0, 1 - (t - t_peak) / dur)
        else:
            amplitude = 1e-6

        # Spectral index: hardens toward peak, then softens
        # SHS pattern: δ goes from 6 → 3 → 6
        if t < t_peak:
            delta = 6.0 - 3.0 * (t / t_peak)
        elif t < t_peak + dur:
            delta = 3.0 + 3.0 * ((t - t_peak) / dur)
        else:
            delta = 6.0

        # Power-law photon spectrum: J(ε) = A * ε^(-δ)
        spectrogram[it] = amplitude * 1e6 * E_photon ** (-delta)

    return t_array, E_photon, spectrogram


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_dynamic_spectrogram(
    t_array: np.ndarray,
    E_photon: np.ndarray,
    spectrogram: np.ndarray,
    *,
    cmap: str = "inferno",
    save: str | None = None,
    dpi: int = 300,
    figsize: tuple[float, float] = (12, 6),
    show_spectral_index: bool = True,
) -> Figure:
    r"""Plot the 2-D dynamic spectrogram as a pcolormesh.

    Parameters
    ----------
    t_array : ndarray, shape (n_time,)
        Time axis in seconds.
    E_photon : ndarray, shape (n_photon,)
        Photon energy axis in keV.
    spectrogram : ndarray, shape (n_time, n_photon)
        Photon flux :math:`J(\varepsilon, t)`.
    cmap : str
        Colourmap name (default: ``'inferno'``).
    save : str | None
        Save figure to this path.
    dpi : int
        Figure resolution.
    figsize : tuple
        Figure size in inches.
    show_spectral_index : bool
        If True, overlay the fitted spectral index γ(t) on a twin axis.

    Returns
    -------
    Figure
    """
    # Floor to avoid log(0)
    spec_pos = np.where(spectrogram > 0, spectrogram, np.nan)
    finite = spec_pos[np.isfinite(spec_pos)]
    if finite.size == 0:
        print("Warning: spectrogram is all zeros. Nothing to plot.")
        return plt.figure()

    vmin = float(np.nanmin(finite))
    vmax = float(np.nanmax(finite))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # pcolormesh: T on x-axis, E_photon on y-axis
    pcm = ax.pcolormesh(
        t_array,
        E_photon,
        spec_pos.T,
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="gouraud",
        rasterized=True,
    )

    ax.set_yscale("log")
    ax.set_xlabel(r"Time  $t$  (s)", fontsize=14)
    ax.set_ylabel(r"Photon Energy  $\varepsilon$  (keV)", fontsize=14)
    ax.set_title(
        r"Dynamic HXR Spectrogram — Soft-Hard-Soft Evolution",
        fontsize=16,
        fontweight="bold",
    )
    ax.tick_params(labelsize=12)

    cbar = fig.colorbar(pcm, ax=ax, pad=0.02, aspect=30)
    cbar.set_label(
        r"$\log_{10}\; J(\varepsilon,\,t)$"
        r"  (photons / cm$^2$ / s / keV)",
        fontsize=13,
    )
    cbar.ax.tick_params(labelsize=11)

    # ── Spectral index overlay ──────────────────────────────────────────
    if show_spectral_index:
        gamma = _fit_spectral_index(t_array, E_photon, spectrogram)
        if gamma is not None:
            ax2 = ax.twinx()
            ax2.plot(t_array, gamma, color="cyan", lw=2.0, ls="--",
                     alpha=0.9, label=r"$\gamma(t)$")
            ax2.set_ylabel(
                r"Spectral Index  $\gamma(t)$",
                fontsize=13,
                color="cyan",
            )
            ax2.tick_params(axis="y", labelcolor="cyan", labelsize=11)
            ax2.legend(loc="upper right", fontsize=11, framealpha=0.7)

    # Grid
    ax.grid(True, which="major", ls=":", lw=0.3, alpha=0.4, color="white")

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight", facecolor="white")
        print(f"Saved: {save}")
        plt.close(fig)

    return fig


def _fit_spectral_index(
    t_array: np.ndarray,
    E_photon: np.ndarray,
    spectrogram: np.ndarray,
    E_fit_range: tuple[float, float] = (20.0, 200.0),
) -> np.ndarray | None:
    r"""Fit a power-law spectral index γ(t) in a photon energy band.

    For each time step, fit :math:`J(\varepsilon) \propto \varepsilon^{-\gamma}`
    in the range ``E_fit_range`` using a least-squares linear fit in
    log-log space.

    Returns
    -------
    gamma : ndarray, shape (n_time,) or None
        Best-fit spectral index at each time step.
    """
    mask = (E_photon >= E_fit_range[0]) & (E_photon <= E_fit_range[1])
    if np.count_nonzero(mask) < 3:
        return None

    log_E = np.log10(E_photon[mask])
    gamma = np.full(spectrogram.shape[0], np.nan)

    for it in range(spectrogram.shape[0]):
        vals = spectrogram[it, mask]
        pos = vals > 0
        if np.count_nonzero(pos) < 3:
            continue
        log_J = np.log10(vals[pos])
        # Linear fit: log_J = -gamma * log_E + c
        coeffs = np.polyfit(log_E[pos], log_J, 1)
        gamma[it] = -coeffs[0]

    # If all NaN, return None
    if np.all(np.isnan(gamma)):
        return None
    return gamma


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Dynamic HXR spectrogram — Soft-Hard-Soft spectral evolution."
        ),
    )
    p.add_argument(
        "--data", type=str, default=None,
        help="Path to a fkrplk output file.  If omitted, uses mock data.",
    )
    p.add_argument(
        "--save", type=str, default=None,
        help="Save figure to this path.",
    )
    p.add_argument(
        "--cmap", type=str, default="inferno",
        help="Colourmap (default: inferno).",
    )
    p.add_argument(
        "--dpi", type=int, default=300,
        help="Figure DPI (default: 300).",
    )
    p.add_argument(
        "--no-gamma", action="store_true",
        help="Do not overlay the spectral index γ(t).",
    )
    p.add_argument(
        "--time-indices", type=str, default=None,
        help="Comma-separated time indices (e.g. '0,1,2,3'). "
             "Defaults to all meaningful snapshots.",
    )
    p.add_argument(
        "--n-photon", type=int, default=60,
        help="Number of photon energy bins (default: 60).",
    )
    p.add_argument(
        "--emin", type=float, default=None,
        help="Minimum photon energy in keV (default: electron grid min).",
    )
    p.add_argument(
        "--emax", type=float, default=None,
        help="Maximum photon energy in keV (default: electron grid max).",
    )
    return p


def main() -> None:
    """Entry point for the CLI."""
    args = _build_parser().parse_args()

    if args.data is not None:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from tdfp_reader import read_fkrplk_output

        data = read_fkrplk_output(args.data)
        n_times = data.phi.shape[0]

        # Detect last meaningful time step
        if args.time_indices is not None:
            t_indices = np.array(
                [int(x) for x in args.time_indices.split(",")]
            )
        else:
            last_good = n_times - 1
            for ti in range(n_times - 1, -1, -1):
                p = data.phi[ti]
                finite = p[np.isfinite(p)]
                if finite.size > 0 and np.max(np.abs(finite)) > 1e-20:
                    last_good = ti
                    break
            t_indices = np.arange(0, last_good + 1)

        # Photon energy grid
        emin = args.emin if args.emin is not None else float(data.ekev[data.ekev > 0].min())
        emax = args.emax if args.emax is not None else float(data.ekev.max())
        E_photon = np.logspace(
            np.log10(emin), np.log10(emax), args.n_photon,
        )

        # Build spectrogram
        spectrogram = compute_spectrogram(
            data.phi, data.ekev, data.amu, data.soft,
            E_photon, n_bg=data.rhomax, t_indices=t_indices,
        )

        # Time axis: use report times if available
        t_array = data.tr[t_indices]

    else:
        # Mock data
        t_array, E_photon, spectrogram = _mock_spectrogram()

    plot_dynamic_spectrogram(
        t_array,
        E_photon,
        spectrogram,
        cmap=args.cmap,
        save=args.save,
        dpi=args.dpi,
        show_spectral_index=not args.no_gamma,
    )

    if args.save is None:
        plt.show()


if __name__ == "__main__":
    main()
