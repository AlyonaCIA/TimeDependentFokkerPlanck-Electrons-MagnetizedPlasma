"""3-D volumetric visualisation of electron density along a coronal loop.

Maps the 1-D spatial coordinate *s* produced by the Fokker-Planck solver
onto a semi-circular coronal-loop geometry and renders it as a glowing
volumetric tube coloured by electron density n_e(s).

The visual style is inspired by the Bifrost/Helen 3-D solar-atmosphere
simulations from the Rosseland Centre for Solar Physics (UiO, 2021).

Usage
-----
Interactive window::

    uv run python scripts/python/plot_3d_coronal_loop.py

Save a high-resolution PNG (headless)::

    uv run python scripts/python/plot_3d_coronal_loop.py --save coronal_loop.png

Load real simulation data::

    uv run python scripts/python/plot_3d_coronal_loop.py \\
        --data experiments/conf_original/fkrplk.test --time-index 3
"""

from __future__ import annotations

# -- Standard library --------------------------------------------------------
import argparse
import sys
from pathlib import Path

# -- Third-party -------------------------------------------------------------
import numpy as np
import pyvista as pv


# ---------------------------------------------------------------------------
# Physics helpers
# ---------------------------------------------------------------------------

def _mock_density(n_points: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """Return a synthetic electron-density profile n_e(s).

    The profile peaks at both footpoints (s = 0 and s = L) and has a
    minimum at the loop apex, which is the canonical behaviour for
    coronal-loop electron distributions.

    Returns
    -------
    s : ndarray, shape (n_points,)
        Arc-length coordinate in Mm.
    n_e : ndarray, shape (n_points,)
        Electron number density in cm⁻³.
    """
    loop_length = 50.0  # Mm (half-circle arc-length = pi * R)
    s = np.linspace(0.0, loop_length, n_points)
    # Normalised position [0, 1]
    t = s / loop_length
    # Footpoint-enhanced density: sum of two Gaussians at each end
    n_e = (
        2.0e10 * np.exp(-((t - 0.0) ** 2) / (2 * 0.04**2))
        + 2.0e10 * np.exp(-((t - 1.0) ** 2) / (2 * 0.04**2))
        + 5.0e8  # tenuous coronal background
        + 3.0e9 * np.exp(-((t - 0.5) ** 2) / (2 * 0.15**2))  # mild apex heating
    )
    return s, n_e


def _density_from_simulation(filepath: str | Path, time_index: int = -1) -> tuple[np.ndarray, np.ndarray]:
    """Integrate f(E, mu, s) from Fortran output to get n_e(s).

    Parameters
    ----------
    filepath : str | Path
        Path to a ``fkrplk`` ASCII output file.
    time_index : int
        Which time snapshot to use (-1 = last available).

    Returns
    -------
    s : ndarray
        Spatial coordinate in Mm.
    n_e : ndarray
        Electron number density in cm⁻³.
    """
    # Import here to avoid hard dependency when using mock data
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tdfp_reader import read_fkrplk_output

    data = read_fkrplk_output(filepath)
    # phi has shape (time, space, pitch, energy)
    phi_t = data.phi[time_index]  # (space, pitch, energy)
    ekev = data.ekev              # (energy,)

    # n_e(s) = integral over energy and pitch angle of f(E, mu, s)
    # Trapezoidal integration in energy, then sum over pitch angles
    # f is in units that when integrated give density
    n_e = np.trapz(
        np.trapz(phi_t, x=ekev, axis=-1),  # integrate over energy
        x=data.amu,                          # integrate over pitch angle
        axis=-1,
    )
    n_e = np.abs(n_e)  # ensure positive

    # Convert spatial coordinate from cm to Mm (1 Mm = 1e8 cm)
    s = data.soft / 1.0e8

    return s, n_e


# ---------------------------------------------------------------------------
# Geometry: map 1-D arc-length onto 3-D semi-circle
# ---------------------------------------------------------------------------

def _loop_geometry(
    s: np.ndarray,
    radius_mm: float | None = None,
) -> np.ndarray:
    """Map arc-length *s* onto a 3-D semi-circular coronal loop.

    The loop lies in the x-z plane with footpoints on the z = 0 surface.

    Parameters
    ----------
    s : ndarray, shape (N,)
        Arc-length coordinate (same units as *radius_mm*).
    radius_mm : float | None
        Loop radius.  Defaults to ``s[-1] / pi``.

    Returns
    -------
    points : ndarray, shape (N, 3)
        Cartesian (x, y, z) coordinates.
    """
    loop_length = s[-1]
    if radius_mm is None:
        radius_mm = loop_length / np.pi  # so that arc = pi * R = L

    theta = s / radius_mm  # angle along the loop [0, pi]
    x = radius_mm * np.cos(theta)
    y = np.zeros_like(s)
    z = radius_mm * np.sin(theta)
    return np.column_stack((x, y, z))


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_coronal_loop(
    s: np.ndarray,
    n_e: np.ndarray,
    *,
    save: str | None = None,
    cmap: str = "inferno",
    tube_radius: float = 0.6,
    window_size: tuple[int, int] = (1920, 1080),
) -> None:
    """Render the coronal loop as a glowing volumetric tube.

    Parameters
    ----------
    s, n_e : ndarray
        Arc-length and electron-density arrays from the mock or real data.
    save : str | None
        If given, save a screenshot to this path (headless).
    cmap : str
        Matplotlib / PyVista colourmap name.
    tube_radius : float
        Tube radius in the same spatial units as *s*.
    window_size : tuple[int, int]
        Render window (width, height) in pixels.
    """
    points = _loop_geometry(s)
    R = s[-1] / np.pi  # loop radius for camera placement

    # -- Build the tube mesh --------------------------------------------------
    n_spline = len(points) * 3
    spline = pv.Spline(points, n_points=n_spline)

    # Interpolate n_e onto the denser spline points via arc-length
    arc = spline.point_data["arc_length"]
    n_e_interp = np.interp(arc, np.linspace(arc.min(), arc.max(), len(n_e)), n_e)

    # Use log10 density for colour mapping (compresses dynamic range so
    # the coronal mid-section is visible against the dark background)
    log_ne = np.log10(np.clip(n_e_interp, 1.0, None))
    spline.point_data["log_ne"] = log_ne

    n_min, n_max = float(log_ne.min()), float(log_ne.max())
    clim = [n_min, n_max]

    tube = spline.tube(radius=tube_radius, n_sides=48)

    # Opacity transfer function: 20 samples from low → high density
    # Ensure the coronal (low-density) part is still clearly visible
    opacity_table = [0.45 + 0.55 * (i / 19.0) ** 0.5 for i in range(20)]

    # Outer glow shell (larger, semi-transparent, emissive)
    glow = spline.tube(radius=tube_radius * 2.5, n_sides=48)

    # Chromospheric footpoint spheres
    fp1 = pv.Sphere(radius=tube_radius * 2.8, center=points[0])
    fp2 = pv.Sphere(radius=tube_radius * 2.8, center=points[-1])

    # Solar surface
    x_extent = float(np.abs(points[:, 0]).max()) * 1.5
    surface = pv.Plane(
        center=(0.0, 0.0, -0.15),
        direction=(0.0, 0.0, 1.0),
        i_size=x_extent * 2.2,
        j_size=x_extent * 1.4,
        i_resolution=4,
        j_resolution=4,
    )

    # -- Plotter --------------------------------------------------------------
    off_screen = save is not None
    pl = pv.Plotter(
        window_size=list(window_size),
        lighting="none",
        off_screen=off_screen,
    )
    pl.set_background("black")

    # Three-point lighting
    pl.add_light(pv.Light(
        position=(0, 0, 50), color=(0.25, 0.27, 0.45),
        intensity=0.7, light_type="scene light",
    ))
    pl.add_light(pv.Light(
        position=(30, 25, 40), focal_point=(0, 0, R * 0.5),
        color=(1.0, 0.95, 0.85), intensity=1.5, light_type="scene light",
    ))
    pl.add_light(pv.Light(
        position=(-25, -20, 35), focal_point=(0, 0, R * 0.5),
        color=(0.4, 0.5, 0.9), intensity=0.9, light_type="scene light",
    ))

    # -- Add meshes -----------------------------------------------------------
    # Outer glow (emissive, semi-transparent)
    pl.add_mesh(
        glow, scalars="log_ne", cmap=cmap, clim=clim,
        opacity=0.15, smooth_shading=True,
        ambient=0.9, diffuse=0.1, show_scalar_bar=False,
    )

    # Main tube — show_scalar_bar=True so we get a correctly-ranged bar
    sbar_args = dict(
        title="log₁₀ n_e  [cm⁻³]", n_labels=5, shadow=True, italic=True,
        fmt="%.1f", color="white",
        title_font_size=16, label_font_size=12,
        position_x=0.72, position_y=0.05, width=0.22, height=0.06,
    )
    pl.add_mesh(
        tube, scalars="log_ne", cmap=cmap, clim=clim,
        opacity=opacity_table, smooth_shading=True,
        ambient=0.5, diffuse=0.5, specular=0.9, specular_power=50,
        scalar_bar_args=sbar_args,
    )

    # Footpoint spheres
    for fp in (fp1, fp2):
        pl.add_mesh(
            fp, color="#FF6030", opacity=0.92, smooth_shading=True,
            ambient=0.7, diffuse=0.3, specular=0.5,
        )

    # Solar chromosphere surface
    pl.add_mesh(
        surface, color="#1A0800", opacity=0.6,
        smooth_shading=True, ambient=0.4,
    )

    # Title text
    pl.add_text(
        "Coronal Loop — Electron Density",
        position="upper_left",
        font_size=14,
        color="white",
        shadow=True,
    )

    # Camera: slightly elevated, looking at apex
    pl.camera_position = [
        (0.0, -R * 3.2, R * 1.4),
        (0.0, 0.0, R * 0.55),
        (0.0, 0.0, 1.0),
    ]
    pl.enable_anti_aliasing("ssaa")

    if save:
        pl.show(screenshot=save)
        print(f"Saved: {save}")
    else:
        pl.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="3-D coronal-loop electron-density visualisation.",
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help="Path to a fkrplk output file.  If omitted, uses mock data.",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        default=-1,
        help="Time-snapshot index to plot (default: last).",
    )
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        help="Save a screenshot to this path instead of showing a window.",
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="inferno",
        help="Colourmap name (default: inferno).",
    )
    parser.add_argument(
        "--tube-radius",
        type=float,
        default=0.6,
        help="Tube radius in Mm (default: 0.6).",
    )
    parser.add_argument(
        "--resolution",
        type=str,
        default="2560x1440",
        help="Window size as WxH (default: 2560x1440).",
    )
    return parser


def main() -> None:
    """Entry point for the CLI."""
    args = _build_parser().parse_args()

    # Parse resolution
    w, h = (int(x) for x in args.resolution.split("x"))

    # Load or generate density profile
    if args.data is not None:
        s, n_e = _density_from_simulation(args.data, time_index=args.time_index)
    else:
        s, n_e = _mock_density()

    render_coronal_loop(
        s,
        n_e,
        save=args.save,
        cmap=args.cmap,
        tube_radius=args.tube_radius,
        window_size=(w, h),
    )


if __name__ == "__main__":
    main()
