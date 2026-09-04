"""Notebook-friendly workflows for bathymetry analysis."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
import importlib
import subprocess
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from scipy.interpolate import splprep, splev

from utils.date_utils import interval_highlight_mask
from bathy_formatter import get_named_colormap
from .bathy import load_clrmap_file
from .plot_style import apply_axes_font, apply_global_style, get_plot_font
from .stratigraphy import (
    StratigraphyConfig,
    Transect,
    compute_stratigraphy,
    datenum_to_datetime64,
    compute_deposit_elev_1d,
    extract_transect_cube,
    load_bathy_cube,
    load_transects_from_shapefiles,
    plot_stacked_stratigraphy_section,
    plot_stratigraphy_stack,
)


def _transect_hits_domain(bathy_cube, transect: Transect) -> bool:
    """Return True if any transect point intersects the bathymetry domain hull."""
    z_last = bathy_cube.z[:, :, -1]
    mask = np.isfinite(z_last)
    if not np.any(mask):
        return False

    x_valid = bathy_cube.x[mask].astype(float)
    y_valid = bathy_cube.y[mask].astype(float)
    pts = np.column_stack([x_valid, y_valid])
    if pts.shape[0] < 3:
        return False

    from scipy.spatial import ConvexHull

    hull = ConvexHull(pts)
    hull_pts = pts[hull.vertices]
    hull_path = MplPath(hull_pts)
    transect_pts = np.column_stack([transect.x, transect.y])
    return bool(np.any(hull_path.contains_points(transect_pts)))


def _unique_xy(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Remove duplicate xy rows while preserving order."""
    coords = np.column_stack([x, y])
    unique_coords = np.unique(coords, axis=0)
    if unique_coords.shape[0] != coords.shape[0]:
        # Stable unique to preserve order.
        _, idx = np.unique(coords, axis=0, return_index=True)
        coords = coords[np.sort(idx)]
    return coords[:, 0], coords[:, 1]


def _smooth_transect_xy(
    x: np.ndarray,
    y: np.ndarray,
    n_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a spline through a polyline and resample it with n_points."""
    if n_points < 2:
        raise ValueError("n_points must be >= 2 for spline resampling")
    if x.size < 4:
        raise ValueError("need at least 4 points for spline smoothing")
    pts = np.column_stack([x, y])
    seg_len = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    s = np.concatenate([[0.0], np.cumsum(seg_len)])
    if s[-1] == 0.0:
        return x, y
    t = s / s[-1]
    tck, _ = splprep([x, y], u=t, s=0.0)
    t_new = np.linspace(0.0, 1.0, n_points)
    x_new, y_new = splev(t_new, tck)
    return np.asarray(x_new), np.asarray(y_new)


def _load_profile_overlay(
    csv_source: str | Path,
    transect_label: str,
    x_column: str | None = None,
    y_column: str | None = None,
    align_start: bool = False,
) -> np.ndarray | None:
    """Load profile-space points from one CSV or a label-matched CSV directory."""
    source = Path(csv_source)
    if source.is_dir():
        label_prefix = f"{transect_label}-{transect_label}".casefold()
        matches = [path for path in sorted(source.glob("*.csv")) if path.stem.casefold().startswith(label_prefix)]
        if not matches:
            return None
        source = matches[0]
    elif not source.is_file():
        raise FileNotFoundError(f"Overlay CSV source does not exist: {source}")

    if (x_column is None) != (y_column is None):
        raise ValueError("overlay_x_column and overlay_y_column must be provided together")
    try:
        frame = pd.read_csv(source, header=None if x_column is None else "infer")
    except pd.errors.EmptyDataError:
        print(f"WARNING: overlay CSV {source} is empty and was skipped.")
        return None

    if x_column is None:
        if frame.shape[1] < 2:
            raise ValueError(f"Overlay CSV must contain at least two columns: {source}")
        x_values = pd.to_numeric(frame.iloc[:, 0], errors="coerce")
        y_values = pd.to_numeric(frame.iloc[:, 1], errors="coerce")
    else:
        missing = [column for column in (x_column, y_column) if column not in frame.columns]
        if missing:
            raise ValueError(f"Overlay CSV {source} is missing column(s): {', '.join(missing)}")
        x_values = pd.to_numeric(frame[x_column], errors="coerce")
        y_values = pd.to_numeric(frame[y_column], errors="coerce")

    points = np.column_stack([x_values.to_numpy(dtype=float), y_values.to_numpy(dtype=float)])
    points = points[np.all(np.isfinite(points), axis=1)]
    if align_start and points.size:
        points[:, 0] -= np.min(points[:, 0])
    return points


def _transect_output_path(directory: Path, filename: str, overlay_enabled: bool) -> Path:
    """Build a transect output path with the optional overlay suffix."""
    path = directory / filename
    if overlay_enabled:
        return path.with_name(f"{path.stem}_overlay{path.suffix}")
    return path

DEFAULT_REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "netCDF4": "netCDF4",
    "pyproj": "pyproj",
    "geopandas": "geopandas",
    "joblib": "joblib",
    "imageio": "imageio",
}

DEFAULT_METRICS = {
    "erosion_m3": "Erosion",
    "accretion_m3": "Accretion",
    "gross_m3": "Gross",
    "net_m3": "Net",
}

DEFAULT_X_VARS = {
    "cum_wave_power_MWh_m": "Total wave power [MWh/m]",
    "cum_wave_power_above_MWh_m": "Wave power (Hs >= 2.0 m) [MWh/m]",
    "cum_wave_power_below_MWh_m": "Wave power (Hs < 2.0 m) [MWh/m]",
    "period_days": "Interval duration [days]",
}


def _label_for_index(index: int) -> str:
    """Convert a 0-based index into A, B, ..., Z, AA, BB, ... labels."""
    repeat = index // 26 + 1
    letter = chr(ord("A") + (index % 26))
    return letter * repeat


def _normalize_highlight_dates(highlight_dates) -> list[str | np.datetime64]:
    """Normalize highlight dates into a list."""
    if highlight_dates is None:
        return []
    if isinstance(highlight_dates, (list, tuple, np.ndarray)):
        return [value for value in highlight_dates if value is not None]
    return [highlight_dates]


def _section_figsize(
    x_range: float,
    y_range: float,
    max_transect_length: float | None,
    max_depth_range: float | None,
    max_fig_width_cm: float,
    max_fig_height_cm: float,
    scale_factor: float = 2.0,
) -> tuple[float, float]:
    """Compute a consistent figure size for cross-section plots."""
    if max_transect_length and max_depth_range and max_transect_length > 0 and max_depth_range > 0:
        max_fig_width_in = max_fig_width_cm / 2.54
        max_fig_height_in = max_fig_height_cm / 2.54
        fig_width = (x_range / max_transect_length) * max_fig_width_in * scale_factor
        fig_height = (y_range / max_depth_range) * max_fig_height_in * scale_factor
        fig_width = max(fig_width, 3.0)
        fig_height = max(fig_height, 2.0)
        return fig_width, fig_height
    return 14.0, 4.0


def _write_stratigraphy_gif(
    slice_cube,
    out_path: Path,
    section_title: str,
    mhw: float,
    mlw: float,
    plot_ylim: tuple[float, float] | None,
    highlight_dates: list[str | np.datetime64] | None,
    max_transect_length: float | None,
    max_depth_range: float | None,
    max_fig_width_cm: float,
    max_fig_height_cm: float,
    fps: int = 6,
    frame_stride: int = 1,
    overlay_points: np.ndarray | None = None,
    overlay_color: str = "#00E5FF",
    overlay_edge_color: str = "black",
    overlay_marker_size: float = 8.0,
) -> None:
    """Create an animated GIF showing stacked stratigraphy over time."""
    try:
        import imageio.v2 as imageio
    except ImportError as exc:
        raise ImportError("imageio is required for GIF export") from exc

    apply_global_style()
    font = get_plot_font()

    x_m = slice_cube.x[0, :]
    z_stack = slice_cube.z[0, :, :].T
    deposit_elev_full = compute_deposit_elev_1d(z_stack)
    nt = z_stack.shape[0]
    colors = plt.cm.viridis(np.linspace(0.15, 0.95, nt))

    min_elev = float(np.nanmin([np.nanmin(z_stack), np.nanmin(deposit_elev_full)]))
    max_elev = float(np.nanmax([np.nanmax(z_stack), np.nanmax(deposit_elev_full)]))
    y_min_plot, y_max_plot = (plot_ylim if plot_ylim is not None else (min_elev, max_elev))
    x_min_plot = float(np.nanmin(x_m)) if len(x_m) else 0.0
    x_max_plot = float(np.nanmax(x_m)) if len(x_m) else 0.0
    if overlay_points is not None and overlay_points.size:
        x_min_plot = min(x_min_plot, float(np.nanmin(overlay_points[:, 0])))
        x_max_plot = max(x_max_plot, float(np.nanmax(overlay_points[:, 0])))
        y_min_plot = min(y_min_plot, float(np.nanmin(overlay_points[:, 1])))
        y_max_plot = max(y_max_plot, float(np.nanmax(overlay_points[:, 1])))
    base = min_elev - 0.01 * abs(min_elev)
    # Clip to the baseline to avoid extra whitespace below the section.
    y_min_plot = max(base, y_min_plot)
    y_max_plot = max(y_max_plot, mhw + 0.05 * abs(mhw))
    y_range = y_max_plot - y_min_plot
    if not np.isfinite(y_range) or y_range <= 0:
        y_range = 1.0

    x_range = x_max_plot - x_min_plot
    fig_w, fig_h = _section_figsize(
        x_range,
        y_range,
        max_transect_length,
        max_depth_range,
        max_fig_width_cm,
        max_fig_height_cm,
    )

    t_dt = datenum_to_datetime64(slice_cube.t)
    highlight_idx = None
    highlight_tag = None
    if highlight_dates:
        # Only one highlight date per GIF; caller should pass a single date.
        highlight_dt = np.datetime64(highlight_dates[0])
        idx = int(np.searchsorted(t_dt, highlight_dt, side="left"))
        if idx >= len(t_dt):
            idx = len(t_dt) - 1
        if idx >= 0:
            highlight_idx = idx
            highlight_tag = str(t_dt[highlight_idx])[:10]
    frames = []
    for idx in range(0, nt, max(frame_stride, 1)):
        # Recompute stratigraphy using surveys up to this frame to preserve erosion timing.
        deposit_elev = compute_deposit_elev_1d(z_stack[: idx + 1, :])
        highlight_layer = highlight_idx if highlight_idx is not None and idx >= highlight_idx else None
        highlight_surface = None
        if highlight_layer is not None:
            intact = np.isfinite(z_stack[highlight_layer, :]) & np.isfinite(deposit_elev[highlight_layer, :])
            intact &= np.isclose(deposit_elev[highlight_layer, :], z_stack[highlight_layer, :], atol=1e-6)
            highlight_surface = np.where(intact, deposit_elev[highlight_layer, :], np.nan)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        base = min_elev - 0.01 * abs(min_elev)
        cumulative = base + np.zeros_like(x_m)
        layer_color = "red" if highlight_layer == 0 else colors[0]
        ax.fill_between(x_m, base, deposit_elev[0, :], color=layer_color, alpha=1.0)
        cumulative = deposit_elev[0, :]
        for tt in range(1, idx + 1):
            layer = np.maximum(0.0, deposit_elev[tt, :] - deposit_elev[tt - 1, :])
            upper = cumulative + layer
            layer_color = "red" if highlight_layer == tt else colors[tt]
            ax.fill_between(x_m, cumulative, upper, color=layer_color, alpha=1.0, zorder=1)
            cumulative = upper
        for tt in range(idx + 1):
            ax.plot(x_m, deposit_elev[tt, :], color="k", linewidth=0.6, alpha=0.8)
        if highlight_surface is not None:
            ax.plot(x_m, highlight_surface, color="red", linewidth=1.6, zorder=3)
        max_surface = np.nanmax(z_stack[: idx + 1, :], axis=0)
        if np.any(np.isfinite(max_surface)):
            ax.plot(x_m, max_surface, color="0.2", linestyle="--", linewidth=0.8, zorder=2)
        ax.plot(x_m, deposit_elev[idx, :], color="k", linewidth=1.6)
        if overlay_points is not None and overlay_points.size:
            ax.scatter(
                overlay_points[:, 0],
                overlay_points[:, 1],
                s=overlay_marker_size,
                c=overlay_color,
                edgecolors=overlay_edge_color,
                linewidths=0.35,
                zorder=4,
            )

        mhw_color = "#0b2e5b"
        dash_style = (0, (6, 3))
        ax.axhline(mhw, linestyle=dash_style, color=mhw_color, linewidth=0.5, zorder=0)
        ax.axhline(mlw, linestyle=dash_style, color=mhw_color, linewidth=0.5, zorder=0)

        date_str = str(t_dt[idx])[:10] if idx < len(t_dt) else ""
        ax.set_title(f"{section_title} | {date_str}")
        ax.set_xlabel("Distance [m]")
        ax.set_ylabel("Elevation [m]")
        if highlight_layer is not None and highlight_tag:
            ax.text(
                0.98,
                0.98,
                f"{highlight_tag}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                color="red",
                fontproperties=font,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.5},
            )
        apply_axes_font(ax, font)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(x_min_plot, x_max_plot)
        ax.set_ylim(y_min_plot, y_max_plot)

        fig.tight_layout()
        fig.canvas.draw()
        # Use buffer_rgba to support newer Matplotlib backends that dropped tostring_rgb.
        image = np.asarray(fig.canvas.buffer_rgba())
        if image.shape[-1] == 4:
            image = image[:, :, :3]
        frames.append(image)
        plt.close(fig)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    imageio.mimsave(out_path, frames, fps=fps, loop=0)


def run_transect_slice_plots(
    bathy_nc_path: str | Path,
    output_root: str | Path,
    transect_rows_km: np.ndarray,
    tick_spacing_m: float = 100.0,
    mhw: float = 0.358,
    mlw: float = -0.590,
    highlight_date: str | np.datetime64 | list[str | np.datetime64] | None = None,
    n_points: int = 400,
    initial_index: int = 0,
    dx: float = 20.0,
    target_crs: str = "EPSG:32618",
    max_fig_width_cm: float = 20.0,
    max_fig_height_cm: float = 5.0,
    make_gif: bool = False,
    gif_fps: int = 6,
    gif_stride: int = 1,
) -> dict[str, object]:
    """Plot 6-panel stratigraphy summaries for a list of transect slices.

    Args:
        bathy_nc_path: Path to the bathymetry cube netCDF file.
        output_root: Output root directory for plots.
        transect_rows_km: Array of transect endpoints in kilometers.
        tick_spacing_m: Tick spacing along transects (meters).
        mhw: Mean high water elevation.
        mlw: Mean low water elevation.
        highlight_date: Optional date string(s) or datetime64(s) to highlight deposits.
        n_points: Number of points used for each transect slice.
        initial_index: Initial stratigraphy index.
        dx: Grid spacing in meters.
        target_crs: CRS string for plot metadata.
        max_fig_width_cm: Maximum plot width for scaled sections (cm).
        max_fig_height_cm: Maximum plot height for scaled sections (cm).
        make_gif: Whether to create a stratigraphy GIF for each transect.
        gif_fps: Frames per second for GIF export.
        gif_stride: Step between frames (e.g., 2 uses every other survey).

    Returns:
        Dict containing the bathy cube, labels, and output directory for reuse.
    """
    labels = [_label_for_index(i) for i in range(len(transect_rows_km))]
    output_root = Path(output_root)
    bathy_cube = load_bathy_cube(bathy_nc_path)
    slice_cfg = StratigraphyConfig(initial_index=initial_index, dx=dx, target_crs=target_crs)
    slice_dir = output_root / "plots" / "python" / "stratigraphy" / "transect_slice"
    slice_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    z_min = np.inf
    z_max = -np.inf
    max_len_m = 0.0

    for (x1, y1, x2, y2), label in zip(transect_rows_km, labels):
        coords_m = np.array([[x1, y1], [x2, y2]], dtype=float) * 1000.0
        length_m = float(np.hypot(coords_m[1, 0] - coords_m[0, 0], coords_m[1, 1] - coords_m[0, 1]))
        max_len_m = max(max_len_m, length_m)
        transect = Transect(name=f"{label}", x=coords_m[:, 0], y=coords_m[:, 1])

        slice_cube = extract_transect_cube(bathy_cube, transect, n_points=n_points, name=transect.name)
        slice_result = compute_stratigraphy(slice_cube, slice_cfg)

        z_min = min(z_min, np.nanmin(slice_cube.z))
        z_max = max(z_max, np.nanmax(slice_cube.z))
        entries.append((label, transect, slice_cube, slice_result, length_m))

    plot_ylim = None
    max_depth_range = None
    if np.isfinite(z_min) and np.isfinite(z_max) and z_max > z_min:
        plot_ylim = (z_min - 0.01 * abs(z_min), z_max + 0.01 * abs(z_max))
        max_depth_range = plot_ylim[1] - plot_ylim[0]

    highlight_dates = _normalize_highlight_dates(highlight_date)

    for label, transect, slice_cube, slice_result, length_m in entries:
        section_title = f"{label}-{label}'"
        xticks_m = np.arange(0.0, length_m + 1e-6, tick_spacing_m)

        if not highlight_dates:
            plot_stratigraphy_stack(
                slice_cube,
                slice_result,
                slice_dir / f"strat_overview_section_{transect.name}.png",
                f"overview_{label}",
                time_vals=slice_cube.t,
                title_prefix=section_title,
            )
            plot_stacked_stratigraphy_section(
                slice_cube,
                slice_dir / f"strat_section_{label}.png",
                f"Cross-section {label}-{label}'",
                plot_xlim=None,
                plot_ylim=plot_ylim,
                x_ticks=xticks_m,
                mhw=mhw,
                mlw=mlw,
                highlight_date=None,
                max_transect_length=max_len_m,
                max_depth_range=max_depth_range,
                max_fig_width_cm=max_fig_width_cm,
                max_fig_height_cm=max_fig_height_cm,
                clip_x_to_data=True,
                clip_y_to_data=True,
            )
        else:
            for date_val in highlight_dates:
                highlight_str = date_val if isinstance(date_val, str) else str(np.datetime64(date_val))
                highlight_tag = highlight_str[:10]
                slice_highlight_dir = slice_dir / f"highlight_{highlight_tag}"
                slice_highlight_dir.mkdir(parents=True, exist_ok=True)
                plot_stratigraphy_stack(
                    slice_cube,
                    slice_result,
                    slice_highlight_dir / f"strat_{transect.name}_highlight_{highlight_tag}.png",
                    f"overview_slice_{label}",
                    time_vals=slice_cube.t,
                    highlight_date=date_val,
                    title_prefix=section_title,
                )
                plot_stacked_stratigraphy_section(
                    slice_cube,
                    slice_highlight_dir / f"strat_section_{label}_highlight_{highlight_tag}.png",
                    f"Cross-section {label}-{label}'",
                    plot_xlim=None,
                    plot_ylim=plot_ylim,
                    x_ticks=xticks_m,
                    mhw=mhw,
                    mlw=mlw,
                    highlight_date=date_val,
                    max_transect_length=max_len_m,
                    max_depth_range=max_depth_range,
                    max_fig_width_cm=max_fig_width_cm,
                    max_fig_height_cm=max_fig_height_cm,
                    clip_x_to_data=True,
                    clip_y_to_data=True,
                )

        if make_gif:
            if not highlight_dates:
                gif_path = slice_dir / f"strat_section_{label}.gif"
                _write_stratigraphy_gif(
                    slice_cube,
                    gif_path,
                    section_title,
                    mhw=mhw,
                    mlw=mlw,
                    plot_ylim=plot_ylim,
                    highlight_dates=None,
                    max_transect_length=max_len_m,
                    max_depth_range=max_depth_range,
                    max_fig_width_cm=max_fig_width_cm,
                    max_fig_height_cm=max_fig_height_cm,
                    fps=gif_fps,
                    frame_stride=gif_stride,
                )
            else:
                for date_val in highlight_dates:
                    highlight_str = date_val if isinstance(date_val, str) else str(np.datetime64(date_val))
                    highlight_tag = highlight_str[:10]
                    slice_highlight_dir = slice_dir / f"highlight_{highlight_tag}"
                    slice_highlight_dir.mkdir(parents=True, exist_ok=True)
                    gif_path = slice_highlight_dir / f"strat_section_{label}_highlight_{highlight_tag}.gif"
                    _write_stratigraphy_gif(
                        slice_cube,
                        gif_path,
                        section_title,
                        mhw=mhw,
                        mlw=mlw,
                        plot_ylim=plot_ylim,
                        highlight_dates=[date_val],
                        max_transect_length=max_len_m,
                        max_depth_range=max_depth_range,
                        max_fig_width_cm=max_fig_width_cm,
                        max_fig_height_cm=max_fig_height_cm,
                        fps=gif_fps,
                        frame_stride=gif_stride,
                    )

    return {
        "bathy_cube": bathy_cube,
        "labels": labels,
        "transect_rows_km": transect_rows_km,
        "output_root": output_root,
        "slice_dir": slice_dir,
        "highlight_dates": highlight_dates,
    }


def run_transect_plots(
    source: str,
    bathy_nc_path: str | Path,
    output_root: str | Path,
    transect_rows_km: np.ndarray | None = None,
    shp_dir: str | Path | None = None,
    tick_spacing_m: float = 100.0,
    mhw: float = 0.358,
    mlw: float = -0.590,
    highlight_date: str | np.datetime64 | list[str | np.datetime64] | None = None,
    n_points: int = 400,
    initial_index: int = 0,
    dx: float = 20.0,
    target_crs: str = "EPSG:32618",
    max_fig_width_cm: float = 20.0,
    max_fig_height_cm: float = 5.0,
    make_gif: bool = False,
    gif_fps: int = 6,
    gif_stride: int = 1,
    source_crs: str | None = None,
    prompt_if_missing_crs: bool = True,
    prefer_xy_columns: bool = False,
    smooth_transects: bool = True,
    spline_points: int = 1000,
    distance_mode: str = "curvy",
    map_tick_length_km: float = 0.02,
    map_output_name: str = "transect_location_plan_shapefiles.png",
    overlay_csv: str | Path | None = None,
    overlay_x_column: str | None = None,
    overlay_y_column: str | None = None,
    overlay_color: str = "#00E5FF",
    overlay_edge_color: str = "black",
    overlay_marker_size: float = 8.0,
    overlay_align_start: bool = False,
) -> dict[str, object]:
    """Run transect plots from endpoints or shapefiles with a unified interface.

    Args:
        source: "endpoints" for transect rows or "shapefiles" for shapefile inputs.
        bathy_nc_path: Path to the bathymetry cube netCDF file.
        output_root: Output root directory for plots.
        transect_rows_km: Array of transect endpoints (required for endpoints).
        shp_dir: Directory containing shapefiles (required for shapefiles).
        tick_spacing_m: Tick spacing along transects (meters).
        mhw: Mean high water elevation.
        mlw: Mean low water elevation.
        highlight_date: Optional highlight date(s).
        n_points: Number of points used for each transect slice.
        initial_index: Initial stratigraphy index.
        dx: Grid spacing in meters.
        target_crs: Target CRS for transect endpoints.
        max_fig_width_cm: Maximum plot width for scaled sections (cm).
        max_fig_height_cm: Maximum plot height for scaled sections (cm).
        make_gif: Whether to create stratigraphy GIFs.
        gif_fps: Frames per second for GIF export.
        gif_stride: Step between frames (e.g., 2 uses every other survey).
        source_crs: Optional CRS for shapefile inputs.
        prompt_if_missing_crs: Prompt if shapefile CRS metadata is missing.
        prefer_xy_columns: Prefer DBF x/y columns over line geometry.
        smooth_transects: Whether to spline-smooth shapefile transects.
        spline_points: Number of spline points for smoothed transects.
        distance_mode: "curvy" uses polylines; "straight" uses endpoints.
        map_tick_length_km: Tick length for transect ticks in maps.
        map_output_name: Filename for shapefile transect maps.
        overlay_csv: Optional profile CSV or directory of label-prefixed CSV files.
        overlay_x_column: Optional CSV column containing profile distance in meters.
        overlay_y_column: Optional CSV column containing elevation in meters.
        overlay_color: Color used for profile overlay points.
        overlay_edge_color: Edge color used for profile overlay points.
        overlay_marker_size: Overlay marker area in points squared.
        overlay_align_start: Shift each overlay so its minimum distance is zero.

    Returns:
        Dict containing outputs from the selected workflow.
    """
    source = source.lower().strip()
    if source == "endpoints":
        if transect_rows_km is None:
            raise ValueError("transect_rows_km is required when source='endpoints'.")
        return run_transect_slice_plots(
            bathy_nc_path=bathy_nc_path,
            output_root=output_root,
            transect_rows_km=transect_rows_km,
            tick_spacing_m=tick_spacing_m,
            mhw=mhw,
            mlw=mlw,
            highlight_date=highlight_date,
            n_points=n_points,
            initial_index=initial_index,
            dx=dx,
            target_crs=target_crs,
            max_fig_width_cm=max_fig_width_cm,
            max_fig_height_cm=max_fig_height_cm,
            make_gif=make_gif,
            gif_fps=gif_fps,
            gif_stride=gif_stride,
        )
    if source == "shapefiles":
        if shp_dir is None:
            raise ValueError("shp_dir is required when source='shapefiles'.")
        return run_shapefile_transect_plots(
            bathy_nc_path=bathy_nc_path,
            output_root=output_root,
            shp_dir=shp_dir,
            target_crs=target_crs,
            source_crs=source_crs,
            prompt_if_missing_crs=prompt_if_missing_crs,
            prefer_xy_columns=prefer_xy_columns,
            smooth_transects=smooth_transects,
            spline_points=spline_points,
            distance_mode=distance_mode,
            n_points=n_points,
            tick_spacing_m=tick_spacing_m,
            map_tick_length_km=map_tick_length_km,
            mhw=mhw,
            mlw=mlw,
            highlight_date=highlight_date,
            initial_index=initial_index,
            dx=dx,
            max_fig_width_cm=max_fig_width_cm,
            max_fig_height_cm=max_fig_height_cm,
            map_output_name=map_output_name,
            make_gif=make_gif,
            gif_fps=gif_fps,
            gif_stride=gif_stride,
            overlay_csv=overlay_csv,
            overlay_x_column=overlay_x_column,
            overlay_y_column=overlay_y_column,
            overlay_color=overlay_color,
            overlay_edge_color=overlay_edge_color,
            overlay_marker_size=overlay_marker_size,
            overlay_align_start=overlay_align_start,
        )
    raise ValueError("source must be 'endpoints' or 'shapefiles'.")


def plot_transect_location_plan(
    bathy_cube,
    transect_rows_km: np.ndarray,
    output_root: str | Path,
    labels: list[str] | None = None,
    mlw: float = 0.0,
    tick_spacing_m: float = 100.0,
    tick_length_km: float = 0.02,
    cmap_name: str = "kg2",
    make_inset: bool = False,
    inset_scale: float = 0.5,
    inset_text_scale: float = 1.4,
) -> Path:
    """Plot transect locations on the latest bathymetry surface.

    Args:
        bathy_cube: BathyCube loaded from the bathymetry stack.
        transect_rows_km: Array of transect endpoints in kilometers.
        output_root: Output root directory for plots.
        labels: Optional list of transect labels.
        mlw: Mean low water contour in meters.
        tick_spacing_m: Tick spacing along transects (meters).
        tick_length_km: Tick length for transect tick marks (km).
        cmap_name: Colormap name for bathymetry rendering.
        make_inset: Whether to save a smaller inset version of the plan.
        inset_scale: Linear scale factor for inset figure size.
        inset_text_scale: Scale factor for inset text sizes.

    Returns:
        Path to the saved plan figure.
    """
    apply_global_style()
    font = get_plot_font()

    def _plot_transect_ticks(ax, x1, y1, x2, y2, spacing_km, tick_len_km) -> None:
        """Draw small perpendicular ticks along a transect line."""
        dx = x2 - x1
        dy = y2 - y1
        length = float(np.hypot(dx, dy))
        if length <= 0.0:
            return
        ux, uy = dx / length, dy / length
        px, py = -uy, ux
        for s in np.arange(0.0, length + 1e-9, spacing_km):
            cx = x1 + ux * s
            cy = y1 + uy * s
            x0 = cx - px * tick_len_km * 0.5
            x1t = cx + px * tick_len_km * 0.5
            y0 = cy - py * tick_len_km * 0.5
            y1t = cy + py * tick_len_km * 0.5
            ax.plot([x0, x1t], [y0, y1t], color="k", linewidth=0.6, zorder=2)

    labels = labels or [_label_for_index(i) for i in range(len(transect_rows_km))]
    output_root = Path(output_root)

    x_km = bathy_cube.x / 1000.0
    y_km = bathy_cube.y / 1000.0
    z_last = bathy_cube.z[:, :, -1]

    fig, ax = plt.subplots(figsize=(9, 7))
    levels = np.arange(-10.0, 5.01, 0.2)
    cmap = get_named_colormap(cmap_name) if not cmap_name.endswith(".clrmap") else load_clrmap_file(cmap_name)
    cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, extend="both")
    ax.contour(x_km, y_km, z_last, levels=[mlw], colors=["0.5"], linewidths=1.0)
    ax.contour(x_km, y_km, z_last, levels=[-6.0], colors="k", linestyles=":", linewidths=0.5)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Easting [km]")
    ax.set_ylabel("Northing [km]")
    ax.set_title("Transect Locations (Latest Bathymetry)")
    apply_axes_font(ax, font)
    ax.grid(True, color="0.5", alpha=0.4)
    ax.set_axisbelow(False)

    cb = fig.colorbar(cf, ax=ax)
    cb.set_label("Depth [m]", fontproperties=font)
    for tick in cb.ax.get_yticklabels():
        tick.set_fontproperties(font)

    tick_spacing_km = tick_spacing_m / 1000.0
    for (x1, y1, x2, y2), label in zip(transect_rows_km, labels):
        ax.plot([x1, x2], [y1, y2], "-k", linewidth=1.0)
        _plot_transect_ticks(ax, x1, y1, x2, y2, tick_spacing_km, tick_length_km)
        ax.scatter([x1, x2], [y1, y2], s=20, c="w", edgecolors="k", zorder=3)
        ax.text(x1 - 0.03, y1 + 0.03, label, fontproperties=font, color="k")
        ax.text(x2 + 0.03, y2 - 0.03, f"{label}'", fontproperties=font, color="k")

    plan_dir = output_root / "plots" / "python" / "stratigraphy" / "transect_slice"
    plan_dir.mkdir(parents=True, exist_ok=True)
    plan_path = plan_dir / "transect_location_plan.png"
    fig.savefig(plan_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)
    print(f"Saved transect plan to: {plan_path}")

    if make_inset:
        inset_figsize = (9 * inset_scale, 7 * inset_scale)
        fig, ax = plt.subplots(figsize=inset_figsize)
        levels = np.arange(-10.0, 5.01, 0.2)
        cmap = get_named_colormap(cmap_name) if not cmap_name.endswith(".clrmap") else load_clrmap_file(cmap_name)
        cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, extend="both")
        ax.contour(x_km, y_km, z_last, levels=[mlw], colors=["0.5"], linewidths=1.0)
        ax.contour(x_km, y_km, z_last, levels=[-6.0], colors="k", linestyles=":", linewidths=0.5)

        ax.set_aspect("equal", adjustable="box")
        base_size = float(plt.rcParams.get("font.size", 10.0))
        inset_size = base_size * inset_text_scale
        ax.set_xlabel("Easting [km]", fontsize=inset_size, fontproperties=font)
        ax.set_ylabel("Northing [km]", fontsize=inset_size, fontproperties=font)
        ax.set_title("Transect Locations (Latest Bathymetry)", fontsize=inset_size, fontproperties=font)
        ax.grid(True, color="0.5", alpha=0.4)
        ax.set_axisbelow(False)
        ax.tick_params(labelsize=inset_size * 0.9)

        cb = fig.colorbar(cf, ax=ax)
        cb.set_label("Depth [m]", fontsize=inset_size, fontproperties=font)
        for tick in cb.ax.get_yticklabels():
            tick.set_fontproperties(font)
            tick.set_fontsize(inset_size * 0.9)

        tick_spacing_km = tick_spacing_m / 1000.0
        for (x1, y1, x2, y2), label in zip(transect_rows_km, labels):
            ax.plot([x1, x2], [y1, y2], "-k", linewidth=1.0)
            _plot_transect_ticks(ax, x1, y1, x2, y2, tick_spacing_km, tick_length_km)
            ax.scatter([x1, x2], [y1, y2], s=20, c="w", edgecolors="k", zorder=3)
            ax.text(x1 - 0.03, y1 + 0.03, label, fontproperties=font, fontsize=inset_size, color="k")
            ax.text(x2 + 0.03, y2 - 0.03, f"{label}'", fontproperties=font, fontsize=inset_size, color="k")

        inset_path = plan_dir / "transect_location_plan_inset.png"
        fig.savefig(inset_path, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved inset transect plan to: {inset_path}")
    return plan_path


def _plot_ticks_along_polyline(
    ax,
    x_km: np.ndarray,
    y_km: np.ndarray,
    spacing_km: float,
    tick_len_km: float,
) -> None:
    """Draw perpendicular tick marks along a polyline at fixed spacing."""
    if spacing_km <= 0 or tick_len_km <= 0 or x_km.size < 2:
        return

    pts = np.column_stack([x_km, y_km])
    seg_len = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    s = np.concatenate([[0.0], np.cumsum(seg_len)])
    if s[-1] <= 0:
        return

    for dist in np.arange(0.0, s[-1] + 1e-9, spacing_km):
        idx = int(np.searchsorted(s, dist, side="right") - 1)
        idx = min(max(idx, 0), len(seg_len) - 1)
        if seg_len[idx] <= 0:
            continue
        frac = (dist - s[idx]) / seg_len[idx]
        x0, y0 = pts[idx]
        x1, y1 = pts[idx + 1]
        cx = x0 + frac * (x1 - x0)
        cy = y0 + frac * (y1 - y0)
        ux = (x1 - x0) / seg_len[idx]
        uy = (y1 - y0) / seg_len[idx]
        px, py = -uy, ux
        tx0 = cx - px * tick_len_km * 0.5
        ty0 = cy - py * tick_len_km * 0.5
        tx1 = cx + px * tick_len_km * 0.5
        ty1 = cy + py * tick_len_km * 0.5
        ax.plot([tx0, tx1], [ty0, ty1], color="k", linewidth=0.6, zorder=2)


def plot_transect_location_plan_from_transects(
    bathy_cube,
    transects: list[Transect],
    output_root: str | Path,
    labels: list[str] | None = None,
    mlw: float = 0.0,
    tick_spacing_m: float = 100.0,
    tick_length_km: float = 0.02,
    cmap_name: str = "kg2",
    map_output_name: str = "transect_location_plan_shapefiles.png",
    make_inset: bool = False,
    inset_scale: float = 0.5,
    inset_text_scale: float = 1.4,
) -> Path:
    """Plot transect polylines on the latest bathymetry surface.

    Args:
        bathy_cube: BathyCube loaded from the bathymetry stack.
        transects: List of transects in projected coordinates (meters).
        output_root: Output root directory for plots.
        labels: Optional list of transect labels (A, B, C...).
        mlw: Mean low water contour in meters.
        tick_spacing_m: Tick spacing along transects (meters).
        tick_length_km: Tick length for transect tick marks (km).
        cmap_name: Colormap name for bathymetry rendering.
        map_output_name: Filename for the saved map.
        make_inset: Whether to save a smaller inset version of the plan.
        inset_scale: Linear scale factor for inset figure size.
        inset_text_scale: Scale factor for inset text sizes.

    Returns:
        Path to the saved plan figure.
    """
    apply_global_style()
    font = get_plot_font()
    labels = labels or [_label_for_index(i) for i in range(len(transects))]
    output_root = Path(output_root)

    x_km = bathy_cube.x / 1000.0
    y_km = bathy_cube.y / 1000.0
    z_last = bathy_cube.z[:, :, -1]

    fig, ax = plt.subplots(figsize=(9, 7))
    levels = np.arange(-10.0, 5.01, 0.2)
    cmap = get_named_colormap(cmap_name) if not cmap_name.endswith(".clrmap") else load_clrmap_file(cmap_name)
    cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, extend="both")
    ax.contour(x_km, y_km, z_last, levels=[mlw], colors=["0.5"], linewidths=1.0)
    ax.contour(x_km, y_km, z_last, levels=[-6.0], colors="k", linestyles=":", linewidths=0.5)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Easting [km]")
    ax.set_ylabel("Northing [km]")
    ax.set_title("Transect Locations (Latest Bathymetry)")
    apply_axes_font(ax, font)
    ax.grid(True, color="0.5", alpha=0.4)
    ax.set_axisbelow(False)

    cb = fig.colorbar(cf, ax=ax)
    cb.set_label("Depth [m]", fontproperties=font)
    for tick in cb.ax.get_yticklabels():
        tick.set_fontproperties(font)

    tick_spacing_km = tick_spacing_m / 1000.0
    for transect, label in zip(transects, labels):
        x_line = np.asarray(transect.x, dtype=float) / 1000.0
        y_line = np.asarray(transect.y, dtype=float) / 1000.0
        ax.plot(x_line, y_line, "-k", linewidth=1.0)
        _plot_ticks_along_polyline(ax, x_line, y_line, tick_spacing_km, tick_length_km)
        ax.scatter([x_line[0], x_line[-1]], [y_line[0], y_line[-1]], s=20, c="w", edgecolors="k", zorder=3)
        ax.text(x_line[0] - 0.03, y_line[0] + 0.03, label, fontproperties=font, color="k")
        ax.text(x_line[-1] + 0.03, y_line[-1] - 0.03, f"{label}'", fontproperties=font, color="k")

    map_dir = output_root / "plots" / "python" / "stratigraphy" / "transect_slice"
    map_dir.mkdir(parents=True, exist_ok=True)
    map_path = map_dir / map_output_name
    fig.savefig(map_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)
    print(f"Saved transect plan to: {map_path}")

    if make_inset:
        inset_figsize = (9 * inset_scale, 7 * inset_scale)
        fig, ax = plt.subplots(figsize=inset_figsize)
        levels = np.arange(-10.0, 5.01, 0.2)
        cmap = get_named_colormap(cmap_name) if not cmap_name.endswith(".clrmap") else load_clrmap_file(cmap_name)
        cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, extend="both")
        ax.contour(x_km, y_km, z_last, levels=[mlw], colors=["0.5"], linewidths=1.0)
        ax.contour(x_km, y_km, z_last, levels=[-6.0], colors="k", linestyles=":", linewidths=0.5)

        ax.set_aspect("equal", adjustable="box")
        base_size = float(plt.rcParams.get("font.size", 10.0))
        inset_size = base_size * inset_text_scale
        ax.set_xlabel("Easting [km]", fontsize=inset_size, fontproperties=font)
        ax.set_ylabel("Northing [km]", fontsize=inset_size, fontproperties=font)
        ax.set_title("Transect Locations (Latest Bathymetry)", fontsize=inset_size, fontproperties=font)
        ax.grid(True, color="0.5", alpha=0.4)
        ax.set_axisbelow(False)
        ax.tick_params(labelsize=inset_size * 0.9)

        cb = fig.colorbar(cf, ax=ax)
        cb.set_label("Depth [m]", fontsize=inset_size, fontproperties=font)
        for tick in cb.ax.get_yticklabels():
            tick.set_fontproperties(font)
            tick.set_fontsize(inset_size * 0.9)

        tick_spacing_km = tick_spacing_m / 1000.0
        for transect, label in zip(transects, labels):
            x_line = np.asarray(transect.x, dtype=float) / 1000.0
            y_line = np.asarray(transect.y, dtype=float) / 1000.0
            ax.plot(x_line, y_line, "-k", linewidth=1.0)
            _plot_ticks_along_polyline(ax, x_line, y_line, tick_spacing_km, tick_length_km)
            ax.scatter([x_line[0], x_line[-1]], [y_line[0], y_line[-1]], s=20, c="w", edgecolors="k", zorder=3)
            ax.text(x_line[0] - 0.03, y_line[0] + 0.03, label, fontproperties=font, fontsize=inset_size, color="k")
            ax.text(x_line[-1] + 0.03, y_line[-1] - 0.03, f"{label}'", fontproperties=font, fontsize=inset_size, color="k")

        inset_name = map_output_name.replace(".png", "_inset.png")
        inset_path = map_dir / inset_name
        fig.savefig(inset_path, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"Saved inset transect plan to: {inset_path}")
    return map_path


def run_shapefile_transect_plots(
    bathy_nc_path: str | Path,
    output_root: str | Path,
    shp_dir: str | Path,
    target_crs: str = "EPSG:32618",
    source_crs: str | None = None,
    prompt_if_missing_crs: bool = True,
    prefer_xy_columns: bool = False,
    smooth_transects: bool = True,
    spline_points: int = 1000,
    distance_mode: str = "curvy",
    n_points: int = 400,
    tick_spacing_m: float = 100.0,
    map_tick_length_km: float = 0.02,
    mhw: float = 0.358,
    mlw: float = -0.590,
    highlight_date: str | np.datetime64 | list[str | np.datetime64] | None = None,
    initial_index: int = 0,
    dx: float = 20.0,
    max_fig_width_cm: float = 20.0,
    max_fig_height_cm: float = 5.0,
    map_output_name: str = "transect_location_plan_shapefiles.png",
    make_gif: bool = False,
    gif_fps: int = 6,
    gif_stride: int = 1,
    overlay_csv: str | Path | None = None,
    overlay_x_column: str | None = None,
    overlay_y_column: str | None = None,
    overlay_color: str = "#00E5FF",
    overlay_edge_color: str = "black",
    overlay_marker_size: float = 8.0,
    overlay_align_start: bool = False,
) -> dict[str, object]:
    """Plot stratigraphy sections and a location map for shapefile transects.

    Args:
        bathy_nc_path: Path to the bathymetry cube netCDF file.
        output_root: Output root directory for plots.
        shp_dir: Directory containing transect shapefiles.
        target_crs: Target CRS for reprojection.
        source_crs: Optional source CRS override.
        prompt_if_missing_crs: Prompt if shapefile CRS metadata is missing.
        prefer_xy_columns: Prefer DBF x/y columns over line geometry.
        smooth_transects: Whether to spline-smooth shapefile transects.
        spline_points: Number of spline points used for smoothing.
        distance_mode: "curvy" for along-transect distances, "straight" for endpoints.
        n_points: Number of samples along each transect.
        tick_spacing_m: Tick spacing along transects (meters).
        map_tick_length_km: Tick length for map ticks (km).
        mhw: Mean high water elevation.
        mlw: Mean low water elevation.
        highlight_date: Optional date string(s) or datetime64(s) to highlight deposits.
        initial_index: Initial stratigraphy index.
        dx: Grid spacing in meters.
        max_fig_width_cm: Maximum plot width for scaled sections (cm).
        max_fig_height_cm: Maximum plot height for scaled sections (cm).
        map_output_name: Filename for the transect location plan.
        make_gif: Whether to create a stratigraphy GIF for each transect.
        gif_fps: Frames per second for GIF export.
        gif_stride: Step between frames (e.g., 2 uses every other survey).
        overlay_csv: Optional profile CSV or directory of label-prefixed CSV files.
        overlay_x_column: Optional CSV column containing profile distance in meters.
        overlay_y_column: Optional CSV column containing elevation in meters.
        overlay_color: Color used for profile overlay points.
        overlay_edge_color: Edge color used for profile overlay points.
        overlay_marker_size: Overlay marker area in points squared.
        overlay_align_start: Shift each overlay so its minimum distance is zero.

    Returns:
        Dict containing bathy cube, transects, and output paths.
    """
    output_root = Path(output_root)
    bathy_cube = load_bathy_cube(bathy_nc_path)
    transects = load_transects_from_shapefiles(
        shp_dir=shp_dir,
        target_crs=target_crs,
        source_crs=source_crs,
        prompt_if_missing_crs=prompt_if_missing_crs,
        prefer_xy_columns=prefer_xy_columns,
    )

    slice_cfg = StratigraphyConfig(initial_index=initial_index, dx=dx, target_crs=target_crs)
    slice_dir = output_root / "plots" / "python" / "stratigraphy" / "transect_slice" / "shapefiles"
    slice_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    valid_transects = []
    skipped_transects = []
    label_map = {}
    z_min = np.inf
    z_max = -np.inf
    max_len_m = 0.0

    label_index = 0
    for transect in transects:
        x_raw, y_raw = _unique_xy(transect.x, transect.y)
        if x_raw.size < 2:
            print(f"WARNING: input file {transect.name}.shp has fewer than 2 points and was ignored.")
            skipped_transects.append(transect.name)
            continue

        if smooth_transects:
            try:
                x_raw, y_raw = _smooth_transect_xy(x_raw, y_raw, spline_points)
            except ValueError as exc:
                print(f"WARNING: input file {transect.name}.shp could not be smoothed ({exc}).")

        if distance_mode.lower() == "straight":
            x_raw = np.array([x_raw[0], x_raw[-1]])
            y_raw = np.array([y_raw[0], y_raw[-1]])

        transect_use = Transect(name=transect.name, x=x_raw, y=y_raw)

        if not _transect_hits_domain(bathy_cube, transect_use):
            print(
                f"WARNING: input file {transect.name}.shp does not intersect the surveyed area and "
                "was therefore ignored."
            )
            skipped_transects.append(transect.name)
            continue

        coords = np.column_stack([transect_use.x, transect_use.y])
        seg_len = np.hypot(np.diff(coords[:, 0]), np.diff(coords[:, 1]))
        length_m = float(np.sum(seg_len))
        max_len_m = max(max_len_m, length_m)

        label = _label_for_index(label_index)
        label_index += 1
        label_map[label] = transect.name
        transect_label = Transect(name=label, x=transect_use.x, y=transect_use.y)

        slice_cube = extract_transect_cube(
            bathy_cube, transect_label, n_points=n_points, name=transect_label.name
        )
        if not np.isfinite(slice_cube.z).any():
            print(
                f"WARNING: input file {transect.name}.shp does not contain valid data within the surveyed "
                "area and was therefore ignored."
            )
            skipped_transects.append(transect.name)
            continue

        slice_result = compute_stratigraphy(slice_cube, slice_cfg)

        z_min = min(z_min, np.nanmin(slice_cube.z))
        z_max = max(z_max, np.nanmax(slice_cube.z))
        entries.append((label, transect_label, slice_cube, slice_result, length_m))
        valid_transects.append(transect_label)

    plot_ylim = None
    max_depth_range = None
    if np.isfinite(z_min) and np.isfinite(z_max) and z_max > z_min:
        plot_ylim = (z_min - 0.01 * abs(z_min), z_max + 0.01 * abs(z_max))
        max_depth_range = plot_ylim[1] - plot_ylim[0]

    highlight_dates = _normalize_highlight_dates(highlight_date)
    overlay_enabled = overlay_csv is not None
    transect_output_dir = slice_dir / "overlay" if overlay_enabled else slice_dir
    transect_output_dir.mkdir(parents=True, exist_ok=True)

    for label, transect, slice_cube, slice_result, length_m in entries:
        section_title = f"{label}-{label}'"
        xticks_m = np.arange(0.0, length_m + 1e-6, tick_spacing_m)
        overlay_points = None
        if overlay_csv is not None:
            overlay_points = _load_profile_overlay(
                overlay_csv,
                label,
                x_column=overlay_x_column,
                y_column=overlay_y_column,
                align_start=overlay_align_start,
            )

        if not highlight_dates:
            plot_stratigraphy_stack(
                slice_cube,
                slice_result,
                _transect_output_path(
                    transect_output_dir, f"strat_overview_section_{label}.png", overlay_enabled
                ),
                f"overview_{label}",
                time_vals=slice_cube.t,
                title_prefix=section_title,
            )
            plot_stacked_stratigraphy_section(
                slice_cube,
                _transect_output_path(transect_output_dir, f"strat_section_{label}.png", overlay_enabled),
                f"Cross-section {section_title}",
                plot_xlim=None,
                plot_ylim=plot_ylim,
                x_ticks=xticks_m,
                mhw=mhw,
                mlw=mlw,
                highlight_date=None,
                max_transect_length=max_len_m,
                max_depth_range=max_depth_range,
                max_fig_width_cm=max_fig_width_cm,
                max_fig_height_cm=max_fig_height_cm,
                clip_x_to_data=True,
                clip_y_to_data=True,
                overlay_points=overlay_points,
                overlay_color=overlay_color,
                overlay_edge_color=overlay_edge_color,
                overlay_marker_size=overlay_marker_size,
            )
        else:
            for date_val in highlight_dates:
                highlight_str = date_val if isinstance(date_val, str) else str(np.datetime64(date_val))
                highlight_tag = highlight_str[:10]
                slice_highlight_dir = transect_output_dir / f"highlight_{highlight_tag}"
                slice_highlight_dir.mkdir(parents=True, exist_ok=True)
                plot_stratigraphy_stack(
                    slice_cube,
                    slice_result,
                    _transect_output_path(
                        slice_highlight_dir,
                        f"strat_{label}_highlight_{highlight_tag}.png",
                        overlay_enabled,
                    ),
                    f"overview_slice_{label}",
                    time_vals=slice_cube.t,
                    highlight_date=date_val,
                    title_prefix=section_title,
                )
                plot_stacked_stratigraphy_section(
                    slice_cube,
                    _transect_output_path(
                        slice_highlight_dir,
                        f"strat_section_{label}_highlight_{highlight_tag}.png",
                        overlay_enabled,
                    ),
                    f"Cross-section {section_title}",
                    plot_xlim=None,
                    plot_ylim=plot_ylim,
                    x_ticks=xticks_m,
                    mhw=mhw,
                    mlw=mlw,
                    highlight_date=date_val,
                    max_transect_length=max_len_m,
                    max_depth_range=max_depth_range,
                    max_fig_width_cm=max_fig_width_cm,
                    max_fig_height_cm=max_fig_height_cm,
                    clip_x_to_data=True,
                    clip_y_to_data=True,
                    overlay_points=overlay_points,
                    overlay_color=overlay_color,
                    overlay_edge_color=overlay_edge_color,
                    overlay_marker_size=overlay_marker_size,
                )

        if make_gif:
            if not highlight_dates:
                gif_path = _transect_output_path(
                    transect_output_dir, f"strat_section_{label}.gif", overlay_enabled
                )
                _write_stratigraphy_gif(
                    slice_cube,
                    gif_path,
                    section_title,
                    mhw=mhw,
                    mlw=mlw,
                    plot_ylim=plot_ylim,
                    highlight_dates=None,
                    max_transect_length=max_len_m,
                    max_depth_range=max_depth_range,
                    max_fig_width_cm=max_fig_width_cm,
                    max_fig_height_cm=max_fig_height_cm,
                    fps=gif_fps,
                    frame_stride=gif_stride,
                    overlay_points=overlay_points,
                    overlay_color=overlay_color,
                    overlay_edge_color=overlay_edge_color,
                    overlay_marker_size=overlay_marker_size,
                )
            else:
                for date_val in highlight_dates:
                    highlight_str = date_val if isinstance(date_val, str) else str(np.datetime64(date_val))
                    highlight_tag = highlight_str[:10]
                    slice_highlight_dir = transect_output_dir / f"highlight_{highlight_tag}"
                    slice_highlight_dir.mkdir(parents=True, exist_ok=True)
                    gif_path = _transect_output_path(
                        slice_highlight_dir,
                        f"strat_section_{label}_highlight_{highlight_tag}.gif",
                        overlay_enabled,
                    )
                    _write_stratigraphy_gif(
                        slice_cube,
                        gif_path,
                        section_title,
                        mhw=mhw,
                        mlw=mlw,
                        plot_ylim=plot_ylim,
                        highlight_dates=[date_val],
                        max_transect_length=max_len_m,
                        max_depth_range=max_depth_range,
                        max_fig_width_cm=max_fig_width_cm,
                        max_fig_height_cm=max_fig_height_cm,
                        fps=gif_fps,
                        frame_stride=gif_stride,
                        overlay_points=overlay_points,
                        overlay_color=overlay_color,
                        overlay_edge_color=overlay_edge_color,
                        overlay_marker_size=overlay_marker_size,
                    )

    if valid_transects:
        map_path = plot_transect_location_plan_from_transects(
            bathy_cube,
            valid_transects,
            output_root=output_root,
            labels=[t.name for t in valid_transects],
            mlw=mlw,
            tick_spacing_m=tick_spacing_m,
            tick_length_km=map_tick_length_km,
            cmap_name="kg2",
            map_output_name=map_output_name,
        )
    else:
        map_path = output_root / "plots" / "python" / "stratigraphy" / "transect_slice" / map_output_name
        print("WARNING: no valid transects found for plotting.")

    return {
        "bathy_cube": bathy_cube,
        "transects": valid_transects,
        "skipped_transects": skipped_transects,
        "output_root": output_root,
        "slice_dir": slice_dir,
        "map_path": map_path,
        "label_map": label_map,
        "highlight_dates": highlight_dates,
    }


def _ensure_local_imports() -> Path:
    """Ensure the local code directory is on ``sys.path``.

    Returns:
        Path to the resolved code directory.
    """
    cwd = Path.cwd()
    code_dir = cwd if (cwd / "bathy_formatter.py").exists() else (cwd / "code")
    if code_dir.exists() and str(code_dir) not in sys.path:
        sys.path.insert(0, str(code_dir))
    return code_dir


def _ensure_packages(required_modules: dict[str, str]) -> None:
    """Install any missing packages and import required modules.

    Args:
        required_modules: Mapping of import name to pip package name.
    """
    missing_packages = []
    for module_name, package_name in required_modules.items():
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError:
            missing_packages.append(package_name)

    if missing_packages:
        missing_packages = sorted(set(missing_packages))
        print(f"Installing missing packages into this kernel: {missing_packages}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", *missing_packages])

    for module_name in required_modules:
        importlib.import_module(module_name)


def startup_check(required_modules: dict[str, str] | None = None) -> dict[str, object]:
    """Ensure local imports, dependencies, and kernel selection."""
    required_modules = required_modules or DEFAULT_REQUIRED_MODULES
    code_dir = _ensure_local_imports()
    _ensure_packages(required_modules)

    exe_norm = sys.executable.replace("\\", "/").lower()
    if "/.venv/" not in exe_norm:
        print("WARNING: Notebook is not using the project .venv kernel.")
        print("Select kernel: Python (.venv bathy2strat)")

    print(f"Python executable: {sys.executable}")
    print("Startup dependency check complete.")
    return {"code_dir": code_dir, "required_modules": required_modules}


def _datenum_to_datetime(datenum_value: float) -> datetime:
    """Convert a MATLAB datenum to a Python ``datetime``.

    Args:
        datenum_value: MATLAB datenum float.

    Returns:
        Converted ``datetime``.
    """
    value = float(np.asarray(datenum_value).reshape(-1)[0])
    return datetime.fromordinal(int(value)) + timedelta(days=value % 1) - timedelta(days=366)


def run_morphodynamics_vs_wave_power(
    wave_power_stats_file: str | Path,
    morpho_csv_path: str | Path,
    stacked_plot_dir: str | Path,
) -> pd.DataFrame:
    """Compare morphodynamic change intervals with wave power and write plots.

    Args:
        wave_power_stats_file: Tab-delimited wave power summary file path.
        morpho_csv_path: CSV path with morphodynamic step statistics.
        stacked_plot_dir: Output directory for stacked plots.

    Returns:
        Merged DataFrame of wave power and morphodynamic step stats.
    """
    stacked_plot_dir = Path(stacked_plot_dir)
    stacked_plot_dir.mkdir(parents=True, exist_ok=True)

    morpho_df = pd.read_csv(morpho_csv_path)
    step_df = morpho_df[morpho_df["scope"] == "step"].copy()
    step_df["end_date"] = step_df["time_datenum"].apply(_datenum_to_datetime)
    step_df["end_date"] = pd.to_datetime(step_df["end_date"]).dt.normalize()

    wave_df = pd.read_csv(wave_power_stats_file, sep="\t")
    wave_df["end_date"] = pd.to_datetime(wave_df["end_date"]).dt.normalize()
    wave_df["start_date"] = pd.to_datetime(wave_df["start_date"]).dt.normalize()

    merged = wave_df.merge(step_df, on="end_date", how="inner")
    if merged.empty:
        raise ValueError("No matching intervals found between morphodynamics and wave power.")
    merged["period_days_safe"] = merged["period_days"].replace(0, np.nan)

    for buoy in sorted(merged["buoy"].unique()):
        buoy_df = merged[merged["buoy"] == buoy].sort_values("end_date")
        x = mdates.date2num(buoy_df["end_date"].to_list())
        if len(x) > 1:
            spacing_days = float(np.nanmedian(np.diff(x)))
        else:
            spacing_days = 30.0
        bar_width = 0.6 * spacing_days

        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        for col, label in DEFAULT_METRICS.items():
            series = buoy_df[col]
            if col == "erosion_m3":
                series = -series
            ax_top.plot(buoy_df["end_date"], series, marker="o", label=label)
        ax_top.set_ylabel("Volume change [m^3]")
        ax_top.set_title(f"Morphodynamics vs Wave Power ({buoy})")
        ax_top.grid(True, alpha=0.3)
        ax_top.legend(fontsize=8)

        ax_bottom.bar(x, buoy_df["cum_wave_power_MWh_m"].to_numpy(), width=bar_width, color="#4477AA")
        ax_bottom.set_ylabel("Cumulative wave power [MWh/m]")
        ax_bottom.set_xlabel("Survey interval end date")
        ax_bottom.grid(True, alpha=0.3, axis="y")
        ax_bottom.xaxis_date()
        ax_bottom.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        fig.autofmt_xdate(rotation=45, ha="right")

        fig.tight_layout()
        out_path = stacked_plot_dir / f"morpho_vs_wave_power_{buoy}.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved stacked plot: {out_path}")

        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        for col, label in DEFAULT_METRICS.items():
            series = buoy_df[col] / buoy_df["period_days_safe"]
            if col == "erosion_m3":
                series = -series
            ax_top.plot(buoy_df["end_date"], series, marker="o", label=label)
        ax_top.set_ylabel("Volume change [m^3/day]")
        ax_top.set_title(f"Morphodynamics vs Wave Power per Day ({buoy})")
        ax_top.grid(True, alpha=0.3)
        ax_top.legend(fontsize=8)

        ax_bottom.bar(
            x,
            (buoy_df["cum_wave_power_MWh_m"] / buoy_df["period_days_safe"]).to_numpy(),
            width=bar_width,
            color="#4477AA",
        )
        ax_bottom.set_ylabel("Cumulative wave power [MWh/m/day]")
        ax_bottom.set_xlabel("Survey interval end date")
        ax_bottom.grid(True, alpha=0.3, axis="y")
        ax_bottom.xaxis_date()
        ax_bottom.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        fig.autofmt_xdate(rotation=45, ha="right")

        fig.tight_layout()
        out_path = stacked_plot_dir / f"morpho_vs_wave_power_per_day_{buoy}.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved stacked plot: {out_path}")

    return merged


def plot_correlation_panels(
    merged_df: pd.DataFrame,
    buoys_list: list[str],
    metrics_map: dict[str, str],
    x_map: dict[str, str],
    output_dir: Path,
    highlight_dates=None,
    outlier_dates=None,
    verbose: bool = False,
) -> None:
    """Create per-buoy correlation panels for each metric vs x-variable.

    Args:
        merged_df: DataFrame of merged wave power and morphodynamics.
        buoys_list: Ordered list of buoy names to plot.
        metrics_map: Mapping of y-variable column name to label.
        x_map: Mapping of x-variable column name to label.
        output_dir: Output directory for plots.
        highlight_dates: Optional list of dates to highlight.
        outlier_dates: Optional list of dates to treat as outliers.
    """
    for buoy in buoys_list:
        buoy_df = merged_df[merged_df["buoy"] == buoy].sort_values("end_date")
        highlight_mask = interval_highlight_mask(
            buoy_df["start_date"], buoy_df["end_date"], highlight_dates
        )
        outlier_mask = interval_highlight_mask(
            buoy_df["start_date"], buoy_df["end_date"], outlier_dates
        )
        for x_col, x_label in x_map.items():
            if x_col not in buoy_df.columns:
                continue
            for col, label in metrics_map.items():
                if col not in buoy_df.columns:
                    continue
                x = buoy_df[x_col].to_numpy()
                y = buoy_df[col].to_numpy()
                valid = np.isfinite(x) & np.isfinite(y)
                if valid.sum() < 2:
                    continue

                x_valid = x[valid]
                y_valid = y[valid]
                highlight_valid = highlight_mask & valid
                clean_mask = valid & ~outlier_mask
                if clean_mask.sum() >= 2:
                    x_clean = x[clean_mask]
                    y_clean = y[clean_mask]
                else:
                    x_clean = x_valid
                    y_clean = y_valid

                r_full = float(np.corrcoef(x_valid, y_valid)[0, 1])
                r_clean = float(np.corrcoef(x_clean, y_clean)[0, 1])
                slope_full, intercept_full = np.polyfit(x_valid, y_valid, 1)
                slope_clean, intercept_clean = np.polyfit(x_clean, y_clean, 1)
                x_fit = np.linspace(x_valid.min(), x_valid.max(), 100)
                y_fit_full = slope_full * x_fit + intercept_full
                y_fit_clean = slope_clean * x_fit + intercept_clean

                fig, ax = plt.subplots(figsize=(6.5, 5))
                ax.scatter(x_valid, y_valid, alpha=0.75, edgecolor="black")
                if highlight_valid.any():
                    ax.scatter(
                        x[highlight_valid],
                        y[highlight_valid],
                        s=140,
                        color="red",
                        edgecolor="black",
                        zorder=5,
                    )
                ax.plot(x_fit, y_fit_full, color="black", linewidth=1.5, label=f"r={r_full:.3f}")
                if outlier_dates:
                    ax.plot(
                        x_fit,
                        y_fit_clean,
                        color="black",
                        linewidth=1.5,
                        linestyle="--",
                        label=f"r_clean={r_clean:.3f}",
                    )
                ax.set_xlabel(x_label)
                ax.set_ylabel(f"{label}")
                title = f"{label} vs {x_label} ({buoy})"
                if outlier_dates:
                    title = f"{title} [r_clean={r_clean:.3f}]"
                ax.set_title(title)
                ax.grid(True, alpha=0.3)
                ax.legend()

                out_path = output_dir / f"corr_{buoy}_{x_col}_{col}.png"
                fig.tight_layout()
                fig.savefig(out_path, dpi=300, bbox_inches="tight")
                plt.close(fig)
                if verbose:
                    print(f"Saved correlation plot: {out_path}")


def plot_correlation_compilation(
    merged_df: pd.DataFrame,
    buoys_list: list[str],
    metrics_map: dict[str, str],
    x_col: str,
    x_label: str,
    output_dir: Path,
    highlight_dates=None,
    outlier_dates=None,
    filename: str | None = None,
    verbose: bool = False,
) -> None:
    """Create a grid of correlations for a single x-variable.

    Args:
        merged_df: DataFrame of merged wave power and morphodynamics.
        buoys_list: Ordered list of buoy names to plot.
        metrics_map: Mapping of y-variable column name to label.
        x_col: Column name to use on the x-axis.
        x_label: Axis label for the x-axis.
        output_dir: Output directory for plots.
        highlight_dates: Optional list of dates to highlight.
        outlier_dates: Optional list of dates to treat as outliers.
        filename: Optional output filename override.
    """
    nrows = max(1, len(buoys_list))
    ncols = max(1, len(metrics_map))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 3.4 * nrows), sharex=False, sharey=False)
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    elif ncols == 1:
        axes = axes.reshape(-1, 1)

    for row_idx, buoy in enumerate(buoys_list):
        buoy_df = merged_df[merged_df["buoy"] == buoy].sort_values("end_date")
        highlight_mask = interval_highlight_mask(
            buoy_df["start_date"], buoy_df["end_date"], highlight_dates
        )
        outlier_mask = interval_highlight_mask(
            buoy_df["start_date"], buoy_df["end_date"], outlier_dates
        )
        for col_idx, (col, label) in enumerate(metrics_map.items()):
            if x_col not in buoy_df.columns or col not in buoy_df.columns:
                axes[row_idx, col_idx].set_visible(False)
                continue
            ax = axes[row_idx, col_idx]
            x = buoy_df[x_col].to_numpy()
            y = buoy_df[col].to_numpy()
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 2:
                ax.set_visible(False)
                continue

            x_valid = x[valid]
            y_valid = y[valid]
            highlight_valid = highlight_mask & valid
            clean_mask = valid & ~outlier_mask
            if clean_mask.sum() >= 2:
                x_clean = x[clean_mask]
                y_clean = y[clean_mask]
            else:
                x_clean = x_valid
                y_clean = y_valid

            r_clean = float(np.corrcoef(x_clean, y_clean)[0, 1])
            slope_full, intercept_full = np.polyfit(x_valid, y_valid, 1)
            slope_clean, intercept_clean = np.polyfit(x_clean, y_clean, 1)
            x_fit = np.linspace(x_valid.min(), x_valid.max(), 100)
            y_fit_full = slope_full * x_fit + intercept_full
            y_fit_clean = slope_clean * x_fit + intercept_clean

            ax.scatter(x_valid, y_valid, alpha=0.75, edgecolor="black", s=20)
            if highlight_valid.any():
                ax.scatter(
                    x[highlight_valid],
                    y[highlight_valid],
                    s=140,
                    color="red",
                    edgecolor="black",
                    zorder=5,
                )
            ax.plot(x_fit, y_fit_full, color="black", linewidth=1.2)
            if outlier_dates:
                ax.plot(x_fit, y_fit_clean, color="black", linewidth=1.2, linestyle="--")
            title = f"{label}\nr={r_clean:.3f}"
            if outlier_dates:
                title = f"{label}\nr_clean={r_clean:.3f}"
            ax.set_title(title, fontsize=9)
            if row_idx == nrows - 1:
                ax.set_xlabel(x_label)
            if col_idx == 0:
                ax.set_ylabel(f"{buoy}\n{label}")
            ax.grid(True, alpha=0.3)

    fig.suptitle(f"Morphodynamics vs {x_label}", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    comp_name = filename or f"correlation_compilation_{x_col}.png"
    comp_path = output_dir / comp_name
    fig.savefig(comp_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    if verbose:
        print(f"Saved compilation plot: {comp_path}")


def run_morpho_wave_correlations(
    merged_df: pd.DataFrame,
    output_dir: str | Path,
    highlight_dates=None,
    outlier_dates=None,
    metrics_map: dict[str, str] | None = None,
    x_map: dict[str, str] | None = None,
    verbose: bool = False,
) -> None:
    """Generate correlation plots for morphodynamics vs wave power metrics.

    Args:
        merged_df: DataFrame from ``run_morphodynamics_vs_wave_power``.
        output_dir: Output directory for plots.
        highlight_dates: Optional list of dates to highlight.
        outlier_dates: Optional list of dates to treat as outliers.
        metrics_map: Mapping of y-variable column name to label.
        x_map: Mapping of x-variable column name to label.
    """
    merged_df = merged_df.copy()
    merged_df["period_days_safe"] = merged_df["period_days"].replace(0, np.nan)
    metrics_map = metrics_map or DEFAULT_METRICS
    x_map = x_map or DEFAULT_X_VARS

    for col in metrics_map.keys():
        merged_df[f"{col}_per_day"] = merged_df[col] / merged_df["period_days_safe"]

    power_cols = [
        "cum_wave_power_MWh_m",
        "cum_wave_power_above_MWh_m",
        "cum_wave_power_below_MWh_m",
    ]
    for col in power_cols:
        if col in merged_df.columns:
            merged_df[f"{col}_per_day"] = merged_df[col] / merged_df["period_days_safe"]

    metrics_rate = {
        f"{col}_per_day": f"{label} per day [m^3/day]"
        for col, label in metrics_map.items()
    }
    x_vars_rate = {
        f"{col}_per_day": label.replace("[MWh/m]", "[MWh/m/day]")
        for col, label in x_map.items()
        if col != "period_days" and f"{col}_per_day" in merged_df.columns
    }

    buoys = sorted(merged_df["buoy"].unique())
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_correlation_panels(
        merged_df,
        buoys,
        metrics_map,
        x_map,
        output_dir,
        highlight_dates=highlight_dates,
        outlier_dates=outlier_dates,
        verbose=verbose,
    )
    for x_col, x_label in x_map.items():
        plot_correlation_compilation(
            merged_df,
            buoys,
            metrics_map,
            x_col,
            x_label,
            output_dir,
            highlight_dates=highlight_dates,
            outlier_dates=outlier_dates,
            verbose=verbose,
        )
    plot_correlation_compilation(
        merged_df,
        buoys,
        metrics_map,
        "period_days",
        "Interval duration [days]",
        output_dir,
        highlight_dates=highlight_dates,
        outlier_dates=outlier_dates,
        filename="correlation_compilation_interval_duration.png",
        verbose=verbose,
    )

    plot_correlation_panels(
        merged_df,
        buoys,
        metrics_rate,
        x_vars_rate,
        output_dir,
        highlight_dates=highlight_dates,
        outlier_dates=outlier_dates,
        verbose=verbose,
    )
    for x_col, x_label in x_vars_rate.items():
        plot_correlation_compilation(
            merged_df,
            buoys,
            metrics_rate,
            x_col,
            x_label,
            output_dir,
            highlight_dates=highlight_dates,
            outlier_dates=outlier_dates,
            filename=f"correlation_compilation_rate_{x_col}.png",
            verbose=verbose,
        )
