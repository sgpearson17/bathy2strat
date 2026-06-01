"""Stratigraphy computation and plotting utilities for bathymetry cubes."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Colormap
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from pyproj import CRS, Transformer
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat

from .bathy import BathyGrid, load_bathy_grid
from .plot_style import apply_axes_font as _plot_apply_axes_font
from .plot_style import get_plot_font as _plot_get_plot_font


MHW = 0.358
MSL = -0.112
MLW = -0.590


def _get_plot_font() -> FontProperties:
    return _plot_get_plot_font()


def _apply_axes_font(ax, font: FontProperties) -> None:
    _plot_apply_axes_font(ax, font)


@dataclass
class BathyCube:
    """Container for gridded bathymetry time series.

    :ivar location: Descriptive location name or identifier.
    :ivar t: MATLAB datenum survey times, shape (nt,).
    :ivar x: 2D grid of x coordinates in meters.
    :ivar y: 2D grid of y coordinates in meters.
    :ivar z: 3D bathymetry cube, shape (ny, nx, nt) or (nx, ny, nt).
    """

    location: str
    t: np.ndarray  # MATLAB datenum, shape (nt,)
    x: np.ndarray  # 2D grid, meters
    y: np.ndarray  # 2D grid, meters
    z: np.ndarray  # 3D cube, shape (ny, nx, nt) or (nx, ny, nt)


@dataclass
class Transect:
    """Transect polyline in projected coordinates.

    :ivar name: Transect label.
    :ivar x: X coordinates in meters.
    :ivar y: Y coordinates in meters.
    """

    name: str
    x: np.ndarray  # meters
    y: np.ndarray  # meters


@dataclass
class StratigraphyConfig:
    """Configuration values for stratigraphy processing and plotting.

    :ivar initial_index: Survey index used as the baseline surface.
    :ivar dx: Grid cell size in meters.
    :ivar target_crs: Target projected CRS for transects (e.g., EPSG:32618).
    :ivar dot_spacing_km: Spacing between map dots along transects.
    :ivar end_dot_size: Marker size for transect endpoints.
    :ivar mid_dot_size: Marker size for intermediate points.
    """

    initial_index: int = 0
    dx: float = 20.0
    target_crs: str = "EPSG:32618"
    dot_spacing_km: float = 0.1
    end_dot_size: float = 40.0
    mid_dot_size: float = 15.0


@dataclass
class StratigraphyResult:
    """Computed stratigraphy outputs and derived metrics.

    :ivar t: Survey datenums.
    :ivar t_full: Monthly datenums spanning the survey years.
    :ivar min_surf: Minimum surface after the baseline index.
    :ivar max_surf: Maximum surface after the baseline index.
    :ivar deposit_elev: Erosion-adjusted deposit elevations.
    :ivar deposit_thk_full: Monthly deposit thickness stack.
    :ivar deposit_per_year: Annual deposit volumes per survey year.
    :ivar total_sed_vol_per_year: Total sediment volume per survey year.
    :ivar theseus_ratio: Ratio of preserved to original deposits.
    """

    t: np.ndarray
    t_full: np.ndarray
    min_surf: np.ndarray
    max_surf: np.ndarray
    deposit_elev: np.ndarray
    deposit_thk_full: np.ndarray
    deposit_per_year: np.ndarray
    total_sed_vol_per_year: np.ndarray
    theseus_ratio: np.ndarray


def datenum_to_datetime64(datenum: np.ndarray) -> np.ndarray:
    """Convert MATLAB datenum values to numpy datetime64.

    :param datenum: Array-like MATLAB datenums.
    :returns: 1D array of numpy datetime64 values.
    """
    import datetime as dt

    out = []
    for value in np.asarray(datenum, dtype=float).reshape(-1):
        ordinal = int(value)
        frac = value - ordinal
        py_dt = dt.datetime.fromordinal(ordinal) + dt.timedelta(days=frac) - dt.timedelta(days=366)
        out.append(np.datetime64(py_dt))
    return np.asarray(out)


def datetime64_to_datenum(values: Iterable[np.datetime64]) -> np.ndarray:
    """Convert datetime64 values to MATLAB datenum floats.

    :param values: Iterable of numpy datetime64 values.
    :returns: Array of MATLAB datenums as floats.
    """
    import datetime as dt

    out = []
    for value in values:
        py_dt = pd.Timestamp(value).to_pydatetime()
        out.append(py_dt.toordinal() + 366 + (py_dt.hour / 24.0) + (py_dt.minute / 1440.0) + (py_dt.second / 86400.0))
    return np.asarray(out, dtype=float)


def _datenum_to_years(t_vals: np.ndarray) -> np.ndarray:
    import datetime as dt

    years = []
    for value in np.asarray(t_vals, dtype=float).reshape(-1):
        ordinal = int(value)
        frac = value - ordinal
        py_dt = dt.datetime.fromordinal(ordinal) + dt.timedelta(days=frac) - dt.timedelta(days=366)
        years.append(py_dt.year + (py_dt.timetuple().tm_yday - 1) / 365.0)
    return np.asarray(years, dtype=float)


def _resolve_time_axis(t_vals: np.ndarray | None, nt: int, start_at_zero: bool) -> tuple[np.ndarray, str]:
    if t_vals is None:
        time = np.arange(nt, dtype=float)
        label = "Years since start" if start_at_zero else "Year"
        return time, label

    years = _datenum_to_years(t_vals)
    if start_at_zero:
        return years - years[0], "Years since start"
    return years, "Year"


def compute_deposit_elev_1d(z_stack: np.ndarray) -> np.ndarray:
    """Apply the erosion rule to build deposit elevations for 1D profiles."""
    deposit_elev = z_stack.copy()
    nt = z_stack.shape[0]
    for tt in range(1, nt):
        dz = z_stack[tt, :] - z_stack[tt - 1, :]
        deposit_elev[tt, :] = z_stack[tt, :]
        erosion = dz < 0
        if np.any(erosion):
            for qq in range(0, tt):
                prev = deposit_elev[qq, :]
                update = erosion & (deposit_elev[tt, :] < prev)
                prev[update] = deposit_elev[tt, :][update]
                deposit_elev[qq, :] = prev
    return deposit_elev


def compute_remaining_volumes(z_stack: np.ndarray, dx: float, base: float) -> tuple[np.ndarray, float]:
    """Compute remaining deposit volumes per deposit year through time (unit width)."""
    nt = z_stack.shape[0]
    remaining = np.zeros((nt, nt), dtype=float)
    initial_total = None
    for t in range(nt):
        dep = compute_deposit_elev_1d(z_stack[: t + 1, :])
        thickness = np.zeros((t + 1, z_stack.shape[1]), dtype=float)
        thickness[0, :] = np.maximum(0.0, dep[0, :] - base)
        for k in range(1, t + 1):
            thickness[k, :] = np.maximum(0.0, dep[k, :] - dep[k - 1, :])
        vols = np.nansum(thickness, axis=1) * dx
        remaining[t, : t + 1] = vols
        if t == 0:
            initial_total = np.nansum(vols)
    if initial_total is None or initial_total == 0:
        initial_total = 1.0
    return remaining, initial_total


def plot_stratigraphy_stack(
    cube: BathyCube,
    result: StratigraphyResult,
    out_path: Path,
    scenario: str,
    time_vals: np.ndarray | None = None,
    start_year_at_zero: bool = True,
    plot_xlim: tuple[float, float] | None = None,
    plot_ylim: tuple[float, float] | None = None,
    highlight_date: str | np.datetime64 | None = None,
    title_prefix: str | None = None,
) -> None:
    """Plot stacked stratigraphy plus preservation metrics for a 1D transect cube."""
    x_m = cube.x[0, :]
    z_stack = cube.z[0, :, :].T
    deposit_elev = compute_deposit_elev_1d(z_stack)
    max_surface = np.nanmax(z_stack, axis=0)
    nt = z_stack.shape[0]
    colors = plt.cm.viridis(np.linspace(0.15, 0.95, nt))
    dx = float(x_m[1] - x_m[0]) if x_m.size > 1 else 1.0
    font = _get_plot_font()

    highlight_layer = None
    highlight_surface = None
    highlight_tag = None
    if highlight_date is not None and cube.t is not None and len(cube.t) > 0:
        t_dt = datenum_to_datetime64(cube.t)
        highlight_dt = np.datetime64(highlight_date)
        idx = np.searchsorted(t_dt, highlight_dt, side="left")
        if idx >= len(t_dt):
            idx = len(t_dt) - 1
        if idx >= 0:
            highlight_layer = int(idx)
            highlight_tag = str(t_dt[highlight_layer])[:10]
            intact = np.isfinite(z_stack[highlight_layer, :]) & np.isfinite(deposit_elev[highlight_layer, :])
            intact &= np.isclose(deposit_elev[highlight_layer, :], z_stack[highlight_layer, :], atol=1e-6)
            highlight_surface = np.where(intact, deposit_elev[highlight_layer, :], np.nan)

    fig = plt.figure(figsize=(14.5, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.7], hspace=0.45, wspace=0.25)
    axes_top = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]
    axes_bottom = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])]

    # Compose subplot titles with optional section names and panel labels.
    def _title(prefix: str, base: str) -> str:
        if title_prefix:
            return f"{prefix} {title_prefix} - {base}"
        return f"{prefix} {base}"

    ax = axes_top[0]
    for tt in range(nt):
        ax.plot(x_m, z_stack[tt, :], color=colors[tt], linewidth=0.9)
    if highlight_layer is not None:
        ax.plot(x_m, z_stack[highlight_layer, :], color="red", linewidth=1.6, zorder=3)
        if highlight_tag:
            legend_label = f"{highlight_tag}"
            legend_handle = Line2D([], [], color="red", linewidth=1.6, label=legend_label)
            legend = ax.legend(
                handles=[legend_handle],
                loc="upper right",
                frameon=True,
                facecolor="white",
                edgecolor="none",
                framealpha=0.5,
            )
            for text in legend.get_texts():
                text.set_color("red")
    ax.set_title(_title("(a)", "Raw surfaces"), fontproperties=font)
    ax.set_xlabel("Distance [m]", fontproperties=font)
    ax.set_ylabel("Elevation [m]", fontproperties=font)
    ax.grid(True, alpha=0.3)

    ax = axes_top[1]
    for tt in range(nt):
        ax.plot(x_m, deposit_elev[tt, :], color=colors[tt], linewidth=0.9)
    if highlight_surface is not None:
        ax.plot(x_m, highlight_surface, color="red", linewidth=1.6, zorder=3)
    ax.set_title(_title("(b)", "After erosion rule"), fontproperties=font)
    ax.set_xlabel("Distance [m]", fontproperties=font)
    ax.grid(True, alpha=0.3)

    ax = axes_top[2]
    min_elev = np.nanmin([np.nanmin(z_stack), np.nanmin(deposit_elev)])
    max_elev = np.nanmax([np.nanmax(z_stack), np.nanmax(deposit_elev)])
    base = min_elev - 0.01 * abs(min_elev)
    cumulative = base + np.zeros_like(x_m)
    ax.fill_between(x_m, base, deposit_elev[0, :], color=colors[0], alpha=1.0)
    cumulative = deposit_elev[0, :]
    for tt in range(1, nt):
        layer = np.maximum(0.0, deposit_elev[tt, :] - deposit_elev[tt - 1, :])
        upper = cumulative + layer
        layer_color = "red" if highlight_layer == tt else colors[tt]
        ax.fill_between(x_m, cumulative, upper, color=layer_color, alpha=1.0)
        cumulative = upper
    for tt in range(nt):
        ax.plot(x_m, deposit_elev[tt, :], color="k", linewidth=0.6, alpha=0.8)
    if highlight_surface is not None:
        ax.plot(x_m, highlight_surface, color="red", linewidth=1.6, zorder=3)
    if np.any(np.isfinite(max_surface)):
        ax.plot(x_m, max_surface, color="0.2", linestyle="--", linewidth=0.8, zorder=2)
    ax.plot(x_m, deposit_elev[-1, :], color="k", linewidth=2.0)
    ax.set_title(_title("(c)", "Stacked stratigraphy"), fontproperties=font)
    ax.set_xlabel("Distance [m]", fontproperties=font)
    ax.grid(True, alpha=0.3)

    y_min = base
    y_max = max_elev + 0.01 * abs(max_elev)
    for ax in axes_top:
        ax.set_ylim(y_min, y_max)

    show_water_labels = True

    if scenario == "bruun_slr":
        if plot_xlim is None:
            plot_xlim = (-300, 800)
        if plot_ylim is None:
            plot_ylim = (-5, 3)

    if plot_xlim is not None:
        for ax in axes_top:
            ax.set_xlim(plot_xlim)
    if plot_ylim is not None:
        for ax in axes_top:
            ax.set_ylim(plot_ylim)

    remaining, initial_total = compute_remaining_volumes(z_stack, dx=dx, base=base)
    time_axis, time_label = _resolve_time_axis(time_vals, nt, start_year_at_zero)
    for k in range(1, nt):
        y = remaining[:, k].astype(float)
        y[:k] = np.nan
        axes_bottom[0].plot(time_axis, y, linewidth=1.2, color=colors[k])
        denom = remaining[k, k] if np.isfinite(remaining[k, k]) else 0.0
        if denom == 0.0:
            y_norm = np.zeros_like(y)
            y_norm[:k] = np.nan
        else:
            y_norm = y / denom
        axes_bottom[1].plot(time_axis, y_norm, linewidth=1.2, color=colors[k])
    if highlight_layer is not None and highlight_layer >= 1:
        y = remaining[:, highlight_layer].astype(float)
        y[:highlight_layer] = np.nan
        axes_bottom[0].plot(time_axis, y, linewidth=2.4, color="red", zorder=3)
        denom = remaining[highlight_layer, highlight_layer] if np.isfinite(remaining[highlight_layer, highlight_layer]) else 0.0
        if denom == 0.0:
            y_norm = np.zeros_like(y)
            y_norm[:highlight_layer] = np.nan
        else:
            y_norm = y / denom
        axes_bottom[1].plot(time_axis, y_norm, linewidth=2.4, color="red", zorder=3)
    axes_bottom[0].set_title(_title("(d)", "Volume preserved (absolute)"), fontproperties=font)
    axes_bottom[0].set_xlabel(time_label, fontproperties=font)
    axes_bottom[0].set_ylabel("Volume (unit width)", fontproperties=font)
    axes_bottom[0].grid(True, alpha=0.3)
    axes_bottom[1].set_title(_title("(e)", "Volume preserved (normalized)"), fontproperties=font)
    axes_bottom[1].set_xlabel(time_label, fontproperties=font)
    axes_bottom[1].set_ylabel("Fraction of initial", fontproperties=font)
    axes_bottom[1].grid(True, alpha=0.3)

    denom0 = remaining[0, 0] if np.isfinite(remaining[0, 0]) and remaining[0, 0] != 0 else np.nan
    ratio = np.clip(remaining[:, 0] / denom0, 0.0, 1.0)
    axes_bottom[2].plot(time_axis, ratio, color="k", linewidth=1.4)
    if highlight_layer is not None and highlight_layer < len(time_axis):
        axes_bottom[2].axvline(time_axis[highlight_layer], color="red", linewidth=1.6, zorder=3)
    axes_bottom[2].set_title(_title("(f)", "Theseus ratio (t0 preserved)"), fontproperties=font)
    axes_bottom[2].set_xlabel(time_label, fontproperties=font)
    axes_bottom[2].set_ylabel("Fraction of initial", fontproperties=font)
    ratio_min = np.nanmin(ratio)
    if np.isfinite(ratio_min):
        axes_bottom[2].set_ylim(ratio_min, 1.0)
    else:
        axes_bottom[2].set_ylim(0.0, 1.0)
    time_max = np.nanmax(time_axis)
    if np.isfinite(time_max):
        if start_year_at_zero:
            axes_bottom[2].set_xlim(0.0, time_max)
        else:
            time_min = np.nanmin(time_axis)
            if np.isfinite(time_min):
                axes_bottom[2].set_xlim(time_min, time_max)
    axes_bottom[2].grid(True, alpha=0.3)

    # if show_water_labels:
    #     x_text = x_m[0] if len(x_m) else 0.0
    #     offset = 0.01 * (y_max - y_min) if np.isfinite(y_max - y_min) and y_max > y_min else 0.15
    #     axes_top[1].text(x_text, MHW + offset, "MHW", fontsize=9, zorder=0, fontproperties=font)
    #     axes_top[1].text(x_text, MLW + offset, "MLW", fontsize=9, zorder=0, fontproperties=font)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    for ax in axes_top + axes_bottom:
        _apply_axes_font(ax, font)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def plot_stacked_stratigraphy_section(
    cube: BathyCube,
    out_path: Path,
    title: str,
    plot_xlim: tuple[float, float] | None = None,
    plot_ylim: tuple[float, float] | None = None,
    x_ticks: np.ndarray | None = None,
    mhw: float | None = None,
    mlw: float | None = None,
    highlight_date: str | np.datetime64 | None = None,
    max_transect_length: float | None = None,
    max_depth_range: float | None = None,
    max_fig_width_cm: float = 20.0,
    max_fig_height_cm: float = 5.0,
    scale_factor: float = 2.0,
    clip_x_to_data: bool = False,
    clip_y_to_data: bool = False,
) -> None:
    """Plot a standalone stacked stratigraphy section with a fixed title.

    :param cube: 1D transect cube.
    :param out_path: Output image path.
    :param title: Plot title.
    :param plot_xlim: Optional x-axis limits in meters.
    :param plot_ylim: Optional y-axis limits in meters.
    :param x_ticks: Optional x-axis tick positions in meters.
    :param mhw: Optional mean high water elevation.
    :param mlw: Optional mean low water elevation.
    :param highlight_date: Optional date string or datetime64 to highlight a deposit.
    :param max_transect_length: Max transect length for scaling (meters).
    :param max_depth_range: Max depth range for scaling (meters).
    :param max_fig_width_cm: Max figure width for scaling (cm).
    :param max_fig_height_cm: Max figure height for scaling (cm).
    :param scale_factor: Multiplier for scaled figure size.
    :param clip_x_to_data: Clip x-limits to the transect data extent.
    :param clip_y_to_data: Clip y-limits to the data extent.
    """
    x_m = cube.x[0, :]
    z_stack = cube.z[0, :, :].T
    deposit_elev = compute_deposit_elev_1d(z_stack)
    max_surface = np.nanmax(z_stack, axis=0)
    nt = z_stack.shape[0]
    colors = plt.cm.viridis(np.linspace(0.15, 0.95, nt))
    font = _get_plot_font()

    min_elev = np.nanmin([np.nanmin(z_stack), np.nanmin(deposit_elev)])
    max_elev = np.nanmax([np.nanmax(z_stack), np.nanmax(deposit_elev)])
    base = min_elev - 0.01 * abs(min_elev)

    highlight_layer = None
    highlight_surface = None
    highlight_tag = None
    if highlight_date is not None and cube.t is not None and len(cube.t) > 0:
        t_dt = datenum_to_datetime64(cube.t)
        highlight_dt = np.datetime64(highlight_date)
        idx = np.searchsorted(t_dt, highlight_dt, side="left")
        if idx >= len(t_dt):
            idx = len(t_dt) - 1
        if idx >= 0:
            highlight_layer = int(idx)
            highlight_tag = str(t_dt[highlight_layer])[:10]
            intact = np.isfinite(z_stack[highlight_layer, :]) & np.isfinite(deposit_elev[highlight_layer, :])
            intact &= np.isclose(deposit_elev[highlight_layer, :], z_stack[highlight_layer, :], atol=1e-6)
            highlight_surface = np.where(intact, deposit_elev[highlight_layer, :], np.nan)

    data_xmin = float(np.nanmin(x_m)) if len(x_m) else 0.0
    data_xmax = float(np.nanmax(x_m)) if len(x_m) else 0.0
    data_ymin = float(np.nanmin([np.nanmin(z_stack), np.nanmin(deposit_elev)]))
    data_ymax = float(np.nanmax([np.nanmax(z_stack), np.nanmax(deposit_elev)]))

    if clip_x_to_data:
        x_range = data_xmax - data_xmin
    elif plot_xlim is not None:
        x_range = float(plot_xlim[1] - plot_xlim[0])
    else:
        x_range = data_xmax - data_xmin

    if clip_y_to_data:
        y_range = data_ymax - data_ymin
    elif plot_ylim is not None:
        y_range = float(plot_ylim[1] - plot_ylim[0])
    else:
        y_range = float(max_elev - min_elev) if np.isfinite(max_elev) and np.isfinite(min_elev) else 0.0
    if max_transect_length and max_depth_range and max_transect_length > 0 and max_depth_range > 0:
        max_fig_width_in = max_fig_width_cm / 2.54
        max_fig_height_in = max_fig_height_cm / 2.54
        fig_width = (x_range / max_transect_length) * max_fig_width_in * scale_factor
        fig_height = (y_range / max_depth_range) * max_fig_height_in * scale_factor
        fig_width = max(fig_width, 3.0)
        fig_height = max(fig_height, 2.0)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    else:
        fig, ax = plt.subplots(figsize=(14.0, 4.0))

    mhw_val = MHW if mhw is None else mhw
    mlw_val = MLW if mlw is None else mlw
    mhw_color = "#0b2e5b"
    dash_style = (0, (6, 3))
    ax.axhline(mhw_val, linestyle=dash_style, color=mhw_color, linewidth=0.5, zorder=0)
    ax.axhline(mlw_val, linestyle=dash_style, color=mhw_color, linewidth=0.5, zorder=0)
    cumulative = base + np.zeros_like(x_m)
    ax.fill_between(x_m, base, deposit_elev[0, :], color=colors[0], alpha=1.0)
    cumulative = deposit_elev[0, :]
    for tt in range(1, nt):
        layer = np.maximum(0.0, deposit_elev[tt, :] - deposit_elev[tt - 1, :])
        upper = cumulative + layer
        layer_color = "red" if highlight_layer == tt else colors[tt]
        ax.fill_between(x_m, cumulative, upper, color=layer_color, alpha=1.0, zorder=1)
        cumulative = upper
    for tt in range(nt):
        ax.plot(x_m, deposit_elev[tt, :], color="k", linewidth=0.6, alpha=0.8)
    if highlight_surface is not None:
        ax.plot(x_m, highlight_surface, color="red", linewidth=1.6, zorder=3)
    if np.any(np.isfinite(max_surface)):
        ax.plot(x_m, max_surface, color="0.2", linestyle="--", linewidth=0.8, zorder=2)
    ax.plot(x_m, deposit_elev[-1, :], color="k", linewidth=2.0)

    ax.set_title(title, fontproperties=font)
    if highlight_tag:
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
    ax.set_xlabel("Distance [m]", fontproperties=font)
    ax.set_ylabel("Elevation [m]", fontproperties=font)
    ax.grid(True, alpha=0.3)

    if clip_y_to_data:
        y_min_plot, y_max_plot = data_ymin, data_ymax
    elif plot_ylim is not None:
        y_min_plot, y_max_plot = plot_ylim
    else:
        y_min_plot, y_max_plot = min_elev, max_elev
    y_range = y_max_plot - y_min_plot
    if not np.isfinite(y_range) or y_range <= 0:
        y_range = 1.0
    fig.canvas.draw()
    font_size = font.get_size_in_points()
    if not font_size:
        font_size = float(plt.rcParams.get("font.size", 10.0))
    text_px = font_size * fig.dpi / 72.0
    text_px += 1.1 * text_px
    bbox_height_px = ax.bbox.height if ax.bbox.height > 0 else 1.0
    data_per_px = y_range / bbox_height_px
    mhw_margin = text_px * data_per_px
    y_max_plot = max(y_max_plot, mhw_val + mhw_margin)
    y_range = y_max_plot - y_min_plot
    offset = 0.01 * y_range if np.isfinite(y_range) and y_range > 0 else 0.15
    x_text = x_m[0] if len(x_m) else 0.0
    if np.isfinite(y_max_plot) and y_max_plot >= mlw_val:
        ax.text(x_text, mhw_val + offset, "MHW", fontsize=9, zorder=0, fontproperties=font, color=mhw_color)
        ax.text(x_text, mlw_val + offset, "MLW", fontsize=9, zorder=0, fontproperties=font, color=mhw_color)

    if clip_x_to_data:
        ax.set_xlim(data_xmin, data_xmax)
    elif plot_xlim is not None:
        ax.set_xlim(plot_xlim)
    else:
        ax.set_xlim(0.0, float(x_m[-1]) if len(x_m) else 0.0)
    if clip_y_to_data or plot_ylim is not None:
        ax.set_ylim(y_min_plot, y_max_plot)
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    _apply_axes_font(ax, font)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def _coerce_mat_struct(value):
    """Normalize MATLAB structs loaded as object arrays.

    :param value: Value returned from ``scipy.io.loadmat``.
    :returns: Unwrapped struct if stored in a 0-d object array.
    """
    if isinstance(value, np.ndarray) and value.dtype == object:
        return value.flat[0]
    return value


def load_bathy_mat_known_structure(path: str | Path) -> BathyCube:
    """Load a bathy MAT file with the expected ``bathy`` struct fields.

    :param path: Path to the MAT file.
    :returns: Parsed :class:`BathyCube`.
    :raises ValueError: If the MAT file does not match the expected structure.
    """
    mat = loadmat(str(path), squeeze_me=True, struct_as_record=False)
    if "bathy" not in mat:
        raise ValueError(f"Expected a 'bathy' struct in {path}")

    b = _coerce_mat_struct(mat["bathy"])
    location = str(getattr(b, "location", Path(path).stem))
    t = np.asarray(getattr(b, "t"), dtype=float).reshape(-1)
    x = np.asarray(getattr(b, "x"), dtype=float)
    y = np.asarray(getattr(b, "y"), dtype=float)
    z = np.asarray(getattr(b, "z"), dtype=float)

    if x.shape != y.shape:
        raise ValueError(f"x and y shapes differ: {x.shape} vs {y.shape}")
    if z.ndim != 3:
        raise ValueError(f"Expected 3D z array. Got shape: {z.shape}")

    # Align z first two dimensions to x/y if needed.
    if z.shape[:2] != x.shape and z.shape[:2] == x.T.shape:
        z = np.transpose(z, (1, 0, 2))
    if z.shape[:2] != x.shape:
        raise ValueError(f"z spatial shape {z.shape[:2]} does not match x/y shape {x.shape}")

    return BathyCube(location=location, t=t, x=x, y=y, z=z)


def _bathygrid_to_cube(grid: BathyGrid) -> BathyCube:
    """Convert a :class:`BathyGrid` to a :class:`BathyCube`.

    :param grid: Regridded bathymetry grid.
    :returns: BathyCube representation.
    """
    t = np.asarray(grid.t, dtype=float).reshape(-1)
    return BathyCube(location=grid.location, t=t, x=grid.x, y=grid.y, z=grid.z)


def load_bathy_cube(path: str | Path) -> BathyCube:
    """Load bathymetry from .mat or .nc into a :class:`BathyCube`.

    :param path: Path to bathymetry file.
    :returns: BathyCube instance.
    :raises ValueError: If the file extension is unsupported.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".mat":
        return load_bathy_mat_known_structure(path)
    if suffix == ".nc":
        grid = load_bathy_grid(path)
        return _bathygrid_to_cube(grid)
    raise ValueError(f"Unsupported bathy file extension: {path.suffix}. Use .mat or .nc")


def _derive_xy_axes(x2d: np.ndarray, y2d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Extract monotonic 1D x/y axes from 2D mesh grids.

    :param x2d: 2D mesh grid of x coordinates.
    :param y2d: 2D mesh grid of y coordinates.
    :returns: Tuple of ``(x_axis, y_axis)``.
    :raises ValueError: If axes are not monotonic.
    """
    x_axis = np.asarray(x2d[0, :], dtype=float)
    y_axis = np.asarray(y2d[:, 0], dtype=float)

    if not (np.all(np.diff(x_axis) > 0) and np.all(np.diff(y_axis) > 0)):
        raise ValueError("Only monotonic regular grids are currently supported for transect interpolation")

    return x_axis, y_axis


def _nanminmax_surface(stack: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute nan-safe min/max surfaces along an axis.

    :param stack: Input array.
    :param axis: Axis along which to compute min/max.
    :returns: Tuple of ``(min_surface, max_surface)``.
    """
    all_nan = np.all(np.isnan(stack), axis=axis)
    stack_safe = np.where(all_nan[..., None], np.inf, stack)
    min_surf = np.nanmin(stack_safe, axis=axis)
    min_surf[all_nan] = np.nan

    stack_safe = np.where(all_nan[..., None], -np.inf, stack)
    max_surf = np.nanmax(stack_safe, axis=axis)
    max_surf[all_nan] = np.nan
    return min_surf, max_surf


def compute_stratigraphy(bathy: BathyCube, config: StratigraphyConfig) -> StratigraphyResult:
    """Compute stratigraphic thickness and volume metrics from a bathy cube.

    :param bathy: Bathymetry cube with surveys over time.
    :param config: Stratigraphy configuration.
    :returns: Stratigraphy results and derived metrics.
    :raises ValueError: If ``initial_index`` is out of bounds.
    """
    initial_idx = config.initial_index
    x = bathy.x
    y = bathy.y
    z = np.asarray(bathy.z, dtype=float).copy()
    t = np.asarray(bathy.t, dtype=float).reshape(-1)

    if initial_idx < 0 or initial_idx >= z.shape[2]:
        raise ValueError(f"initial_index out of bounds: {initial_idx}")

    # If any survey at a cell is NaN, keep it NaN across all surveys.
    any_nan = np.any(np.isnan(z), axis=2)
    z[any_nan] = np.nan

    min_surf, max_surf = _nanminmax_surface(z[:, :, initial_idx:], axis=2)

    deposit_elev = np.full_like(z, np.nan)
    deposit_elev[:, :, initial_idx] = z[:, :, initial_idx]

    nt = z.shape[2]
    deposit_per_year = np.zeros((nt, nt), dtype=float)
    t_full_global = np.array([], dtype=float)
    deposit_thk_full_global = np.zeros((z.shape[0], z.shape[1], 0), dtype=float)

    for t_deposit_year in range(initial_idx + 1, nt):
        for tt in range(initial_idx + 1, t_deposit_year + 1):
            dz = z[:, :, tt] - z[:, :, tt - 1]
            deposit_elev[:, :, tt] = z[:, :, tt]

            erosion_mask = dz < 0
            if np.any(erosion_mask):
                cur = deposit_elev[:, :, tt]
                for k in range(initial_idx, tt):
                    prev = deposit_elev[:, :, k]
                    # Push erosion downward to earlier preserved surfaces.
                    update = erosion_mask & ~np.isnan(prev) & (cur < prev)
                    prev[update] = cur[update]
                    deposit_elev[:, :, k] = prev

        deposit_thk = np.diff(deposit_elev, axis=2)
        if initial_idx > 0:
            deposit_thk = deposit_thk[:, :, initial_idx - 1 :]
        deposit_thk = np.concatenate([min_surf[:, :, None], deposit_thk], axis=2)

        survey_datetimes = datenum_to_datetime64(t)
        start_year = pd.Timestamp(survey_datetimes[0]).year
        end_year = pd.Timestamp(survey_datetimes[-1]).year

        t_full = []
        for year in range(start_year, end_year + 1):
            for month in range(1, 13):
                t_full.append(np.datetime64(f"{year:04d}-{month:02d}-01"))
        t_full = np.asarray(t_full)
        t_full_dn = datetime64_to_datenum(t_full)

        deposit_thk_full = np.zeros((x.shape[0], x.shape[1], len(t_full_dn)), dtype=float)

        survey_count = 0
        for tt in range(len(t_full_dn) - 1):
            if survey_count >= len(t):
                break
            if t_full_dn[tt] <= t[survey_count] <= t_full_dn[tt + 1]:
                if survey_count < deposit_thk.shape[2]:
                    deposit_thk_full[:, :, tt] = deposit_thk[:, :, survey_count]
                if survey_count >= len(t) - 1:
                    break
                survey_count += 1

        deposit_per_year[:, t_deposit_year] = np.nansum(deposit_thk * (config.dx**2), axis=(0, 1))

        t_full_global = t_full_dn
        deposit_thk_full_global = deposit_thk_full

    total_sed_vol_per_year, deposit_per_year, theseus_ratio = compute_preservation_potential(
        z=z,
        initial_idx=initial_idx,
        dx=config.dx,
        deposit_per_year=deposit_per_year,
    )

    return StratigraphyResult(
        t=t,
        t_full=t_full_global,
        min_surf=min_surf,
        max_surf=max_surf,
        deposit_elev=deposit_elev,
        deposit_thk_full=deposit_thk_full_global,
        deposit_per_year=deposit_per_year,
        total_sed_vol_per_year=total_sed_vol_per_year,
        theseus_ratio=theseus_ratio,
    )


def compute_preservation_potential(
    z: np.ndarray,
    initial_idx: int,
    dx: float,
    deposit_per_year: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute preservation potential and Theseus ratio metrics.

    :param z: Bathymetry stack, shape (ny, nx, nt).
    :param initial_idx: Baseline survey index.
    :param dx: Grid spacing in meters.
    :param deposit_per_year: Deposit volume matrix to update.
    :returns: Tuple ``(total_sed_vol_per_year, deposit_per_year, theseus_ratio)``.
    """
    z_min = np.nanmin(z[:, :, initial_idx:])
    total_sed_vol_per_year = np.nansum((z - z_min) * (dx**2), axis=(0, 1))

    nt = z.shape[2]
    if nt > 0:
        deposit_per_year[0, 0] = total_sed_vol_per_year[0]

    theseus_ratio = np.full((nt, nt), np.nan, dtype=float)
    for tt in range(nt):
        deposit_per_year[0, tt] = total_sed_vol_per_year[tt] - np.nansum(deposit_per_year[1:, tt])
        denom = deposit_per_year[tt, tt]
        if np.isfinite(denom) and denom != 0:
            theseus_ratio[tt, :] = deposit_per_year[tt, :] / denom

    return total_sed_vol_per_year, deposit_per_year, theseus_ratio


def _extract_line_xy(gdf) -> tuple[np.ndarray, np.ndarray] | None:
    """Extract coordinates from LineString or MultiLineString geometry.

    :param gdf: GeoDataFrame containing geometries.
    :returns: ``(x, y)`` arrays or ``None`` if no line geometry is found.
    """
    for geom in gdf.geometry:
        if geom is None:
            continue
        if geom.geom_type == "LineString":
            coords = np.asarray(geom.coords)
            return coords[:, 0], coords[:, 1]
        if geom.geom_type == "MultiLineString":
            all_coords = []
            for line in geom.geoms:
                all_coords.extend(list(line.coords))
            coords = np.asarray(all_coords)
            return coords[:, 0], coords[:, 1]
    return None


def _extract_xy_columns(gdf) -> tuple[np.ndarray, np.ndarray] | None:
    """Extract coordinates from likely x/y attribute columns.

    :param gdf: GeoDataFrame with attribute columns.
    :returns: ``(x, y)`` arrays or ``None`` if no columns match.
    """
    lower_to_col = {c.lower(): c for c in gdf.columns}

    x_keys = ["lon1", "x", "easting", "lon", "longitude"]
    y_keys = ["lat1", "y", "northing", "lat", "latitude"]

    x_col = next((lower_to_col[k] for k in x_keys if k in lower_to_col), None)
    y_col = next((lower_to_col[k] for k in y_keys if k in lower_to_col), None)

    if x_col is None or y_col is None:
        return None

    x = pd.to_numeric(gdf[x_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(gdf[y_col], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    return x[valid], y[valid]


def _resolve_source_crs(gdf, source_crs: str | None, prompt_if_missing: bool, shp_name: str) -> CRS:
    """Resolve a CRS from arguments, GeoDataFrame, or user prompt.

    :param gdf: GeoDataFrame with optional CRS metadata.
    :param source_crs: CRS override string (e.g., ``EPSG:4326``).
    :param prompt_if_missing: Whether to prompt the user if CRS is missing.
    :param shp_name: Shapefile name for error messaging.
    :returns: Resolved :class:`pyproj.CRS`.
    :raises ValueError: If CRS is missing and prompting is disabled.
    """
    if source_crs is not None:
        return CRS.from_user_input(source_crs)
    if gdf.crs is not None:
        return CRS.from_user_input(gdf.crs)
    if not prompt_if_missing:
        raise ValueError(f"Shapefile {shp_name} has no CRS and source_crs was not provided")

    response = input(
        f"Shapefile {shp_name} has no CRS metadata. Enter source CRS (e.g., EPSG:4326 or EPSG:32618): "
    ).strip()
    if not response:
        raise ValueError(f"No CRS provided for shapefile {shp_name}")
    return CRS.from_user_input(response)


def load_transects_from_shapefiles(
    shp_dir: str | Path,
    target_crs: str,
    source_crs: str | None = None,
    prompt_if_missing_crs: bool = True,
    prefer_xy_columns: bool = False,
) -> list[Transect]:
    """Load and reproject transects from shapefiles in a directory.

    :param shp_dir: Directory containing shapefiles.
    :param target_crs: Target CRS for reprojection.
    :param source_crs: Optional source CRS override.
    :param prompt_if_missing_crs: Prompt for CRS if shapefile lacks metadata.
    :param prefer_xy_columns: Prefer attribute columns (e.g., Lon/Lat) over geometry.
    :returns: List of :class:`Transect` objects.
    :raises FileNotFoundError: If no shapefiles are found.
    :raises ValueError: If transect coordinates cannot be derived.
    """
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise ImportError("geopandas is required for shapefile transect loading") from exc

    shp_paths = sorted(Path(shp_dir).glob("*.shp"))
    transects: list[Transect] = []

    if not shp_paths:
        raise FileNotFoundError(f"No shapefiles found in {shp_dir}")

    dst_crs = CRS.from_user_input(target_crs)

    for shp in shp_paths:
        gdf = gpd.read_file(shp)

        if prefer_xy_columns:
            xy = _extract_xy_columns(gdf)
            if xy is None:
                xy = _extract_line_xy(gdf)
        else:
            xy = _extract_line_xy(gdf)
            if xy is None:
                xy = _extract_xy_columns(gdf)
        if xy is None:
            raise ValueError(
                f"Could not derive transect coordinates from {shp.name}. "
                "Expected LineString geometry or columns like Lon1/Lat1 or x/y."
            )

        x_raw, y_raw = xy
        if len(x_raw) < 2:
            continue

        src_crs_obj = _resolve_source_crs(gdf, source_crs, prompt_if_missing_crs, shp.name)
        transformer = Transformer.from_crs(src_crs_obj, dst_crs, always_xy=True)
        x_m, y_m = transformer.transform(x_raw, y_raw)

        transects.append(Transect(name=shp.stem, x=np.asarray(x_m, dtype=float), y=np.asarray(y_m, dtype=float)))

    return transects


def build_manual_transects(transect_points: list[np.ndarray], names: list[str] | None = None) -> list[Transect]:
    """Build transects from manual Nx2 arrays of points.

    :param transect_points: List of arrays shaped (n, 2).
    :param names: Optional list of transect names.
    :returns: List of :class:`Transect` objects.
    :raises ValueError: If any transect array is not shape (n, 2).
    """
    out: list[Transect] = []
    for i, points in enumerate(transect_points):
        arr = np.asarray(points, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("Each transect in transect_points must be an Nx2 array")
        name = names[i] if names is not None and i < len(names) else f"T{i + 1}"
        out.append(Transect(name=name, x=arr[:, 0], y=arr[:, 1]))
    return out


def pick_interactive_transect(ax, name: str = "GUI") -> Transect:
    """Collect a transect polyline from interactive matplotlib clicks.

    :param ax: Matplotlib axes for plotting feedback.
    :param name: Transect name.
    :returns: Transect created from user clicks.
    :raises ValueError: If fewer than two points are selected.
    """
    print("Left-click to add points, press Enter to finish transect")
    pts = np.asarray(plt.ginput(n=-1, timeout=0), dtype=float)
    if pts.shape[0] < 2:
        raise ValueError("At least two points are required for a transect")
    ax.plot(pts[:, 0], pts[:, 1], "xk")
    return Transect(name=name, x=pts[:, 0], y=pts[:, 1])


def _resample_polyline(x: np.ndarray, y: np.ndarray, n_points: int = 1000) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resample a polyline to uniform distance spacing.

    :param x: X coordinates of the polyline.
    :param y: Y coordinates of the polyline.
    :param n_points: Number of output samples.
    :returns: Tuple ``(x_q, y_q, s_q)`` of resampled coordinates and distances.
    :raises ValueError: If the polyline length is zero.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    keep = np.ones_like(x, dtype=bool)
    keep[1:] = (np.diff(x) != 0) | (np.diff(y) != 0)
    x = x[keep]
    y = y[keep]

    dx = np.diff(x)
    dy = np.diff(y)
    seg_len = np.hypot(dx, dy)
    s = np.concatenate([[0.0], np.cumsum(seg_len)])

    if s[-1] <= 0:
        raise ValueError("Transect has zero length")

    s_q = np.linspace(0.0, s[-1], n_points)
    x_q = np.interp(s_q, s, x)
    y_q = np.interp(s_q, s, y)
    return x_q, y_q, s_q


def _interp_stack_along_transect(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    stack: np.ndarray,
    x_q: np.ndarray,
    y_q: np.ndarray,
) -> np.ndarray:
    """Interpolate a 3D stack along a transect polyline.

    :param x_axis: 1D x-axis for the grid.
    :param y_axis: 1D y-axis for the grid.
    :param stack: 3D array to sample, shape (ny, nx, nt).
    :param x_q: Transect x coordinates to sample.
    :param y_q: Transect y coordinates to sample.
    :returns: Array of shape (nt, n_points).
    """
    out = np.full((stack.shape[2], len(x_q)), np.nan, dtype=float)
    pts = np.column_stack([y_q, x_q])

    for k in range(stack.shape[2]):
        interp = RegularGridInterpolator(
            (y_axis, x_axis),
            stack[:, :, k],
            bounds_error=False,
            fill_value=np.nan,
        )
        out[k, :] = interp(pts)

    return out


def _interp_surface_along_transect(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    surface: np.ndarray,
    x_q: np.ndarray,
    y_q: np.ndarray,
) -> np.ndarray:
    """Interpolate a 2D surface along a transect polyline.

    :param x_axis: 1D x-axis for the grid.
    :param y_axis: 1D y-axis for the grid.
    :param surface: 2D surface to sample.
    :param x_q: Transect x coordinates to sample.
    :param y_q: Transect y coordinates to sample.
    :returns: Interpolated values along the transect.
    """
    interp = RegularGridInterpolator((y_axis, x_axis), surface, bounds_error=False, fill_value=np.nan)
    pts = np.column_stack([y_q, x_q])
    return interp(pts)


def extract_transect_cube(
    bathy: BathyCube,
    transect: Transect,
    n_points: int = 400,
    name: str | None = None,
) -> BathyCube:
    """Sample a 3D bathy cube along a transect into a 1D cube.

    :param bathy: Source bathymetry cube.
    :param transect: Transect polyline to sample.
    :param n_points: Number of samples along the transect.
    :param name: Optional name for the output cube.
    :returns: 1D BathyCube with shape (1, n_points, nt).
    """
    x_axis, y_axis = _derive_xy_axes(bathy.x, bathy.y)
    x_q, y_q, s_q = _resample_polyline(transect.x, transect.y, n_points=n_points)
    z_transect = _interp_stack_along_transect(x_axis, y_axis, bathy.z, x_q, y_q)

    x_1d = np.asarray(s_q, dtype=float)
    x_2d = x_1d[None, :]
    y_2d = np.zeros_like(x_2d)
    z_3d = z_transect.T[None, :, :]

    location = name or f"{bathy.location}_{transect.name}_slice"
    return BathyCube(location=location, t=bathy.t, x=x_2d, y=y_2d, z=z_3d)


def _default_cmap(n: int) -> Colormap:
    """Return a discrete matplotlib colormap with at least two bins.

    :param n: Number of bins.
    :returns: Matplotlib colormap.
    """
    return plt.get_cmap("viridis", max(n, 2))


def plot_transect_location_map(
    bathy: BathyCube,
    transects: list[Transect],
    out_path: str | Path,
    config: StratigraphyConfig,
) -> None:
    """Plot bathymetry map with transect polylines and markers.

    :param bathy: Bathymetry cube.
    :param transects: Transects to render.
    :param out_path: Output image path.
    :param config: Plotting configuration.
    """
    fig, ax = plt.subplots(figsize=(12, 8), dpi=150)
    font = _get_plot_font()
    x_km = bathy.x / 1000.0
    y_km = bathy.y / 1000.0

    cont = ax.contourf(x_km, y_km, bathy.z[:, :, -1], levels=np.arange(-30, 5.2, 0.2), cmap="viridis")
    ax.contour(x_km, y_km, bathy.z[:, :, -1], levels=[MLW], colors=[(0.5, 0.5, 0.5)], linewidths=1.0)
    ax.contour(x_km, y_km, bathy.z[:, :, -1], levels=[-6], colors="k", linestyles=":", linewidths=0.5)

    for tr in transects:
        xq, yq, dq = _resample_polyline(tr.x, tr.y, n_points=1000)
        xq_km = xq / 1000.0
        yq_km = yq / 1000.0
        dq_km = dq / 1000.0

        ax.plot(xq_km, yq_km, "-k", lw=1.0)
        ax.scatter([xq_km[0]], [yq_km[0]], s=config.end_dot_size, c="w", edgecolors="k", zorder=3)
        ax.scatter([xq_km[-1]], [yq_km[-1]], s=config.end_dot_size, c="w", edgecolors="k", zorder=3)

        prev_bin = -1
        for i in range(len(dq_km)):
            cur_bin = int(np.floor(dq_km[i] / config.dot_spacing_km))
            if cur_bin > prev_bin:
                ax.scatter([xq_km[i]], [yq_km[i]], s=config.mid_dot_size, c="w", edgecolors="k", zorder=3)
                prev_bin = cur_bin

        ax.text(
            xq_km[0] + 0.03,
            yq_km[0] - 0.03,
            tr.name,
            fontsize=10,
            fontweight="bold",
            fontproperties=font,
        )
        ax.text(
            xq_km[-1] - 0.08,
            yq_km[-1] + 0.03,
            f"{tr.name}'",
            fontsize=10,
            fontweight="bold",
            fontproperties=font,
        )

    cb = fig.colorbar(cont, ax=ax)
    cb.set_label("Elevation [m NAVD88]", fontproperties=font)

    survey_year = pd.Timestamp(datenum_to_datetime64(np.asarray([bathy.t[-1]]))[0]).year
    ax.set_title(f"Transect Locations ({survey_year} Bathymetry)", fontproperties=font)
    ax.set_xlabel("Easting [km]", fontproperties=font)
    ax.set_ylabel("Northing [km]", fontproperties=font)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
    _apply_axes_font(ax, font)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_cross_sections(
    bathy: BathyCube,
    result: StratigraphyResult,
    transects: list[Transect],
    out_dir: str | Path,
) -> None:
    """Plot stratigraphic cross-sections for each transect.

    :param bathy: Bathymetry cube.
    :param result: Stratigraphy results.
    :param transects: Transects to render.
    :param out_dir: Output directory for images.
    """
    x_axis, y_axis = _derive_xy_axes(bathy.x, bathy.y)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmap = _default_cmap(result.deposit_thk_full.shape[2])
    font = _get_plot_font()
    colors = cmap(np.linspace(0, 1, max(result.deposit_thk_full.shape[2], 2)))

    for tr in transects:
        xq, yq, dq_m = _resample_polyline(tr.x, tr.y, n_points=1000)
        dq_km = dq_m / 1000.0

        z_min = _interp_surface_along_transect(x_axis, y_axis, result.min_surf, xq, yq)
        z_max = _interp_surface_along_transect(x_axis, y_axis, result.max_surf, xq, yq)
        z_now = _interp_surface_along_transect(x_axis, y_axis, bathy.z[:, :, -1], xq, yq)
        dep_xs = _interp_stack_along_transect(x_axis, y_axis, result.deposit_thk_full, xq, yq)

        fig, ax = plt.subplots(figsize=(14, 4), dpi=150)

        baseline = -35.0 * np.ones_like(dq_km)
        running = baseline.copy()
        for k in range(dep_xs.shape[0]):
            layer = np.nan_to_num(dep_xs[k, :], nan=0.0)
            next_running = running + layer
            # Stack deposit thickness layers from a fixed baseline.
            ax.fill_between(dq_km, running, next_running, color=colors[k], linewidth=0.0)
            running = next_running

        ax.plot(dq_km, z_min, "-k", lw=1.0)
        ax.plot(dq_km, z_max, ":k", lw=1.0)
        ax.plot(dq_km, z_now, "-k", lw=1.5)

        ax.axhline(MLW, linestyle="--", color=(0.1, 0.2, 0.5), linewidth=0.6)
        ax.axhline(MHW, linestyle="--", color=(0.1, 0.2, 0.5), linewidth=0.6)
        ax.text(dq_km[0] + 0.02, MHW + 0.3, "MHW", color=(0.1, 0.2, 0.5), fontsize=9, fontproperties=font)
        ax.text(dq_km[0] + 0.02, MLW + 0.3, "MLW", color=(0.1, 0.2, 0.5), fontsize=9, fontproperties=font)

        ax.set_xlabel("Distance [km]", fontproperties=font)
        ax.set_ylabel("Elevation [m NAVD88]", fontproperties=font)
        ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
        ax.set_xlim(0.0, np.nanmax(dq_km))
        _apply_axes_font(ax, font)

        y_min = min(np.nanmin(z_min) - 1.0, np.nanmin(baseline) - 1.0)
        y_max = np.nanmax([np.nanmax(z_now) + 1.0, 3.0])
        if np.isfinite(y_min) and np.isfinite(y_max) and y_max > y_min:
            ax.set_ylim(y_min, y_max)

        fig.tight_layout()
        fig.savefig(out_dir / f"Cross-section {tr.name}-{tr.name}'.png")
        plt.close(fig)


def plot_theseus_ratio(result: StratigraphyResult, out_path: str | Path) -> None:
    """Plot Theseus ratio curves in linear and log-log space.

    :param result: Stratigraphy results.
    :param out_path: Output image path.
    """
    t_dt = datenum_to_datetime64(result.t)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=150)
    font = _get_plot_font()

    for tt in range(result.theseus_ratio.shape[0]):
        y = result.theseus_ratio[tt, tt:]
        x = t_dt[tt:]
        valid = np.isfinite(y)
        if np.any(valid):
            ax1.plot(x[valid], y[valid], linewidth=1.0)

    ax1.set_xlabel("Time", fontproperties=font)
    ax1.set_ylabel("Theta (Fraction Preserved) [-]", fontproperties=font)
    ax1.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)

    for tt in range(result.theseus_ratio.shape[0]):
        y = result.theseus_ratio[tt, tt:]
        x = (result.t[tt:] - result.t[tt]) / 10.0
        valid = np.isfinite(y) & np.isfinite(x) & (x > 0) & (y > 0)
        if np.any(valid):
            ax2.plot(x[valid], y[valid], linewidth=1.0)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Time", fontproperties=font)
    ax2.set_ylabel("Theta (Fraction Preserved) [-]", fontproperties=font)
    ax2.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
    _apply_axes_font(ax1, font)
    _apply_axes_font(ax2, font)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def save_metrics(result: StratigraphyResult, out_dir: str | Path) -> None:
    """Write stratigraphy metrics and grids to disk.

    :param result: Stratigraphy results.
    :param out_dir: Output directory.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(result.deposit_per_year).to_csv(out_dir / "deposit_per_year.csv", index=False)
    pd.DataFrame(result.theseus_ratio).to_csv(out_dir / "theseus_ratio.csv", index=False)
    pd.DataFrame({"total_sed_vol_per_year": result.total_sed_vol_per_year}).to_csv(
        out_dir / "total_sed_vol_per_year.csv", index=False
    )
    pd.DataFrame({"survey_datenum": result.t}).to_csv(out_dir / "survey_datenum.csv", index=False)
    pd.DataFrame({"monthly_datenum": result.t_full}).to_csv(out_dir / "monthly_datenum.csv", index=False)

    np.savez_compressed(out_dir / "deposit_thk_full.npz", deposit_thk_full=result.deposit_thk_full)
    np.savez_compressed(out_dir / "min_surf.npz", min_surf=result.min_surf)
    np.savez_compressed(out_dir / "max_surf.npz", max_surf=result.max_surf)


def run_stratigraphy_workflow(
    output_root: str | Path,
    transect_mode: str,
    bathy_nc_path: str | Path | None = None,
    bathy_mat_path: str | Path | None = None,
    shp_dir: str | Path | None = None,
    manual_transects: list[np.ndarray] | None = None,
    source_crs: str | None = None,
    config: StratigraphyConfig | None = None,
) -> StratigraphyResult:
    """Run the full stratigraphy workflow including plotting and exports.

    :param bathy_nc_path: Path to the bathy .nc file (preferred).
    :param bathy_mat_path: Optional path to the bathy MAT file (legacy).
    :param output_root: Root output directory.
    :param transect_mode: ``shapefile``, ``manual``, or ``gui``.
    :param shp_dir: Shapefile directory when using ``shapefile`` mode.
    :param manual_transects: Manual transect point arrays.
    :param source_crs: Optional source CRS override for shapefiles.
    :param config: Optional stratigraphy configuration.
    :returns: Stratigraphy results.
    :raises ValueError: If required transect inputs are missing.
    """
    config = config or StratigraphyConfig()
    bathy_path = bathy_nc_path or bathy_mat_path
    if bathy_path is None:
        raise ValueError("bathy_nc_path is required when bathy_mat_path is not provided")
    bathy = load_bathy_cube(bathy_path)
    result = compute_stratigraphy(bathy, config)

    output_root = Path(output_root)
    strat_dir = output_root / "plots" / "python" / "stratigraphy"
    metrics_dir = output_root / "metrics" / "stratigraphy"

    transects: list[Transect] = []
    if transect_mode == "shapefile":
        if shp_dir is None:
            raise ValueError("shp_dir is required when transect_mode='shapefile'")
        transects = load_transects_from_shapefiles(
            shp_dir=shp_dir,
            target_crs=config.target_crs,
            source_crs=source_crs,
            prompt_if_missing_crs=True,
        )
    elif transect_mode == "manual":
        if not manual_transects:
            raise ValueError("manual_transects is required when transect_mode='manual'")
        transects = build_manual_transects(manual_transects)
    elif transect_mode == "gui":
        fig, ax = plt.subplots(figsize=(12, 8), dpi=150)
        ax.contourf(bathy.x / 1000.0, bathy.y / 1000.0, bathy.z[:, :, -1], levels=np.arange(-30, 5.2, 0.2))
        ax.set_title("Click transect points, then press Enter")
        tr = pick_interactive_transect(ax=ax, name="GUI")
        plt.close(fig)
        transects = [tr]
    else:
        raise ValueError("transect_mode must be one of: shapefile, manual, gui")

    if transects:
        plot_transect_location_map(
            bathy,
            transects,
            strat_dir / "Stratigraphic Transect Map (Python).png",
            config,
        )
        plot_cross_sections(bathy, result, transects, strat_dir)

    plot_theseus_ratio(result, strat_dir / "TheseusRatio.png")
    save_metrics(result, metrics_dir)

    return result


def _build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser.

    :returns: Configured argument parser.
    """
    parser = argparse.ArgumentParser(description="Compute stratigraphy from bathymetry input")
    parser.add_argument("--bathy-nc", default=None, help="Path to bathy netCDF file")
    parser.add_argument("--bathy-mat", default=None, help="Path to bathy MAT file with known structure")
    parser.add_argument("--output-root", default=".", help="Root output folder")
    parser.add_argument(
        "--transect-mode",
        required=True,
        choices=["shapefile", "manual", "gui"],
        help="Transect source mode",
    )
    parser.add_argument("--shp-dir", default=None, help="Folder of shapefiles when --transect-mode shapefile")
    parser.add_argument("--source-crs", default=None, help="Optional source CRS override, e.g. EPSG:4326")
    parser.add_argument("--target-crs", default="EPSG:32618", help="Target projected CRS")
    parser.add_argument("--initial-index", type=int, default=0, help="Initial survey index used as baseline")
    parser.add_argument("--dx", type=float, default=20.0, help="Grid cell size [m]")
    return parser


def main() -> None:
    """CLI entry point."""
    parser = _build_parser()
    args = parser.parse_args()

    if args.bathy_nc is None and args.bathy_mat is None:
        parser.error("--bathy-nc or --bathy-mat is required")

    config = StratigraphyConfig(
        initial_index=args.initial_index,
        dx=args.dx,
        target_crs=args.target_crs,
    )

    run_stratigraphy_workflow(
        bathy_nc_path=args.bathy_nc,
        bathy_mat_path=args.bathy_mat,
        output_root=args.output_root,
        transect_mode=args.transect_mode,
        shp_dir=args.shp_dir,
        source_crs=args.source_crs,
        config=config,
    )


if __name__ == "__main__":
    main()
