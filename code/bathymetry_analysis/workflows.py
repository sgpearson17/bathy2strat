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
from utils.colormaps import bathymetry_colormap
from bathy_formatter import get_named_colormap
from .bathy import load_clrmap_file
from .plot_style import apply_axes_font, apply_global_style, get_plot_font
from .stratigraphy import (
    BathyCube,
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


def _resolve_date_index(dates: np.ndarray, requested: str | np.datetime64) -> tuple[int, np.datetime64]:
    """Resolve a requested date to the nearest available survey date at or after it."""
    if dates.size == 0:
        raise ValueError("No survey dates are available in the bathymetry cube")
    requested_dt = np.datetime64(requested)
    idx = int(np.searchsorted(dates, requested_dt, side="left"))
    if idx >= dates.size:
        idx = dates.size - 1
    return idx, dates[idx]


def _deposit_layer_thickness(result, layer_idx: int) -> np.ndarray:
    """Return preserved thickness for a specific deposit layer index."""
    layer_idx = int(layer_idx)
    if layer_idx < 0 or layer_idx >= result.deposit_elev.shape[2]:
        raise IndexError(f"Deposit layer index out of bounds: {layer_idx}")
    if layer_idx == 0:
        thickness = np.maximum(0.0, result.deposit_elev[:, :, 0] - result.min_surf)
    else:
        thickness = np.maximum(0.0, result.deposit_elev[:, :, layer_idx] - result.deposit_elev[:, :, layer_idx - 1])
    valid = np.isfinite(result.deposit_elev[:, :, layer_idx])
    return np.where(valid, thickness, np.nan)


def _deposit_layer_thickness_stack(result) -> np.ndarray:
    """Return preserved thickness stack for all deposit layers."""
    dep = np.asarray(result.deposit_elev, dtype=float)
    ny, nx, nt = dep.shape
    thk = np.full((ny, nx, nt), np.nan, dtype=float)
    if nt == 0:
        return thk
    thk[:, :, 0] = np.maximum(0.0, dep[:, :, 0] - np.asarray(result.min_surf, dtype=float))
    for k in range(1, nt):
        thk[:, :, k] = np.maximum(0.0, dep[:, :, k] - dep[:, :, k - 1])
    valid = np.isfinite(dep)
    return np.where(valid, thk, np.nan)


def _auto_contour_levels(surface: np.ndarray, interval_m: float) -> np.ndarray:
    """Build contour levels from a 2D surface using a target interval."""
    valid = np.asarray(surface, dtype=float)
    finite = np.isfinite(valid)
    if not np.any(finite):
        return np.array([-10.0, -5.0, 0.0, 5.0], dtype=float)

    step = float(interval_m) if np.isfinite(interval_m) and interval_m > 0 else 1.0
    zmin = float(np.nanmin(valid[finite]))
    zmax = float(np.nanmax(valid[finite]))
    lo = np.floor(zmin / step) * step
    hi = np.ceil(zmax / step) * step
    levels = np.arange(lo, hi + 0.5 * step, step, dtype=float)
    if levels.size > 40:
        stride = int(np.ceil(levels.size / 40.0))
        levels = levels[::stride]
    if levels.size < 2:
        levels = np.array([lo, lo + step], dtype=float)
    return levels


def _nice_age_levels(*maps: np.ndarray, target_bins: int = 8) -> np.ndarray:
    """Build rounded age color levels shared across one or more maps."""
    max_age = 0.0
    for arr in maps:
        if arr is None:
            continue
        if np.any(np.isfinite(arr)):
            max_age = max(max_age, float(np.nanmax(arr)))

    if not np.isfinite(max_age) or max_age <= 0.0:
        return np.array([0.0, 1.0], dtype=float)

    raw_step = max_age / max(float(target_bins), 1.0)
    magnitude = 10.0 ** np.floor(np.log10(raw_step))
    scaled = raw_step / magnitude
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if scaled <= m:
            step = m * magnitude
            break
    else:
        step = 10.0 * magnitude

    upper = np.ceil(max_age / step) * step
    levels = np.arange(0.0, upper + 0.5 * step, step, dtype=float)
    if levels.size < 2:
        levels = np.array([0.0, step], dtype=float)
    return levels


def _year_fraction_from_datenum(t_vals: np.ndarray) -> np.ndarray:
    """Convert MATLAB datenums to fractional calendar years."""
    t_arr = np.asarray(t_vals, dtype=float).reshape(-1)
    if t_arr.size == 0:
        return np.array([], dtype=float)
    t_dt = pd.DatetimeIndex(pd.to_datetime(datenum_to_datetime64(t_arr), errors="coerce"))
    return t_dt.year.to_numpy(dtype=float) + (t_dt.dayofyear.to_numpy(dtype=float) - 1.0) / 365.25


def _youngest_preserved_layer_index(thickness_stack: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """Return index of youngest preserved layer per cell, or -1 where none exist."""
    valid = np.isfinite(thickness_stack) & (thickness_stack > eps)
    has = np.any(valid, axis=2)
    rev_idx = np.argmax(valid[:, :, ::-1], axis=2)
    nt = thickness_stack.shape[2]
    out = (nt - 1 - rev_idx).astype(int)
    out[~has] = -1
    return out


def _top_layer_average_age_map(
    thickness_stack: np.ndarray,
    layer_age_years: np.ndarray,
    top_layer_thickness_m: float,
    eps: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-cell average age of the top-thickness layer and used thickness."""
    ny, nx, nt = thickness_stack.shape
    age_map = np.full((ny, nx), np.nan, dtype=float)
    used_thickness = np.zeros((ny, nx), dtype=float)
    target = float(top_layer_thickness_m)
    if not np.isfinite(target) or target <= 0.0:
        return age_map, used_thickness

    for iy in range(ny):
        for ix in range(nx):
            remaining = target
            weighted = 0.0
            used = 0.0
            for k in range(nt - 1, -1, -1):
                thk = thickness_stack[iy, ix, k]
                if not np.isfinite(thk) or thk <= eps:
                    continue
                take = min(float(thk), remaining)
                if take <= 0.0:
                    continue
                weighted += take * float(layer_age_years[k])
                used += take
                remaining -= take
                if remaining <= eps:
                    break
            if used > eps:
                age_map[iy, ix] = weighted / used
                used_thickness[iy, ix] = used
    return age_map, used_thickness


def plot_highlight_deposit_thickness_maps(
    bathy_source: str | Path | BathyCube,
    output_root: str | Path,
    highlight_dates: str | np.datetime64 | list[str | np.datetime64],
    initial_index: int = 0,
    dx: float = 20.0,
    contour_interval_m: float = 1.0,
    contour_levels: np.ndarray | list[float] | None = None,
    thickness_cmap: str = "magma_r",
) -> dict[str, object]:
    """Plot local preserved deposit thickness maps for highlighted dates.

    Args:
        bathy_source: Bathymetry source path or already-loaded BathyCube.
        output_root: Output root directory.
        highlight_dates: One or more dates to map, one map per date.
        initial_index: Baseline stratigraphy index.
        dx: Grid spacing in meters.
        contour_interval_m: Contour interval for latest bathymetry overlay.
        contour_levels: Optional explicit contour levels (meters).
        thickness_cmap: Colormap used for local preserved deposit thickness.

    Returns:
        Dict with output directory and per-date map metadata.
    """
    dates_to_plot = _normalize_highlight_dates(highlight_dates)
    if not dates_to_plot:
        raise ValueError("highlight_dates must include at least one date")

    bathy = _coerce_bathy_cube(bathy_source)
    config = StratigraphyConfig(initial_index=initial_index, dx=dx)
    result = compute_stratigraphy(bathy, config)
    survey_dates = datenum_to_datetime64(result.t)
    latest_surface = np.asarray(bathy.z[:, :, -1], dtype=float)
    contour_tag = str(survey_dates[-1])[:10] if survey_dates.size > 0 else "latest"

    x_km = np.asarray(bathy.x, dtype=float) / 1000.0
    y_km = np.asarray(bathy.y, dtype=float) / 1000.0

    out_dir = Path(output_root) / "plots" / "python" / "stratigraphy" / "spatial_deposit_thickness"
    out_dir.mkdir(parents=True, exist_ok=True)

    if contour_levels is None:
        levels = _auto_contour_levels(latest_surface, contour_interval_m)
    else:
        levels = np.asarray(contour_levels, dtype=float)
        levels = levels[np.isfinite(levels)]
        levels = np.unique(levels)
        if levels.size < 2:
            levels = _auto_contour_levels(latest_surface, contour_interval_m)

    outputs: list[dict[str, object]] = []
    for requested in dates_to_plot:
        idx, resolved = _resolve_date_index(survey_dates, requested)
        thk = _deposit_layer_thickness(result, idx)

        vmax = float(np.nanmax(thk)) if np.any(np.isfinite(thk)) else 0.0
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0
        thk_levels = np.linspace(0.0, vmax, 21)

        apply_global_style()
        font = get_plot_font()
        fig, ax = plt.subplots(figsize=(9.5, 7.5))
        cf = ax.contourf(
            x_km,
            y_km,
            thk,
            levels=thk_levels,
            cmap=thickness_cmap,
            extend="max",
        )
        cs = ax.contour(
            x_km,
            y_km,
            latest_surface,
            levels=levels,
            colors="k",
            linewidths=0.55,
            alpha=0.45,
        )
        if len(cs.levels) > 0:
            ax.clabel(cs, cs.levels[::2], inline=True, fontsize=7, fmt="%.0f")

        cb = fig.colorbar(cf, ax=ax)
        cb.set_label("Preserved deposit thickness [m]", fontproperties=font)
        for tick in cb.ax.get_yticklabels():
            tick.set_fontproperties(font)

        tag = str(resolved)[:10]
        ax.set_title(f"Local preserved deposit thickness ({tag}) with {contour_tag} depth contours")
        ax.set_xlabel("Easting [km]")
        ax.set_ylabel("Northing [km]")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, color="0.5", alpha=0.35)
        apply_axes_font(ax, font)

        file_path = out_dir / f"deposit_thickness_map_{tag}.png"
        fig.tight_layout()
        fig.savefig(file_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        outputs.append(
            {
                "requested_date": str(np.datetime64(requested))[:10],
                "resolved_survey_date": tag,
                "layer_index": int(idx),
                "map_path": file_path,
            }
        )

    return {
        "output_dir": out_dir,
        "maps": outputs,
    }


def plot_domain_surface_age_maps(
    bathy_source: str | Path | BathyCube,
    output_root: str | Path,
    initial_index: int = 0,
    dx: float = 20.0,
    top_layer_thickness_m: float = 0.5,
    contour_interval_m: float = 1.0,
    contour_levels: np.ndarray | list[float] | None = None,
    surface_age_cmap: str = "viridis_r",
    top_layer_age_cmap: str | None = None,
) -> dict[str, object]:
    """Plot full-domain surface age and average age of the top layer.

    Args:
        bathy_source: Bathymetry source path or already-loaded BathyCube.
        output_root: Output root directory.
        initial_index: Baseline stratigraphy index.
        dx: Grid spacing in meters.
        top_layer_thickness_m: Thickness of the top layer used for average-age calculation.
        contour_interval_m: Contour interval for latest bathymetry overlay.
        contour_levels: Optional explicit contour levels (meters).
        surface_age_cmap: Colormap for surface age map.
        top_layer_age_cmap: Colormap for top-layer average age map (defaults to surface_age_cmap).

    Returns:
        Dict with figure paths and summary statistics.
    """
    bathy = _coerce_bathy_cube(bathy_source)
    config = StratigraphyConfig(initial_index=initial_index, dx=dx)
    result = compute_stratigraphy(bathy, config)

    x_km = np.asarray(bathy.x, dtype=float) / 1000.0
    y_km = np.asarray(bathy.y, dtype=float) / 1000.0
    latest_surface = np.asarray(bathy.z[:, :, -1], dtype=float)

    if contour_levels is None:
        levels = _auto_contour_levels(latest_surface, contour_interval_m)
    else:
        levels = np.asarray(contour_levels, dtype=float)
        levels = levels[np.isfinite(levels)]
        levels = np.unique(levels)
        if levels.size < 2:
            levels = _auto_contour_levels(latest_surface, contour_interval_m)

    survey_dates = datenum_to_datetime64(result.t)
    contour_tag = str(survey_dates[-1])[:10] if survey_dates.size > 0 else "latest"
    layer_year = _year_fraction_from_datenum(result.t)
    if layer_year.size == 0:
        raise ValueError("No survey dates available to compute age metrics")
    latest_year = float(layer_year[-1])
    layer_age_years = np.maximum(0.0, latest_year - layer_year)

    thickness_stack = _deposit_layer_thickness_stack(result)
    youngest_idx = _youngest_preserved_layer_index(thickness_stack)
    surface_age_years = np.full(youngest_idx.shape, np.nan, dtype=float)
    has_surface = youngest_idx >= 0
    if np.any(has_surface):
        surface_age_years[has_surface] = layer_age_years[youngest_idx[has_surface]]
    surface_age_years = np.where(np.isfinite(latest_surface), surface_age_years, np.nan)

    top_age_map, top_used_thickness = _top_layer_average_age_map(
        thickness_stack,
        layer_age_years,
        top_layer_thickness_m=top_layer_thickness_m,
    )
    top_age_map = np.where(np.isfinite(latest_surface), top_age_map, np.nan)
    domain_avg_top_age_years = float(np.nanmean(top_age_map)) if np.any(np.isfinite(top_age_map)) else np.nan

    if top_layer_age_cmap is None:
        top_layer_age_cmap = surface_age_cmap
    age_levels = _nice_age_levels(surface_age_years, top_age_map)

    apply_global_style()
    font = get_plot_font()
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 6.0), constrained_layout=True)

    # Surface age map.
    ax = axes[0]
    cf0 = ax.contourf(x_km, y_km, surface_age_years, levels=age_levels, cmap=surface_age_cmap, extend="max")
    cs0 = ax.contour(x_km, y_km, latest_surface, levels=levels, colors="k", linewidths=0.55, alpha=0.45)
    if len(cs0.levels) > 0:
        ax.clabel(cs0, cs0.levels[::2], inline=True, fontsize=7, fmt="%.0f")
    cb0 = fig.colorbar(cf0, ax=ax)
    cb0.set_label("Surface age [years]", fontproperties=font)
    for tick in cb0.ax.get_yticklabels():
        tick.set_fontproperties(font)
    ax.set_title(f"Surface age with {contour_tag} depth contours")
    ax.set_xlabel("Easting [km]")
    ax.set_ylabel("Northing [km]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.5", alpha=0.35)
    apply_axes_font(ax, font)

    # Top-layer average age map.
    ax = axes[1]
    cf1 = ax.contourf(x_km, y_km, top_age_map, levels=age_levels, cmap=top_layer_age_cmap, extend="max")
    cs1 = ax.contour(x_km, y_km, latest_surface, levels=levels, colors="k", linewidths=0.55, alpha=0.45)
    if len(cs1.levels) > 0:
        ax.clabel(cs1, cs1.levels[::2], inline=True, fontsize=7, fmt="%.0f")
    cb1 = fig.colorbar(cf1, ax=ax)
    cb1.set_label(f"Average age of top {top_layer_thickness_m:.2f} m [years]", fontproperties=font)
    for tick in cb1.ax.get_yticklabels():
        tick.set_fontproperties(font)
    ax.set_title(f"Top-layer average age with {contour_tag} depth contours")
    ax.set_xlabel("Easting [km]")
    ax.set_ylabel("Northing [km]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.5", alpha=0.35)
    apply_axes_font(ax, font)

    out_dir = Path(output_root) / "plots" / "python" / "stratigraphy" / "surface_age"
    out_dir.mkdir(parents=True, exist_ok=True)
    map_path = out_dir / "surface_age_and_top_layer_age.png"
    fig.savefig(map_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Coverage map to show where requested top-layer thickness was available.
    fig_cov, ax_cov = plt.subplots(1, 1, figsize=(7.1, 6.0), constrained_layout=True)
    vmax_cov = max(float(top_layer_thickness_m), 1e-6)
    lv_cov = np.linspace(0.0, vmax_cov, 21)
    cfc = ax_cov.contourf(x_km, y_km, top_used_thickness, levels=lv_cov, cmap="Greys", extend="max")
    csc = ax_cov.contour(x_km, y_km, latest_surface, levels=levels, colors="k", linewidths=0.55, alpha=0.45)
    if len(csc.levels) > 0:
        ax_cov.clabel(csc, csc.levels[::2], inline=True, fontsize=7, fmt="%.0f")
    cbc = fig_cov.colorbar(cfc, ax=ax_cov)
    cbc.set_label("Used top-layer thickness [m]", fontproperties=font)
    for tick in cbc.ax.get_yticklabels():
        tick.set_fontproperties(font)
    ax_cov.set_title("Top-layer thickness used in age average")
    ax_cov.set_xlabel("Easting [km]")
    ax_cov.set_ylabel("Northing [km]")
    ax_cov.set_aspect("equal", adjustable="box")
    ax_cov.grid(True, color="0.5", alpha=0.35)
    apply_axes_font(ax_cov, font)
    coverage_path = out_dir / "top_layer_age_thickness_used.png"
    fig_cov.savefig(coverage_path, dpi=300, bbox_inches="tight")
    plt.close(fig_cov)

    return {
        "age_map_path": map_path,
        "top_layer_coverage_path": coverage_path,
        "output_dir": out_dir,
        "top_layer_thickness_m": float(top_layer_thickness_m),
        "domain_avg_top_age_years": domain_avg_top_age_years,
        "surface_age_years": surface_age_years,
        "top_layer_avg_age_years_map": top_age_map,
        "top_layer_used_thickness_map": top_used_thickness,
    }


def _domain_time_axis_from_datenum(t_vals: np.ndarray, start_year_at_zero: bool) -> tuple[np.ndarray, str]:
    """Return a plottable time axis and axis label from MATLAB datenums."""
    t_arr = np.asarray(t_vals, dtype=float).reshape(-1)
    if t_arr.size == 0:
        return np.array([], dtype=float), "Time step"

    t_dt = pd.DatetimeIndex(pd.to_datetime(datenum_to_datetime64(t_arr), errors="coerce"))
    if start_year_at_zero:
        delta_days = ((t_dt - t_dt[0]) / np.timedelta64(1, "D")).to_numpy(dtype=float)
        return delta_days, "Time since first survey [days]"

    year_frac = t_dt.year.to_numpy(dtype=float) + (t_dt.dayofyear.to_numpy(dtype=float) - 1.0) / 365.25
    return year_frac, "Year"


def plot_domain_preservation_metrics(
    bathy_source: str | Path | BathyCube,
    output_root: str | Path,
    initial_index: int = 0,
    dx: float = 20.0,
    start_year_at_zero: bool = True,
    highlight_date: str | np.datetime64 | None = None,
) -> dict[str, object]:
    """Plot full-domain preserved-volume metrics (absolute, normalized, and Theseus ratio).

    Args:
        bathy_source: Bathymetry source path or already-loaded BathyCube.
        output_root: Output root directory.
        initial_index: Baseline stratigraphy index.
        dx: Grid spacing in meters.
        start_year_at_zero: If True, x-axis is days since first survey.
        highlight_date: Optional date to highlight in red.

    Returns:
        Dict containing output figure path and metadata.
    """
    bathy = _coerce_bathy_cube(bathy_source)
    config = StratigraphyConfig(initial_index=initial_index, dx=dx)
    result = compute_stratigraphy(bathy, config)

    nt = int(result.deposit_per_year.shape[0])
    time_axis, time_label = _domain_time_axis_from_datenum(result.t, start_year_at_zero)
    if time_axis.size == 0:
        time_axis = np.arange(nt, dtype=float)
        time_label = "Time step"

    survey_dates = datenum_to_datetime64(result.t)
    highlight_idx = None
    highlight_tag = None
    if highlight_date is not None and survey_dates.size > 0:
        idx, resolved = _resolve_date_index(survey_dates, highlight_date)
        highlight_idx = int(idx)
        highlight_tag = str(resolved)[:10]

    apply_global_style()
    font = get_plot_font()
    colors = plt.cm.viridis(np.linspace(0.15, 0.95, max(nt, 2)))
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), constrained_layout=True)

    ax_abs, ax_norm, ax_ratio = axes
    for k in range(1, nt):
        y_abs = np.asarray(result.deposit_per_year[k, :], dtype=float)
        y_abs[:k] = np.nan
        ax_abs.plot(time_axis, y_abs, linewidth=1.2, color=colors[k])

        denom = float(result.deposit_per_year[k, k]) if np.isfinite(result.deposit_per_year[k, k]) else np.nan
        y_norm = np.full_like(y_abs, np.nan, dtype=float)
        if np.isfinite(denom) and denom != 0.0:
            y_norm[k:] = y_abs[k:] / denom
        else:
            y_norm[k:] = 0.0
        ax_norm.plot(time_axis, y_norm, linewidth=1.2, color=colors[k])

    if highlight_idx is not None and highlight_idx >= 1:
        y_abs = np.asarray(result.deposit_per_year[highlight_idx, :], dtype=float)
        y_abs[:highlight_idx] = np.nan
        ax_abs.plot(time_axis, y_abs, linewidth=2.4, color="red", zorder=3)

        denom = float(result.deposit_per_year[highlight_idx, highlight_idx]) if np.isfinite(result.deposit_per_year[highlight_idx, highlight_idx]) else np.nan
        y_norm = np.full_like(y_abs, np.nan, dtype=float)
        if np.isfinite(denom) and denom != 0.0:
            y_norm[highlight_idx:] = y_abs[highlight_idx:] / denom
        else:
            y_norm[highlight_idx:] = 0.0
        ax_norm.plot(time_axis, y_norm, linewidth=2.4, color="red", zorder=3)

    denom0 = float(result.deposit_per_year[0, 0]) if np.isfinite(result.deposit_per_year[0, 0]) else np.nan
    ratio_t0 = np.full(nt, np.nan, dtype=float)
    if np.isfinite(denom0) and denom0 != 0.0:
        ratio_t0 = np.clip(np.asarray(result.deposit_per_year[0, :], dtype=float) / denom0, 0.0, 1.0)
    ax_ratio.plot(time_axis, ratio_t0, color="k", linewidth=1.4)
    if highlight_idx is not None and highlight_idx < len(time_axis):
        ax_ratio.axvline(time_axis[highlight_idx], color="red", linewidth=1.6, zorder=3)

    if highlight_tag:
        fig.suptitle(f"Domain preservation metrics (highlight: {highlight_tag})", fontproperties=font)
    else:
        fig.suptitle("Domain preservation metrics", fontproperties=font)

    ax_abs.set_title("(d) Volume preserved (absolute)", fontproperties=font)
    ax_abs.set_xlabel(time_label, fontproperties=font)
    ax_abs.set_ylabel("Volume [m^3]", fontproperties=font)
    ax_abs.grid(True, alpha=0.3)

    ax_norm.set_title("(e) Volume preserved (normalized)", fontproperties=font)
    ax_norm.set_xlabel(time_label, fontproperties=font)
    ax_norm.set_ylabel("Fraction of initial", fontproperties=font)
    ax_norm.grid(True, alpha=0.3)

    ax_ratio.set_title("(f) Theseus ratio (t0 preserved)", fontproperties=font)
    ax_ratio.set_xlabel(time_label, fontproperties=font)
    ax_ratio.set_ylabel("Fraction of initial", fontproperties=font)
    if np.any(np.isfinite(ratio_t0)):
        ratio_min = float(np.nanmin(ratio_t0))
        if np.isfinite(ratio_min) and ratio_min < 1.0:
            ax_ratio.set_ylim(ratio_min, 1.0)
        else:
            ax_ratio.set_ylim(0.0, 1.0)
    else:
        ax_ratio.set_ylim(0.0, 1.0)
    ax_ratio.grid(True, alpha=0.3)

    if time_axis.size > 0 and np.any(np.isfinite(time_axis)):
        t_max = float(np.nanmax(time_axis))
        if start_year_at_zero:
            t_min = 0.0
            if t_max <= 0.0:
                t_max = 1.0
        else:
            t_min = float(np.nanmin(time_axis))
            if not np.isfinite(t_min) or t_max <= t_min:
                t_min, t_max = 0.0, 1.0
        for ax in axes:
            ax.set_xlim(t_min, t_max)

    for ax in axes:
        apply_axes_font(ax, font)

    out_dir = Path(output_root) / "plots" / "python" / "stratigraphy" / "domain_preservation_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "domain_preservation_metrics.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return {
        "metrics_path": out_path,
        "output_dir": out_dir,
        "highlight_index": highlight_idx,
        "highlight_tag": highlight_tag,
    }


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


def _resolve_bathy_colormap(cmap_name: str):
    """Resolve bathymetry colormap names from project and SedTRAILS palettes."""
    if cmap_name.endswith(".clrmap"):
        return load_clrmap_file(cmap_name), None
    try:
        return bathymetry_colormap(cmap_name)
    except ValueError:
        return get_named_colormap(cmap_name), None


def _bathy_contour_levels(norm, fallback: tuple[float, float] = (-20.0, 10.0), step: float = 0.2) -> np.ndarray:
    """Return contour levels matching a bathymetry colormap normalization."""
    vmin = getattr(norm, "vmin", None)
    vmax = getattr(norm, "vmax", None)
    if vmin is None or vmax is None or not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = fallback
    return np.arange(float(vmin), float(vmax) + 0.5 * step, step)


def _coerce_bathy_cube(source: str | Path | BathyCube) -> BathyCube:
    """Return a BathyCube from a path or an already-converted cube."""
    if isinstance(source, BathyCube):
        return source
    return load_bathy_cube(source)


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
    base = min_elev - 0.01 * abs(min_elev)
    # Clip to the baseline to avoid extra whitespace below the section.
    y_min_plot = max(base, y_min_plot)
    y_max_plot = max(y_max_plot, mhw + 0.05 * abs(mhw))
    y_range = y_max_plot - y_min_plot
    if not np.isfinite(y_range) or y_range <= 0:
        y_range = 1.0

    x_range = float(np.nanmax(x_m) - np.nanmin(x_m)) if len(x_m) else 0.0
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
        ax.set_xlim(0.0, float(x_m[-1]) if len(x_m) else 0.0)
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
    bathy_nc_path: str | Path | BathyCube,
    output_root: str | Path,
    transect_rows_km: np.ndarray,
    tick_spacing_m: float = 100.0,
    mhw: float = 0.358,
    mlw: float = -0.590,
    highlight_date: str | np.datetime64 | list[str | np.datetime64] | None = None,
    n_points: int = 400,
    initial_index: int = 0,
    dx: float | None = 20.0,
    target_crs: str = "EPSG:32618",
    max_fig_width_cm: float = 20.0,
    max_fig_height_cm: float = 5.0,
    plot_xlim: tuple[float, float] | None = None,
    plot_ylim: tuple[float, float] | None = None,
    clip_x_to_data: bool = True,
    clip_y_to_data: bool = True,
    make_gif: bool = False,
    gif_fps: int = 6,
    gif_stride: int = 1,
) -> dict[str, object]:
    """Plot 6-panel stratigraphy summaries for a list of transect slices.

    Args:
        bathy_nc_path: Path to a bathymetry cube file or an already-converted BathyCube.
        output_root: Output root directory for plots.
        transect_rows_km: Array of transect endpoints in kilometers.
        tick_spacing_m: Tick spacing along transects (meters).
        mhw: Mean high water elevation.
        mlw: Mean low water elevation.
        highlight_date: Optional date string(s) or datetime64(s) to highlight deposits.
        n_points: Number of points used for each transect slice.
        initial_index: Initial stratigraphy index.
        dx: Grid spacing in meters. If None, use each extracted transect's sample spacing.
        target_crs: CRS string for plot metadata.
        max_fig_width_cm: Maximum plot width for scaled sections (cm).
        max_fig_height_cm: Maximum plot height for scaled sections (cm).
        plot_xlim: Optional x-axis limits for section plots in along-transect meters.
        plot_ylim: Optional elevation limits for section plots in meters.
        clip_x_to_data: Clip section x-limits to data extent when plot_xlim is not provided.
        clip_y_to_data: Clip section y-limits to data extent when plot_ylim is not provided.
        make_gif: Whether to create a stratigraphy GIF for each transect.
        gif_fps: Frames per second for GIF export.
        gif_stride: Step between frames (e.g., 2 uses every other survey).

    Returns:
        Dict containing the bathy cube, labels, and output directory for reuse.
    """
    labels = [_label_for_index(i) for i in range(len(transect_rows_km))]
    output_root = Path(output_root)
    bathy_cube = _coerce_bathy_cube(bathy_nc_path)
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
        slice_dx = float(dx) if dx is not None else float(np.nanmedian(np.diff(slice_cube.x[0, :])))
        slice_cfg = StratigraphyConfig(initial_index=initial_index, dx=slice_dx, target_crs=target_crs)
        slice_result = compute_stratigraphy(slice_cube, slice_cfg)

        z_min = min(z_min, np.nanmin(slice_cube.z))
        z_max = max(z_max, np.nanmax(slice_cube.z))
        entries.append((label, transect, slice_cube, slice_result, length_m))

    resolved_plot_ylim = plot_ylim
    max_depth_range = None
    if resolved_plot_ylim is None and np.isfinite(z_min) and np.isfinite(z_max) and z_max > z_min:
        resolved_plot_ylim = (z_min - 0.01 * abs(z_min), z_max + 0.01 * abs(z_max))
    if resolved_plot_ylim is not None:
        max_depth_range = resolved_plot_ylim[1] - resolved_plot_ylim[0]

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
                plot_xlim=plot_xlim,
                plot_ylim=resolved_plot_ylim,
                title_prefix=section_title,
            )
            plot_stacked_stratigraphy_section(
                slice_cube,
                slice_dir / f"strat_section_{label}.png",
                f"Cross-section {label}-{label}'",
                plot_xlim=plot_xlim,
                plot_ylim=resolved_plot_ylim,
                x_ticks=xticks_m,
                mhw=mhw,
                mlw=mlw,
                highlight_date=None,
                max_transect_length=max_len_m,
                max_depth_range=max_depth_range,
                max_fig_width_cm=max_fig_width_cm,
                max_fig_height_cm=max_fig_height_cm,
                clip_x_to_data=clip_x_to_data,
                clip_y_to_data=clip_y_to_data,
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
                    plot_xlim=plot_xlim,
                    plot_ylim=resolved_plot_ylim,
                    highlight_date=date_val,
                    title_prefix=section_title,
                )
                plot_stacked_stratigraphy_section(
                    slice_cube,
                    slice_highlight_dir / f"strat_section_{label}_highlight_{highlight_tag}.png",
                    f"Cross-section {label}-{label}'",
                    plot_xlim=plot_xlim,
                    plot_ylim=resolved_plot_ylim,
                    x_ticks=xticks_m,
                    mhw=mhw,
                    mlw=mlw,
                    highlight_date=date_val,
                    max_transect_length=max_len_m,
                    max_depth_range=max_depth_range,
                    max_fig_width_cm=max_fig_width_cm,
                    max_fig_height_cm=max_fig_height_cm,
                    clip_x_to_data=clip_x_to_data,
                    clip_y_to_data=clip_y_to_data,
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
                    plot_ylim=resolved_plot_ylim,
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
                        plot_ylim=resolved_plot_ylim,
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
    cmap_name: str = "SEAWAD",
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
    cmap, norm = _resolve_bathy_colormap(cmap_name)
    levels = _bathy_contour_levels(norm)
    cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, norm=norm, extend="both")
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
        cmap, norm = _resolve_bathy_colormap(cmap_name)
        levels = _bathy_contour_levels(norm)
        cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, norm=norm, extend="both")
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
    cmap_name: str = "SEAWAD",
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
    cmap, norm = _resolve_bathy_colormap(cmap_name)
    levels = _bathy_contour_levels(norm)
    cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, norm=norm, extend="both")
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
        cmap, norm = _resolve_bathy_colormap(cmap_name)
        levels = _bathy_contour_levels(norm)
        cf = ax.contourf(x_km, y_km, z_last, levels=levels, cmap=cmap, norm=norm, extend="both")
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

    for label, transect, slice_cube, slice_result, length_m in entries:
        section_title = f"{label}-{label}'"
        xticks_m = np.arange(0.0, length_m + 1e-6, tick_spacing_m)

        if not highlight_dates:
            plot_stratigraphy_stack(
                slice_cube,
                slice_result,
                slice_dir / f"strat_overview_section_{label}.png",
                f"overview_{label}",
                time_vals=slice_cube.t,
                title_prefix=section_title,
            )
            plot_stacked_stratigraphy_section(
                slice_cube,
                slice_dir / f"strat_section_{label}.png",
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
                    slice_highlight_dir / f"strat_{label}_highlight_{highlight_tag}.png",
                    f"overview_slice_{label}",
                    time_vals=slice_cube.t,
                    highlight_date=date_val,
                    title_prefix=section_title,
                )
                plot_stacked_stratigraphy_section(
                    slice_cube,
                    slice_highlight_dir / f"strat_section_{label}_highlight_{highlight_tag}.png",
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

    if valid_transects:
        map_path = plot_transect_location_plan_from_transects(
            bathy_cube,
            valid_transects,
            output_root=output_root,
            labels=[t.name for t in valid_transects],
            mlw=mlw,
            tick_spacing_m=tick_spacing_m,
            tick_length_km=map_tick_length_km,
            cmap_name="SEAWAD",
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
