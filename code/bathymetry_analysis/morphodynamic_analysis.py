"""Morphodynamic change analysis for regridded bathymetry."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bathy import BathyGrid, BathyProcessResult, load_for_analysis


@dataclass
class MorphodynamicResult:
    """Container for morphodynamic change outputs.

    :ivar bathy: Source bathymetry grid.
    :ivar change_grids: 3D array of stepwise elevation changes.
    :ivar step_labels: Human-readable step labels (start -> end).
    :ivar step_times: MATLAB datenums for step endpoints.
    :ivar rel_times: MATLAB datenums for relative-to-first steps.
    :ivar step_stats: DataFrame of stepwise volume metrics.
    :ivar rel_stats: DataFrame of cumulative volume metrics.
    :ivar grid_area_m2: Grid cell area in square meters.
    :ivar vmax: Symmetric plotting limit for change maps.
    :ivar overlap_mask: Optional overlap mask used for plotting.
    """
    bathy: BathyGrid
    change_grids: np.ndarray
    step_labels: list[str]
    step_times: np.ndarray
    rel_times: np.ndarray
    step_stats: pd.DataFrame
    rel_stats: pd.DataFrame
    grid_area_m2: float
    vmax: float
    overlap_mask: np.ndarray | None


def _datenum_to_year_month(datenum_value: float) -> tuple[int, int]:
    """Convert MATLAB datenum to year and month.

    Args:
        datenum_value: MATLAB datenum float.

    Returns:
        Tuple of ``(year, month)``.
    """
    value = float(np.asarray(datenum_value).reshape(-1)[0])
    dt = datetime.fromordinal(int(value)) + timedelta(days=value % 1) - timedelta(days=366)
    return dt.year, dt.month


def _datenum_to_datetime(datenum_value: float) -> datetime:
    """Convert MATLAB datenum to ``datetime``.

    Args:
        datenum_value: MATLAB datenum float.

    Returns:
        Converted ``datetime``.
    """
    value = float(np.asarray(datenum_value).reshape(-1)[0])
    return datetime.fromordinal(int(value)) + timedelta(days=value % 1) - timedelta(days=366)


def _format_step_label(t0: float, t1: float) -> str:
    """Format a date interval label from two MATLAB datenums."""
    y0, m0 = _datenum_to_year_month(t0)
    y1, m1 = _datenum_to_year_month(t1)
    return f"{y0:04d}-{m0:02d} -> {y1:04d}-{m1:02d}"


def _format_time_label(t_value: float) -> str:
    """Format a year-month label from a MATLAB datenum."""
    y, m = _datenum_to_year_month(t_value)
    return f"{y:04d}-{m:02d}"


def _grid_spacing_from_mesh(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Infer grid spacing from 1D or 2D coordinate arrays.

    Args:
        x: X coordinates (1D or meshgrid).
        y: Y coordinates (1D or meshgrid).

    Returns:
        Tuple of ``(dx, dy)`` in coordinate units.

    Raises:
        ValueError: If spacing cannot be determined.
    """
    if x.ndim == 2:
        dx_vals = np.diff(x[0, :])
    else:
        dx_vals = np.diff(x)
    if y.ndim == 2:
        dy_vals = np.diff(y[:, 0])
    else:
        dy_vals = np.diff(y)
    dx = float(np.nanmedian(np.abs(dx_vals)))
    dy = float(np.nanmedian(np.abs(dy_vals)))
    if not np.isfinite(dx) or dx <= 0:
        raise ValueError("Could not determine valid grid spacing in x")
    if not np.isfinite(dy) or dy <= 0:
        raise ValueError("Could not determine valid grid spacing in y")
    return dx, dy


def _compute_volume_stats(delta: np.ndarray, cell_area: float) -> tuple[float, float, float, float]:
    """Compute erosion/accretion/gross/net volumes for a change grid.

    Args:
        delta: 2D change grid.
        cell_area: Grid cell area in square meters.

    Returns:
        Tuple of ``(erosion, accretion, gross, net)`` in cubic meters.
    """
    valid = np.isfinite(delta)
    if not np.any(valid):
        return 0.0, 0.0, 0.0, 0.0
    delta_valid = delta[valid]
    accretion = float(np.sum(delta_valid[delta_valid > 0.0]) * cell_area)
    erosion = float(-np.sum(delta_valid[delta_valid < 0.0]) * cell_area)
    gross = float(np.sum(np.abs(delta_valid)) * cell_area)
    net = float(np.sum(delta_valid) * cell_area)
    return erosion, accretion, gross, net


def compute_morphodynamic_stats(
    source: BathyGrid | BathyProcessResult | str | Path,
    percentile: float = 98.0,
) -> MorphodynamicResult:
    """Compute morphodynamic change grids and volume statistics.

    Args:
        source: Bathymetry source (grid, process result, or file path).
        percentile: Percentile for symmetric change-map scaling.

    Returns:
        ``MorphodynamicResult`` with change grids and summary tables.

    Raises:
        ValueError: If the bathymetry stack is invalid.
    """
    overlap_mask = None
    if isinstance(source, BathyProcessResult):
        overlap_mask = source.overlap_mask
    elif hasattr(source, "overlap_mask"):
        overlap_mask = getattr(source, "overlap_mask")

    bathy = load_for_analysis(source)
    z = np.asarray(bathy.z, dtype=float)
    if z.ndim != 3 or z.shape[-1] < 2:
        raise ValueError("bathy.z must be a 3D array with at least two timesteps")

    t_vals = np.asarray(bathy.t, dtype=float).reshape(-1)
    if t_vals.size != z.shape[-1]:
        raise ValueError("bathy.t length does not match bathy.z time dimension")

    dx, dy = _grid_spacing_from_mesh(bathy.x, bathy.y)
    cell_area = dx * dy

    change_grids = z[:, :, 1:] - z[:, :, :-1]
    nan_mask = ~np.isfinite(z[:, :, 1:]) | ~np.isfinite(z[:, :, :-1])
    change_grids = np.where(nan_mask, np.nan, change_grids)

    abs_change = np.abs(change_grids)
    vmax = float(np.nanpercentile(abs_change, percentile)) if np.any(np.isfinite(abs_change)) else 0.0
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0

    step_rows = []
    step_labels: list[str] = []
    step_times = t_vals[1:].copy()

    for i in range(change_grids.shape[-1]):
        delta = change_grids[:, :, i]
        erosion, accretion, gross, net = _compute_volume_stats(delta, cell_area)
        label = _format_step_label(t_vals[i], t_vals[i + 1])
        step_labels.append(label)
        step_rows.append(
            {
                "scope": "step",
                "index": i,
                "time_datenum": float(step_times[i]),
                "time_label": label,
                "erosion_m3": erosion,
                "accretion_m3": accretion,
                "gross_m3": gross,
                "net_m3": net,
            }
        )

    rel_rows = []
    rel_times = t_vals.copy()
    step_gross = np.array([row["gross_m3"] for row in step_rows], dtype=float)
    cumulative_gross = np.concatenate(([0.0], np.cumsum(step_gross)))
    for i in range(z.shape[-1]):
        delta_rel = z[:, :, i] - z[:, :, 0]
        nan_mask_rel = ~np.isfinite(z[:, :, i]) | ~np.isfinite(z[:, :, 0])
        delta_rel = np.where(nan_mask_rel, np.nan, delta_rel)
        erosion, accretion, _gross, net = _compute_volume_stats(delta_rel, cell_area)
        label = _format_time_label(t_vals[i])
        rel_rows.append(
            {
                "scope": "relative_to_first",
                "index": i,
                "time_datenum": float(rel_times[i]),
                "time_label": label,
                "erosion_m3": erosion,
                "accretion_m3": accretion,
                "gross_m3": float(cumulative_gross[i]),
                "net_m3": net,
            }
        )

    step_stats = pd.DataFrame(step_rows)
    rel_stats = pd.DataFrame(rel_rows)

    return MorphodynamicResult(
        bathy=bathy,
        change_grids=change_grids,
        step_labels=step_labels,
        step_times=step_times,
        rel_times=rel_times,
        step_stats=step_stats,
        rel_stats=rel_stats,
        grid_area_m2=cell_area,
        vmax=vmax,
        overlap_mask=overlap_mask,
    )


def save_morphodynamic_csv(result: MorphodynamicResult, csv_path: str | Path) -> Path:
    """Save combined step and cumulative statistics to CSV.

    Args:
        result: Morphodynamic result to serialize.
        csv_path: Output CSV path.

    Returns:
        Path to the saved CSV.
    """
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    combined = pd.concat([result.step_stats, result.rel_stats], ignore_index=True)
    combined.to_csv(csv_path, index=False)
    return csv_path


def plot_change_maps(
    result: MorphodynamicResult,
    plot_dir: str | Path,
    cmap_name: str = "RdBu_r",
    levels: int = 41,
) -> list[Path]:
    """Plot stepwise change maps for each survey interval.

    Args:
        result: Morphodynamic result to plot.
        plot_dir: Output directory for images.
        cmap_name: Matplotlib colormap name.
        levels: Number of contour levels.

    Returns:
        List of saved plot paths.
    """
    out_dir = Path(plot_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    x_km = result.bathy.x / 1000.0
    y_km = result.bathy.y / 1000.0

    level_values = np.linspace(-result.vmax, result.vmax, levels)
    paths: list[Path] = []

    t_vals = np.asarray(result.bathy.t, dtype=float).reshape(-1)

    for i, label in enumerate(result.step_labels):
        fig, ax = plt.subplots(figsize=(12, 8))
        z_plot = np.ma.masked_invalid(result.change_grids[:, :, i])
        cf = ax.contourf(x_km, y_km, z_plot, levels=level_values, cmap=cmap_name, extend="both")

        ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=14, weight="bold")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        if result.overlap_mask is not None:
            ax.contour(
                x_km,
                y_km,
                result.overlap_mask.astype(float),
                levels=[0.5],
                colors=[(0.0, 0.0, 0.0)],
                linewidths=1.2,
            )

        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label("Elevation change [m]")

        fig.tight_layout()

        y0, m0 = _datenum_to_year_month(t_vals[i])
        y1, m1 = _datenum_to_year_month(t_vals[i + 1])
        fname = f"{result.bathy.location}_MorphChange_{y0:04d}_{m0:02d}_to_{y1:04d}_{m1:02d}_py.png"
        out_path = out_dir / fname
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        paths.append(out_path)

    return paths


def plot_volume_timeseries(
    result: MorphodynamicResult,
    plot_dir: str | Path,
    filename: str = "morphodynamics_volume_timeseries.png",
) -> Path:
    """Plot cumulative and stepwise volume change time series.

    Args:
        result: Morphodynamic result to plot.
        plot_dir: Output directory for the image.
        filename: Output filename.

    Returns:
        Path to the saved plot.
    """
    out_dir = Path(plot_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    import matplotlib.dates as mdates

    rel_stats = result.rel_stats.copy()
    step_stats = result.step_stats.copy()

    t_vals = np.asarray(result.bathy.t, dtype=float).reshape(-1)
    rel_dates = [_datenum_to_datetime(t) for t in result.rel_times]
    step_dates = [_datenum_to_datetime(t) for t in result.step_times]

    rel_nums = mdates.date2num(rel_dates)
    step_nums = mdates.date2num(step_dates)

    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(12, 8), sharex=False)

    ax_top.plot(rel_nums, rel_stats["gross_m3"], marker="o", color="#4477AA", label="Gross vs first")
    ax_top.plot(rel_nums, rel_stats["net_m3"], marker="o", color="#228833", label="Net vs first")
    ax_top.set_ylabel("Volume change vs first [m^3]")
    ax_top.grid(True, alpha=0.3)
    ax_top.set_title("Morphodynamic volume change")
    ax_top.legend()
    ax_top.xaxis_date()
    ax_top.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))

    if len(step_nums) > 1:
        median_step = float(np.nanmedian(np.diff(step_nums)))
    else:
        median_step = 30.0

    bar_width = 0.6 * median_step
    offset = 0.18 * median_step
    net_width = 0.28 * median_step

    erosion_plot = -step_stats["erosion_m3"].to_numpy()
    accretion_plot = step_stats["accretion_m3"].to_numpy()
    net_plot = step_stats["net_m3"].to_numpy()

    ax_bottom.bar(step_nums - offset, erosion_plot, width=bar_width * 0.6, label="Erosion", color="#4477AA")
    ax_bottom.bar(step_nums + offset, accretion_plot, width=bar_width * 0.6, label="Accretion", color="#EE6677")
    ax_bottom.bar(step_nums, net_plot, width=net_width, label="Net", color="#777777", alpha=0.6)
    ax_bottom.set_ylabel("Volume [m^3]")
    ax_bottom.grid(True, axis="y", alpha=0.3)
    ax_bottom.legend()
    ax_bottom.xaxis_date()
    ax_bottom.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    fig.autofmt_xdate(rotation=45, ha="right")

    fig.tight_layout()
    out_path = out_dir / filename
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return out_path


def run_morphodynamic_analysis(
    source: BathyGrid | BathyProcessResult | str | Path,
    plot_dir: str | Path,
    csv_path: str | Path | None = None,
    cmap_name: str = "RdBu_r",
    percentile: float = 98.0,
) -> MorphodynamicResult:
    """Run the full morphodynamic workflow and write plots/CSV.

    Args:
        source: Bathymetry source (grid, process result, or file path).
        plot_dir: Output directory for plots.
        csv_path: Optional CSV output path.
        cmap_name: Matplotlib colormap name.
        percentile: Percentile for change-map scaling.

    Returns:
        Morphodynamic result with computed statistics.
    """
    result = compute_morphodynamic_stats(source, percentile=percentile)
    plot_morphodynamic_results(result, plot_dir, csv_path=csv_path, cmap_name=cmap_name)
    return result


def plot_morphodynamic_results(
    result: MorphodynamicResult,
    plot_dir: str | Path,
    csv_path: str | Path | None = None,
    cmap_name: str = "RdBu_r",
) -> Path:
    """Write morphodynamic plots and summary CSV.

    Args:
        result: Morphodynamic result to serialize/plot.
        plot_dir: Output directory for plots.
        csv_path: Optional CSV output path.
        cmap_name: Matplotlib colormap name.

    Returns:
        Path to the written CSV file.
    """
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    csv_out = csv_path or (plot_dir / "morphodynamics_summary.csv")
    save_morphodynamic_csv(result, csv_out)
    plot_change_maps(result, plot_dir, cmap_name=cmap_name)
    plot_volume_timeseries(result, plot_dir)

    return csv_out
