"""Notebook-friendly workflows for bathymetry analysis."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
import importlib
import subprocess
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils.date_utils import interval_highlight_mask

DEFAULT_REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "netCDF4": "netCDF4",
    "pyproj": "pyproj",
    "geopandas": "geopandas",
    "joblib": "joblib",
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
    print(f"Saved compilation plot: {comp_path}")


def run_morpho_wave_correlations(
    merged_df: pd.DataFrame,
    output_dir: str | Path,
    highlight_dates=None,
    outlier_dates=None,
    metrics_map: dict[str, str] | None = None,
    x_map: dict[str, str] | None = None,
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
    )

    plot_correlation_panels(
        merged_df,
        buoys,
        metrics_rate,
        x_vars_rate,
        output_dir,
        highlight_dates=highlight_dates,
        outlier_dates=outlier_dates,
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
        )
