"""Tests for notebook workflow helpers."""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

import bathymetry_analysis.workflows as workflows
from bathymetry_analysis.stratigraphy import BathyCube


def _matlab_datenum(dt: datetime) -> float:
    """Convert a Python datetime to MATLAB datenum."""
    midnight = datetime(dt.year, dt.month, dt.day)
    return dt.toordinal() + 366 + (dt - midnight).total_seconds() / 86400.0


def test_startup_check_no_install(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure startup checks skip installation for empty requirements."""
    called = []

    def fake_check_call(*_args, **_kwargs):
        called.append(True)

    monkeypatch.setattr(workflows.subprocess, "check_call", fake_check_call)
    info = workflows.startup_check(required_modules={})

    assert "code_dir" in info
    assert called == []


def test_run_morphodynamics_vs_wave_power_writes_plots(tmp_path: Path) -> None:
    """Write stacked wave power plots for synthetic inputs."""
    morpho_csv = tmp_path / "morpho.csv"
    wave_stats = tmp_path / "wave_stats.txt"
    plot_dir = tmp_path / "plots"

    dates = [datetime(2020, 1, 15), datetime(2020, 2, 15)]
    morpho_df = pd.DataFrame(
        {
            "scope": ["step", "step"],
            "time_datenum": [_matlab_datenum(d) for d in dates],
            "erosion_m3": [1.0, 2.0],
            "accretion_m3": [3.0, 4.0],
            "gross_m3": [4.0, 6.0],
            "net_m3": [2.0, 2.0],
        }
    )
    morpho_df.to_csv(morpho_csv, index=False)

    wave_df = pd.DataFrame(
        {
            "buoy": ["B1", "B1"],
            "start_date": ["2019-12-15", "2020-01-15"],
            "end_date": ["2020-01-15", "2020-02-15"],
            "period_days": [31, 31],
            "cum_wave_power_MWh_m": [10.0, 20.0],
        }
    )
    wave_df.to_csv(wave_stats, sep="\t", index=False)

    merged = workflows.run_morphodynamics_vs_wave_power(
        wave_power_stats_file=wave_stats,
        morpho_csv_path=morpho_csv,
        stacked_plot_dir=plot_dir,
    )

    assert not merged.empty
    assert (plot_dir / "morpho_vs_wave_power_B1.png").exists()
    assert (plot_dir / "morpho_vs_wave_power_per_day_B1.png").exists()


def test_run_morpho_wave_correlations_creates_outputs(tmp_path: Path) -> None:
    """Generate correlation plots from a minimal merged dataframe."""
    merged = pd.DataFrame(
        {
            "buoy": ["B1", "B1"],
            "start_date": ["2019-12-15", "2020-01-15"],
            "end_date": ["2020-01-15", "2020-02-15"],
            "period_days": [31, 31],
            "cum_wave_power_MWh_m": [10.0, 20.0],
            "erosion_m3": [1.0, 2.0],
            "accretion_m3": [3.0, 4.0],
            "gross_m3": [4.0, 6.0],
            "net_m3": [2.0, 2.0],
        }
    )

    output_dir = tmp_path / "corr"
    metrics_map = {
        "erosion_m3": "Erosion",
        "net_m3": "Net",
    }
    x_map = {
        "cum_wave_power_MWh_m": "Total wave power [MWh/m]",
        "period_days": "Interval duration [days]",
    }

    workflows.run_morpho_wave_correlations(
        merged,
        output_dir,
        highlight_dates=[],
        outlier_dates=[],
        metrics_map=metrics_map,
        x_map=x_map,
    )

    assert (output_dir / "corr_B1_cum_wave_power_MWh_m_erosion_m3.png").exists()
    assert (output_dir / "correlation_compilation_cum_wave_power_MWh_m.png").exists()
    assert (output_dir / "correlation_compilation_interval_duration.png").exists()
    assert (output_dir / "correlation_compilation_rate_cum_wave_power_MWh_m_per_day.png").exists()


def test_run_transect_slice_plots_accepts_bathy_cube(tmp_path: Path) -> None:
    """Run the common transect workflow from an already-converted BathyCube."""
    x_axis = np.linspace(0.0, 100.0, 6)
    y_axis = np.linspace(0.0, 40.0, 3)
    x, y = np.meshgrid(x_axis, y_axis)
    t = np.array([_matlab_datenum(datetime(2020, 1, 1)), _matlab_datenum(datetime(2020, 2, 1))])
    z = np.empty((len(y_axis), len(x_axis), len(t)), dtype=float)
    z[:, :, 0] = -2.0 + 0.01 * x
    z[:, :, 1] = z[:, :, 0] + 0.2
    cube = BathyCube(location="WorkflowCube", t=t, x=x, y=y, z=z)

    outputs = workflows.run_transect_slice_plots(
        bathy_nc_path=cube,
        output_root=tmp_path,
        transect_rows_km=np.array([[0.0, 0.02, 0.1, 0.02]], dtype=float),
        n_points=20,
        dx=None,
        plot_xlim=(10.0, 80.0),
        plot_ylim=(-3.0, 0.0),
        clip_x_to_data=False,
        clip_y_to_data=False,
    )

    assert outputs["bathy_cube"] is cube
    assert (outputs["slice_dir"] / "strat_section_A.png").exists()
