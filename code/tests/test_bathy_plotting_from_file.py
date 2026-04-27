from __future__ import annotations

from pathlib import Path
from datetime import datetime, timedelta

import matplotlib
import numpy as np
import pytest

from bathy_formatter import BathyGrid, plot_bathy_formatter_outputs, save_bathy_grid_mat, save_bathy_grid_netcdf

matplotlib.use("Agg")


def _sample_bathy() -> BathyGrid:
    x_1d = np.array([1000.0, 1020.0, 1040.0], dtype=float)
    y_1d = np.array([2000.0, 2020.0], dtype=float)
    x_grid, y_grid = np.meshgrid(x_1d, y_1d)

    # MATLAB-style datenums as column vector.
    t = np.array([[738000.0], [738001.0]], dtype=float)

    # z is stored as (y, x, time).
    z = np.empty((2, 3, 2), dtype=float)
    z[:, :, 0] = np.array([[-1.0, -0.8, -0.6], [-1.2, -1.0, -0.9]], dtype=float)
    z[:, :, 1] = np.array([[-0.9, -0.7, -0.5], [-1.1, -0.95, -0.85]], dtype=float)

    return BathyGrid(location="UnitTestInlet", t=t, x=x_grid, y=y_grid, z=z)


def _datenum_to_year_month(datenum_value: float) -> tuple[int, int]:
    dt = datetime.fromordinal(int(datenum_value)) + timedelta(days=datenum_value % 1) - timedelta(days=366)
    return dt.year, dt.month


@pytest.mark.parametrize("fmt", ["mat", "nc"])
def test_plot_bathy_formatter_outputs_from_file_path(tmp_path: Path, capsys: pytest.CaptureFixture[str], fmt: str) -> None:
    bathy = _sample_bathy()
    plot_dir = tmp_path / "plots"

    if fmt == "mat":
        data_path = tmp_path / "bathy_test.mat"
        save_bathy_grid_mat(bathy, data_path)
    else:
        data_path = tmp_path / "bathy_test.nc"
        save_bathy_grid_netcdf(bathy, data_path)

    plot_bathy_formatter_outputs(
        process_result=data_path,
        plot_dir=plot_dir,
        bathy_cmap_name="kg2",
        extent_boundary_method="mask",
        mask_nan_plots=True,
    )

    captured = capsys.readouterr()
    assert "Loaded from file path: skipping raw extents/raw survey plots" in captured.out

    raw_dir = plot_dir / "raw"
    regridded_dir = plot_dir / "regridded"

    assert raw_dir.exists()
    assert regridded_dir.exists()

    # Raw survey plots are not produced in file-path mode.
    assert not list(raw_dir.glob("*.png"))

    # Extents plot is not produced in file-path mode.
    assert not list(plot_dir.glob("*_bathyExtents_*_py.png"))

    expected_plot_count = len({_datenum_to_year_month(float(v)) for v in bathy.t.flatten()})
    regrid_plots = list(regridded_dir.glob("*_ReGridBathy_*_py.png"))
    assert len(regrid_plots) == expected_plot_count
