from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from bathy_formatter import (
    RawSurvey,
    _datenum_to_date,
    _matlab_datenum,
    compute_domain_extents,
    compute_minimum_overlap_mask,
    process_bathy_formatter,
    regrid_surveys,
)


def _make_raw_survey(
    location: str,
    year: int,
    month: int,
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    z_vals: np.ndarray,
) -> RawSurvey:
    # Synthetic survey builder for small, deterministic grids.
    x_raw, y_raw = np.meshgrid(x_vals, y_vals)
    return RawSurvey(
        location=location,
        datenum=_matlab_datenum(year, month, 1),
        vertical_datum="NAVD88",
        horizontal_datum=None,
        ncfile=f"{location}_{year:04d}_{month:02d}_NAVD88.nc",
        x_raw=x_raw,
        y_raw=y_raw,
        z_raw=z_vals,
    )


def _write_nc(path: Path, x_vals: np.ndarray, y_vals: np.ndarray, z_vals: np.ndarray) -> None:
    # Minimal NetCDF with x/y/z variables to exercise loader paths.
    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("x", len(x_vals))
        nc.createDimension("y", len(y_vals))
        x_var = nc.createVariable("x", "f8", ("x",))
        y_var = nc.createVariable("y", "f8", ("y",))
        z_var = nc.createVariable("z", "f8", ("y", "x"))
        x_var[:] = x_vals
        y_var[:] = y_vals
        z_var[:] = z_vals


def test_datenum_roundtrip_simple_date() -> None:
    # Date conversion should preserve year/month/day for whole-day datenums.
    dn = _matlab_datenum(2020, 5, 15)
    y, m, d = _datenum_to_date(dn)
    assert (y, m, d) == (2020, 5, 15)


def test_compute_domain_extents_ignores_nans() -> None:
    # Overlap extents should be based on valid data, not NaNs.
    x_vals = np.array([0.0, 1000.0, 2000.0], dtype=float)
    y_vals = np.array([0.0, 1000.0], dtype=float)
    z = np.array([[np.nan, -1.0, -2.0], [np.nan, -1.0, -2.0]], dtype=float)
    s1 = _make_raw_survey("Unit", 2005, 1, x_vals, y_vals, z)

    z2 = np.array([[-3.0, -2.0, -1.0], [-3.0, -2.0, -1.0]], dtype=float)
    s2 = _make_raw_survey("Unit", 2006, 1, x_vals + 1000.0, y_vals + 1000.0, z2)

    x_lims, y_lims, x_min, y_min = compute_domain_extents([s1, s2])

    assert np.allclose(x_lims, [1.0, 3.0])
    assert np.allclose(y_lims, [0.0, 2.0])
    assert np.allclose(x_min, [1.0, 2.0])
    assert np.allclose(y_min, [1.0, 1.0])


def test_regrid_surveys_returns_overlap_mask() -> None:
    # Regridding should return the overlap mask and use it to NaN-out invalid cells.
    x_vals = np.arange(0.0, 80.0, 20.0)
    y_vals = np.arange(0.0, 80.0, 20.0)
    x_raw, y_raw = np.meshgrid(x_vals, y_vals)
    z = -1.0 * np.ones_like(x_raw)

    surveys = [
        _make_raw_survey("Unit", 2005, 1, x_vals, y_vals, z),
        _make_raw_survey("Unit", 2006, 1, x_vals, y_vals, z + 0.1),
    ]

    x_lims, y_lims, x_min, y_min = compute_domain_extents(surveys)
    bathy, overlap_mask = regrid_surveys(
        surveys,
        dx=20.0,
        x_min=x_min,
        y_min=y_min,
        overlap_erosion_cells=1,
    )

    assert overlap_mask.shape == bathy.x.shape
    assert np.any(~overlap_mask)
    for i in range(bathy.z.shape[2]):
        assert np.all(np.isnan(bathy.z[:, :, i])[~overlap_mask])


def test_process_bathy_formatter_marks_dropped_dates(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # Dropped surveys should be flagged in printed date list and excluded from processing.
    nc_dir = tmp_path / "nc"
    out_dir = tmp_path / "out"
    nc_dir.mkdir()
    out_dir.mkdir()

    x_vals = np.array([0.0, 20.0], dtype=float)
    y_vals = np.array([0.0, 20.0], dtype=float)
    z_vals = np.array([[-1.0, -1.1], [-1.2, -1.3]], dtype=float)

    _write_nc(nc_dir / "Unit_2007_01_NAVD88.nc", x_vals, y_vals, z_vals)
    _write_nc(nc_dir / "Unit_2008_01_NAVD88.nc", x_vals, y_vals, z_vals)

    result = process_bathy_formatter(
        nc_dir=nc_dir,
        out_dir=out_dir,
        dx=20.0,
        drop_survey=["2007-01-01"],
        print_raw_survey_dates=True,
        output_format="netcdf",
    )

    captured = capsys.readouterr()
    assert "Raw survey dates:" in captured.out
    assert "2007-01-01 [dropped]" in captured.out
    assert "2008-01-01" in captured.out
    assert len(result.surveys_processed) == 1


def test_compute_minimum_overlap_mask_nonempty() -> None:
    # Overlap mask for a single valid survey should be non-empty.
    x_vals = np.arange(0.0, 60.0, 20.0)
    y_vals = np.arange(0.0, 60.0, 20.0)
    x_raw, y_raw = np.meshgrid(x_vals, y_vals)
    z = -1.0 * np.ones_like(x_raw)

    survey = _make_raw_survey("Unit", 2005, 1, x_vals, y_vals, z)
    grid_x, grid_y = np.meshgrid(x_vals, y_vals)

    mask = compute_minimum_overlap_mask([survey], grid_x, grid_y)
    assert mask.shape == grid_x.shape
    assert np.any(mask)
