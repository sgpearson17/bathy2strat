from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from bathy_formatter import (
    BathyGrid,
    BathyProcessResult,
    load_bathy_grid,
    load_bathy_grid_mat,
    load_bathy_grid_netcdf,
    load_for_analysis,
    save_bathy_grid_mat,
    save_bathy_grid_netcdf,
)


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


def _assert_bathy_equal(expected: BathyGrid, actual: BathyGrid) -> None:
    assert actual.location == expected.location
    np.testing.assert_allclose(actual.t, expected.t)
    np.testing.assert_allclose(actual.x, expected.x)
    np.testing.assert_allclose(actual.y, expected.y)
    np.testing.assert_allclose(actual.z, expected.z)


def test_mat_roundtrip_loaders(tmp_path: Path) -> None:
    bathy = _sample_bathy()
    mat_path = tmp_path / "bathy_test.mat"

    save_bathy_grid_mat(bathy, mat_path)

    loaded_specific = load_bathy_grid_mat(mat_path)
    loaded_generic = load_bathy_grid(mat_path)

    _assert_bathy_equal(bathy, loaded_specific)
    _assert_bathy_equal(bathy, loaded_generic)


def test_netcdf_roundtrip_loaders(tmp_path: Path) -> None:
    bathy = _sample_bathy()
    nc_path = tmp_path / "bathy_test.nc"

    save_bathy_grid_netcdf(bathy, nc_path)

    loaded_specific = load_bathy_grid_netcdf(nc_path)
    loaded_generic = load_bathy_grid(nc_path)

    _assert_bathy_equal(bathy, loaded_specific)
    _assert_bathy_equal(bathy, loaded_generic)


def test_load_for_analysis_accepts_grid_result_and_path(tmp_path: Path) -> None:
    bathy = _sample_bathy()
    nc_path = tmp_path / "bathy_test.nc"
    save_bathy_grid_netcdf(bathy, nc_path)

    process_result = BathyProcessResult(
        surveys_all=[],
        surveys_processed=[],
        bathy=bathy,
        x_lims=np.array([0.0, 1.0], dtype=float),
        y_lims=np.array([0.0, 1.0], dtype=float),
        x_min=np.array([0.0, 1.0], dtype=float),
        y_min=np.array([0.0, 1.0], dtype=float),
        output_base="UnitTestInlet",
        min_year_processed=2020,
        max_year_processed=2021,
    )

    from_grid = load_for_analysis(bathy)
    from_result = load_for_analysis(process_result)
    from_path = load_for_analysis(nc_path)

    _assert_bathy_equal(bathy, from_grid)
    _assert_bathy_equal(bathy, from_result)
    _assert_bathy_equal(bathy, from_path)


def test_load_bathy_grid_rejects_unsupported_extension(tmp_path: Path) -> None:
    txt_path = tmp_path / "not_supported.txt"
    txt_path.write_text("placeholder", encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported bathy file extension"):
        load_bathy_grid(txt_path)
