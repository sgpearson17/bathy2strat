"""Tests for XBeach model-output loading helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from bathymetry_analysis.xbeach import load_xbeach_bathy_cube, median_grid_spacing


def _write_xbeach_file(path: Path) -> None:
    """Write a minimal XBeach-style output file."""
    x_axis = np.array([30.0, 20.0, 10.0], dtype=float)
    y_axis = np.array([200.0, 100.0], dtype=float)
    x, y = np.meshgrid(x_axis, y_axis)
    t = np.array([0.0, 3600.0], dtype=float)

    zb = np.empty((2, 2, 3), dtype=float)
    zb[0, :, :] = np.array([[-1.0, -0.8, -0.6], [-1.2, -1.0, -0.9]], dtype=float)
    zb[1, :, :] = zb[0, :, :] + 0.1

    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("globaltime", len(t))
        nc.createDimension("ny", len(y_axis))
        nc.createDimension("nx", len(x_axis))
        nc.createVariable("globaltime", "f8", ("globaltime",))[:] = t
        nc.createVariable("globalx", "f8", ("ny", "nx"))[:, :] = x
        nc.createVariable("globaly", "f8", ("ny", "nx"))[:, :] = y
        nc.createVariable("zb", "f8", ("globaltime", "ny", "nx"))[:, :, :] = zb


def test_load_xbeach_bathy_cube_flips_decreasing_axes(tmp_path: Path) -> None:
    """XBeach output loads as an increasing-axis BathyCube."""
    nc_path = tmp_path / "xboutput.nc"
    _write_xbeach_file(nc_path)

    cube = load_xbeach_bathy_cube(nc_path, start_time="2024-01-01", location="UnitXBeach")

    assert cube.location == "UnitXBeach"
    assert cube.z.shape == (2, 3, 2)
    assert np.all(np.diff(cube.x[0, :]) > 0)
    assert np.all(np.diff(cube.y[:, 0]) > 0)
    np.testing.assert_allclose(cube.z[0, 0, :], [-0.9, -0.8])
    np.testing.assert_allclose(cube.t[1] - cube.t[0], 3600.0 / 86400.0)
    assert median_grid_spacing(cube) == (10.0, 100.0)
