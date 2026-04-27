from __future__ import annotations

import numpy as np

from bathymetry_analysis.stratigraphy import BathyCube, StratigraphyConfig, compute_stratigraphy


def _toy_cube() -> BathyCube:
    x_axis = np.array([1000.0, 1020.0, 1040.0], dtype=float)
    y_axis = np.array([2000.0, 2020.0], dtype=float)
    x, y = np.meshgrid(x_axis, y_axis)

    # Three survey times (MATLAB datenums for simple synthetic test).
    t = np.array([738000.0, 738365.0, 738730.0], dtype=float)

    z = np.empty((2, 3, 3), dtype=float)
    z[:, :, 0] = np.array([[-1.0, -0.9, -0.8], [-1.2, -1.1, -1.0]], dtype=float)
    z[:, :, 1] = np.array([[-0.8, -1.1, -0.6], [-1.0, -1.3, -0.9]], dtype=float)
    z[:, :, 2] = np.array([[-0.7, -1.0, -0.5], [-0.9, -1.25, -0.8]], dtype=float)

    return BathyCube(location="UnitTest", t=t, x=x, y=y, z=z)


def test_compute_stratigraphy_shapes_and_finiteness() -> None:
    cube = _toy_cube()
    cfg = StratigraphyConfig(initial_index=0, dx=20.0)

    result = compute_stratigraphy(cube, cfg)

    assert result.min_surf.shape == cube.x.shape
    assert result.max_surf.shape == cube.x.shape
    assert result.deposit_elev.shape == cube.z.shape
    assert result.deposit_per_year.shape == (len(cube.t), len(cube.t))
    assert result.theseus_ratio.shape == (len(cube.t), len(cube.t))

    assert result.deposit_thk_full.ndim == 3
    assert result.deposit_thk_full.shape[0:2] == cube.x.shape
    assert result.deposit_thk_full.shape[2] >= 12

    assert np.isfinite(result.total_sed_vol_per_year).all()
    assert np.isfinite(result.deposit_per_year[0, :]).all()


def test_compute_stratigraphy_handles_persistent_nan_cells() -> None:
    cube = _toy_cube()
    cube.z[0, 0, 1] = np.nan

    result = compute_stratigraphy(cube, StratigraphyConfig(initial_index=0, dx=20.0))

    # One NaN in one survey should propagate to all surveys at that cell.
    assert np.isnan(result.deposit_elev[0, 0, :]).all()
