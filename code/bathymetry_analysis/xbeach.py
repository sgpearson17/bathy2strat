"""XBeach model-output adapters for bathy2strat workflows."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from .stratigraphy import BathyCube


def _datetime_to_datenum(value: datetime | str | np.datetime64 | pd.Timestamp) -> float:
    """Convert a datetime-like value to MATLAB datenum."""
    timestamp = pd.Timestamp(value).to_pydatetime()
    midnight = datetime(timestamp.year, timestamp.month, timestamp.day)
    return timestamp.toordinal() + 366 + (timestamp - midnight).total_seconds() / 86400.0


def _require_regular_mesh(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return 1D axes from a regular XBeach mesh."""
    if x.ndim != 2 or y.ndim != 2:
        raise ValueError(f"Expected 2D x/y coordinate arrays. Got {x.shape} and {y.shape}.")
    if x.shape != y.shape:
        raise ValueError(f"x/y coordinate shapes differ: {x.shape} vs {y.shape}.")

    if not np.allclose(x, x[0:1, :], equal_nan=True):
        raise ValueError("XBeach x coordinates must vary by column only for existing transect interpolation.")
    if not np.allclose(y, y[:, 0:1], equal_nan=True):
        raise ValueError("XBeach y coordinates must vary by row only for existing transect interpolation.")

    return np.asarray(x[0, :], dtype=float), np.asarray(y[:, 0], dtype=float)


def load_xbeach_bathy_cube(
    path: str | Path,
    *,
    bed_var: str = "zb",
    x_var: str = "globalx",
    y_var: str = "globaly",
    time_var: str = "globaltime",
    start_time: datetime | str | np.datetime64 | pd.Timestamp = "2020-01-01",
    location: str | None = None,
) -> BathyCube:
    """Load XBeach bed-level output into the package ``BathyCube`` format.

    XBeach stores bed levels as ``(time, y, x)``. The stratigraphy tools expect
    ``(y, x, time)``, with monotonic increasing local x/y coordinates.
    Decreasing axes are flipped together with the bed stack.
    """
    path = Path(path)
    with Dataset(path) as ds:
        missing = [name for name in (bed_var, x_var, y_var, time_var) if name not in ds.variables]
        if missing:
            raise KeyError(f"Missing required XBeach variable(s): {', '.join(missing)}")

        x = np.asarray(ds.variables[x_var][:], dtype=float)
        y = np.asarray(ds.variables[y_var][:], dtype=float)
        t_seconds = np.asarray(ds.variables[time_var][:], dtype=float).reshape(-1)
        z_time_y_x = np.asarray(ds.variables[bed_var][:], dtype=float)

    if z_time_y_x.ndim != 3:
        raise ValueError(f"Expected {bed_var!r} to be 3D (time, y, x). Got {z_time_y_x.shape}.")
    if z_time_y_x.shape[1:] != x.shape:
        raise ValueError(f"{bed_var!r} spatial shape {z_time_y_x.shape[1:]} does not match x/y shape {x.shape}.")
    if z_time_y_x.shape[0] != t_seconds.size:
        raise ValueError(f"{bed_var!r} time length {z_time_y_x.shape[0]} does not match {time_var!r} length {t_seconds.size}.")

    x_axis, y_axis = _require_regular_mesh(x, y)
    z = np.transpose(z_time_y_x, (1, 2, 0))

    if np.all(np.diff(x_axis) < 0):
        x = x[:, ::-1]
        x_axis = x_axis[::-1]
        z = z[:, ::-1, :]
    if np.all(np.diff(y_axis) < 0):
        y = y[::-1, :]
        y_axis = y_axis[::-1]
        z = z[::-1, :, :]

    if not np.all(np.diff(x_axis) > 0):
        raise ValueError("XBeach x axis is not monotonic after loading.")
    if not np.all(np.diff(y_axis) > 0):
        raise ValueError("XBeach y axis is not monotonic after loading.")

    t0 = _datetime_to_datenum(start_time)
    t = t0 + t_seconds / 86400.0
    cube_location = location or path.stem
    return BathyCube(location=cube_location, t=t, x=x, y=y, z=z)


def median_grid_spacing(cube: BathyCube) -> tuple[float, float]:
    """Return median ``(dx, dy)`` spacing from a regular ``BathyCube`` grid."""
    x_axis = np.asarray(cube.x[0, :], dtype=float)
    y_axis = np.asarray(cube.y[:, 0], dtype=float)
    return float(np.nanmedian(np.diff(x_axis))), float(np.nanmedian(np.diff(y_axis)))
