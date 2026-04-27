"""
bathy_formatter.py

Convert bathymetry netCDF files into MATLAB-compatible .mat outputs,
create QC plots, and interpolate all surveys onto a common grid.

First-pass Python conversion of MATLAB bathyFormatter_v004.m.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
from netCDF4 import Dataset
from scipy.io import loadmat, savemat
from scipy.interpolate import CloughTocher2DInterpolator, LinearNDInterpolator, NearestNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import binary_erosion
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from scipy.spatial import KDTree


MHW = 0.358
MSL = -0.112
MLW = -0.590


KG2_SGP_CLR = [
    (0.0000, (0, 67, 143)),
    (0.3333, (13, 182, 255)),
    (0.5000, (255, 255, 255)),
    (0.6000, (199, 181, 181)),
    (0.7333, (158, 144, 144)),
    (0.7667, (29, 89, 74)),
    (1.0000, (29, 89, 74)),
]

VINTAGE_SGP_CLR = [
    (0.0000, (27, 126, 129)),
    (0.1333, (41, 155, 151)),
    (0.2667, (56, 170, 164)),
    (0.4000, (130, 199, 180)),
    (0.5000, (220, 231, 194)),
    (0.5667, (255, 240, 196)),
    (0.6267, (244, 214, 176)),
    (0.6667, (217, 188, 146)),
    (0.7133, (255, 221, 146)),
    (1.0000, (226, 129, 61)),
]


@dataclass
class RawSurvey:
    location: str
    datenum: float
    vertical_datum: str
    horizontal_datum: str | None
    ncfile: str
    x_raw: np.ndarray
    y_raw: np.ndarray
    z_raw: np.ndarray


@dataclass
class BathyGrid:
    location: str
    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray


@dataclass
class BathyProcessResult:
    surveys_all: list[RawSurvey]
    surveys_processed: list[RawSurvey]
    bathy: BathyGrid
    x_lims: np.ndarray
    y_lims: np.ndarray
    x_min: np.ndarray
    y_min: np.ndarray
    output_base: str
    min_year_processed: int
    max_year_processed: int


def _matlab_datenum(year: int, month: int, day: int = 1) -> float:
    """Convert a calendar date to MATLAB-style datenum."""
    # MATLAB day 1 == 0000-01-01, Python ordinal day 1 == 0001-01-01
    # Offset = 366 days.
    from datetime import datetime

    # Force invalid/zero months to January for cleaner file naming.
    # Example: datenum(2005,0,1) -> 2005-01-01.
    if month < 1 or month > 12:
        month = 1

    dt = datetime(year, month, day)
    return dt.toordinal() + 366

def _datenum_to_year_month(datenum_value: float) -> tuple[int, int]:
    from datetime import datetime, timedelta

    dt = datetime.fromordinal(int(datenum_value)) + timedelta(days=datenum_value % 1) - timedelta(days=366)
    return dt.year, dt.month


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _clr_stops_to_cmap(stops: list[tuple[float, tuple[int, int, int]]], name: str) -> LinearSegmentedColormap:
    color_list: list[tuple[float, tuple[float, float, float]]] = []
    for position, rgb in stops:
        color_list.append((position, (rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0)))
    return LinearSegmentedColormap.from_list(name, color_list)


def load_clrmap_file(clrmap_path: str | Path) -> LinearSegmentedColormap:
    """Load a MATLAB-style .clrmap file into a matplotlib colormap."""
    path = Path(clrmap_path)
    if not path.exists():
        raise FileNotFoundError(f"Colormap file not found: {path}")

    stops: list[tuple[float, tuple[int, int, int]]] = []
    cmap_name = path.stem

    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("NAME="):
                cmap_name = line.split("=", 1)[1].strip() or cmap_name
                continue
            if line.startswith("COLORMAP") or line.startswith("SPACE="):
                continue

            parts = line.split()
            if len(parts) != 4:
                continue

            pos = float(parts[0])
            r = int(parts[1])
            g = int(parts[2])
            b = int(parts[3])
            stops.append((pos, (r, g, b)))

    if not stops:
        raise ValueError(f"No color stops found in: {path}")

    return _clr_stops_to_cmap(stops, cmap_name)


def get_named_colormap(name: str, cmap_file: str | Path | None = None):
    """Get plotting colormap by name or custom .clrmap file."""
    if cmap_file is not None:
        return load_clrmap_file(cmap_file)

    key = name.lower()
    if key == "kg2":
        return _clr_stops_to_cmap(KG2_SGP_CLR, "kg2_sgp")
    if key == "vintage":
        return _clr_stops_to_cmap(VINTAGE_SGP_CLR, "vintage_sgp")
    if key == "turbo":
        return plt.get_cmap("turbo")
    if key == "viridis":
        return plt.get_cmap("viridis")
    raise ValueError("Unsupported colormap name. Choose from: kg2, vintage, turbo, viridis")


def load_raw_surveys(nc_dir: str | Path) -> list[RawSurvey]:
    """Load all .nc bathymetry surveys from a directory."""
    nc_path = Path(nc_dir)
    nc_files = sorted(nc_path.glob("*.nc"))
    if not nc_files:
        raise FileNotFoundError(f"No .nc files found in: {nc_path}")

    surveys: list[RawSurvey] = []

    for file_path in nc_files:
        with Dataset(file_path) as ds:
            x_temp = np.array(ds.variables["x"][:], dtype=float)
            y_temp = np.array(ds.variables["y"][:], dtype=float)
            z_temp = np.array(ds.variables["z"][:], dtype=float)

            missing_value = None
            z_var = ds.variables["z"]
            if hasattr(z_var, "missing_value"):
                missing_value = float(getattr(z_var, "missing_value"))
            elif hasattr(z_var, "_FillValue"):
                missing_value = float(getattr(z_var, "_FillValue"))

            if missing_value is not None:
                z_temp[z_temp == missing_value] = np.nan

            # US Survey Feet to meters: 1200/3937.
            z_raw = z_temp * 1200.0 / 3937.0
            x_raw = np.repeat(x_temp[:, None], len(y_temp), axis=1)
            y_raw = np.repeat(y_temp[None, :], len(x_temp), axis=0)

            # Align z orientation to match x_raw/y_raw. Many netCDFs store z as
            # [y, x] while MATLAB workflows often operate on [x, y].
            if z_raw.shape != x_raw.shape:
                if z_raw.T.shape == x_raw.shape:
                    z_raw = z_raw.T
                else:
                    raise ValueError(
                        f"Grid shape mismatch for {file_path.name}: "
                        f"x/y shape {x_raw.shape}, z shape {z_raw.shape}"
                    )

            name_parts = file_path.stem.split("_")
            if len(name_parts) < 4:
                raise ValueError(
                    f"Filename format not recognized for metadata parse: {file_path.name}. "
                    "Expected e.g. Location_YYYY_MM_Datum.nc"
                )

            location = name_parts[0]
            year = int(name_parts[1])
            month = int(name_parts[2])
            vertical_datum = name_parts[3]
            datenum_value = _matlab_datenum(year, month, 1)

            horizontal_datum = None
            if hasattr(z_var, "esri_pe_string"):
                horizontal_datum = str(getattr(z_var, "esri_pe_string"))

        surveys.append(
            RawSurvey(
                location=location,
                datenum=datenum_value,
                vertical_datum=vertical_datum,
                horizontal_datum=horizontal_datum,
                ncfile=file_path.name,
                x_raw=x_raw,
                y_raw=y_raw,
                z_raw=z_raw,
            )
        )

    surveys.sort(key=lambda s: s.datenum)
    return surveys


def save_raw_surveys_mat(surveys: list[RawSurvey], out_path: str | Path) -> None:
    """Save raw survey list to MATLAB .mat file."""
    recs: list[dict[str, Any]] = []
    for s in surveys:
        recs.append(
            {
                "location": s.location,
                "datenum": s.datenum,
                "verticalDatum": s.vertical_datum,
                "horizontalDatum": s.horizontal_datum or "",
                "ncfile": s.ncfile,
                "xRaw": s.x_raw,
                "yRaw": s.y_raw,
                "zRaw": s.z_raw,
            }
        )

    savemat(str(out_path), {"bathyRaw": np.array(recs, dtype=object)})


def remove_year(surveys: list[RawSurvey], year_to_remove: int) -> list[RawSurvey]:
    """Return surveys excluding any that match the target year."""
    filtered: list[RawSurvey] = []
    for s in surveys:
        y, _ = _datenum_to_year_month(s.datenum)
        if y != year_to_remove:
            filtered.append(s)
    return filtered


def compute_domain_extents(surveys: list[RawSurvey]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute max extents and overlap extents in km."""
    if not surveys:
        raise ValueError("No surveys provided for extent computation")

    first = surveys[0]
    x_vec = first.x_raw.reshape(-1)
    y_vec = first.y_raw.reshape(-1)
    z_vec = first.z_raw.reshape(-1)
    mask = ~np.isnan(z_vec)
    xv = x_vec[mask] / 1000.0
    yv = y_vec[mask] / 1000.0

    x_lims = np.array([np.min(xv), np.max(xv)], dtype=float)
    y_lims = np.array([np.min(yv), np.max(yv)], dtype=float)
    x_min = x_lims.copy()
    y_min = y_lims.copy()

    for s in surveys[1:]:
        x_vec = s.x_raw.reshape(-1)
        y_vec = s.y_raw.reshape(-1)
        z_vec = s.z_raw.reshape(-1)
        mask = ~np.isnan(z_vec)

        xv = x_vec[mask] / 1000.0
        yv = y_vec[mask] / 1000.0

        xl = np.array([np.min(xv), np.max(xv)])
        yl = np.array([np.min(yv), np.max(yv)])

        x_lims[0] = min(x_lims[0], xl[0])
        x_lims[1] = max(x_lims[1], xl[1])
        y_lims[0] = min(y_lims[0], yl[0])
        y_lims[1] = max(y_lims[1], yl[1])

        x_min[0] = max(x_min[0], xl[0])
        x_min[1] = min(x_min[1], xl[1])
        y_min[0] = max(y_min[0], yl[0])
        y_min[1] = min(y_min[1], yl[1])

    return x_lims, y_lims, x_min, y_min


def _plot_extent_outline(ax, survey: RawSurvey, color, method: str = "mask") -> None:
    """Plot one survey outline for the extents figure.

    Parameters
    ----------
    method : str
        "mask" traces the valid-data boundary (tighter fit).
        "convex" uses a convex hull (broader fit).
    """
    method_key = method.lower()

    if method_key == "mask":
        valid = np.isfinite(survey.z_raw).astype(float)
        if np.count_nonzero(valid) < 3:
            return
        # Contour of validity mask provides a tighter boundary than convex hull.
        ax.contour(
            survey.x_raw / 1000.0,
            survey.y_raw / 1000.0,
            valid,
            levels=[0.5],
            colors=[color],
            linewidths=1.2,
        )
        return

    if method_key == "convex":
        x_vec = survey.x_raw.reshape(-1)
        y_vec = survey.y_raw.reshape(-1)
        z_vec = survey.z_raw.reshape(-1)
        mask = ~np.isnan(z_vec)
        pts = np.column_stack((x_vec[mask] / 1000.0, y_vec[mask] / 1000.0))

        if len(pts) >= 3:
            hull = ConvexHull(pts)
            loop = np.append(hull.vertices, hull.vertices[0])
            ax.plot(pts[loop, 0], pts[loop, 1], color=color, linewidth=1.2)
        return

    raise ValueError("extent_boundary_method must be one of: mask, convex")


def plot_extents(
    surveys: list[RawSurvey],
    plot_path: str | Path,
    extents_cmap,
    extent_boundary_method: str = "mask",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Plot spatial extents of all raw surveys."""
    x_lims, y_lims, x_min, y_min = compute_domain_extents(surveys)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor("white")

    colors = extents_cmap(np.linspace(0, 1, len(surveys)))

    for i, s in enumerate(surveys):
        _plot_extent_outline(ax, s, colors[i], method=extent_boundary_method)

    min_year, _ = _datenum_to_year_month(min(s.datenum for s in surveys))
    max_year, _ = _datenum_to_year_month(max(s.datenum for s in surveys))
    ax.set_title(f"{surveys[0].location} Bathymetry Extents {min_year}-{max_year}")
    ax.set_xlabel("Easting [km UTM18]")
    ax.set_ylabel("Northing [km UTM18]")
    ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim((float(x_lims[0]), float(x_lims[1])))
    ax.set_ylim((float(y_lims[0]), float(y_lims[1])))

    sm = plt.cm.ScalarMappable(cmap=extents_cmap, norm=Normalize(min_year, max_year))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Time [year]")

    fig.tight_layout()
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return x_lims, y_lims, x_min, y_min


def plot_raw_surveys(
    surveys: list[RawSurvey],
    x_lims: np.ndarray,
    y_lims: np.ndarray,
    out_dir: str | Path,
    bathy_cmap,
    mask_nan: bool = True,
) -> None:
    """Plot each raw bathymetry survey.
    
    Parameters
    ----------
    mask_nan : bool
        If True, mask out NaN values so they appear white instead of filled.
    """
    out_path = Path(out_dir)
    _ensure_dir(out_path)

    for s in surveys:
        fig, ax = plt.subplots(figsize=(12, 8))
        levels = np.arange(-30, 5.2, 0.2)
        z_plot = np.ma.masked_where(np.isnan(s.z_raw), s.z_raw) if mask_nan else s.z_raw
        cf = ax.contourf(s.x_raw / 1000.0, s.y_raw / 1000.0, z_plot, levels=levels, cmap=bathy_cmap)
        ax.contour(s.x_raw / 1000.0, s.y_raw / 1000.0, s.z_raw, levels=[MLW], colors=[(0.5, 0.5, 0.5)], linewidths=1.0)
        ax.contour(s.x_raw / 1000.0, s.y_raw / 1000.0, s.z_raw, levels=[-6], colors="k", linestyles=":", linewidths=0.5)

        y, m = _datenum_to_year_month(s.datenum)
        ax.set_title(f"{s.location} {y:04d}-{m:02d}")
        ax.set_xlabel("Easting [km UTM18]")
        ax.set_ylabel("Northing [km UTM18]")
        ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim((float(x_lims[0]), float(x_lims[1])))
        ax.set_ylim((float(y_lims[0]), float(y_lims[1])))
        cf.set_clim(-10, 5)

        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label("Elevation [m NAVD88]")

        fig.tight_layout()
        fig.savefig(out_path / f"{s.location}_RawBathy_{y:04d}_{m:02d}_py.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def regrid_surveys(
    surveys: list[RawSurvey],
    dx: float,
    x_min: np.ndarray,
    y_min: np.ndarray,
    interp_method: str = "linear",
    boundary_tightness: float = 3.0,
    overlap_erosion_cells: int = 0,
) -> BathyGrid:
    """Interpolate all raw surveys onto a common overlapping grid."""
    x_grid_lim = x_min * 1000.0
    y_grid_lim = y_min * 1000.0

    x_grid_lim[0] = dx * np.floor(x_grid_lim[0] / dx)
    x_grid_lim[1] = dx * np.ceil(x_grid_lim[1] / dx)
    y_grid_lim[0] = dx * np.floor(y_grid_lim[0] / dx)
    y_grid_lim[1] = dx * np.ceil(y_grid_lim[1] / dx)

    gx, gy = np.meshgrid(np.arange(x_grid_lim[0], x_grid_lim[1] + dx, dx), np.arange(y_grid_lim[0], y_grid_lim[1] + dx, dx))

    nt = len(surveys)
    gz = np.full((gx.shape[0], gx.shape[1], nt), np.nan)
    t = np.zeros((nt, 1), dtype=float)

    for i, s in enumerate(surveys):
        t[i, 0] = s.datenum

        x = s.x_raw.reshape(-1)
        y = s.y_raw.reshape(-1)
        z = s.z_raw.reshape(-1)
        mask = ~np.isnan(z)

        points = np.column_stack((x[mask], y[mask]))
        method = interp_method.lower()
        if method == "linear":
            interp = LinearNDInterpolator(points, z[mask], fill_value=np.nan)
        elif method == "nearest":
            interp = NearestNDInterpolator(points, z[mask])
        elif method == "cubic":
            interp = CloughTocher2DInterpolator(points, z[mask], fill_value=np.nan)
        else:
            raise ValueError("interp_method must be one of: linear, nearest, cubic")
        gz[:, :, i] = interp(gx, gy)

    # Remove non-overlapping cells by enforcing NaN where any survey has NaN.
    nan_mask = np.any(np.isnan(gz), axis=2)
    gz[nan_mask, :] = np.nan

    # Enforce minimum-overlap footprint from raw-survey validity masks so the
    # regridded domain matches the extents-mask concept (not convex hull fill).
    overlap_mask = compute_minimum_overlap_mask(
        surveys,
        gx,
        gy,
        boundary_tightness=boundary_tightness,
        overlap_erosion_cells=overlap_erosion_cells,
    )
    gz[~overlap_mask, :] = np.nan

    return BathyGrid(
        location=surveys[0].location,
        t=t,
        x=gx,
        y=gy,
        z=gz,
    )


def compute_minimum_overlap_mask(
    raw_surveys: list[RawSurvey],
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    boundary_tightness: float = 3.0,
    overlap_erosion_cells: int = 0,
) -> np.ndarray:
    """Compute a global overlap mask from all raw-survey validity footprints.

    Returns
    -------
    np.ndarray
        Boolean array with True where all surveys have valid-data footprint.
    """
    if not raw_surveys:
        raise ValueError("No raw surveys provided for overlap-mask computation")

    overlap_mask = np.ones_like(grid_x, dtype=bool)

    for raw in raw_surveys:
        survey_mask = _compute_single_survey_boundary_mask(
            raw,
            grid_x,
            grid_y,
            boundary_tightness=boundary_tightness,
        )
        overlap_mask &= survey_mask

    pre_erosion_mask = overlap_mask.copy()

    if overlap_erosion_cells > 0:
        overlap_mask = np.asarray(
            binary_erosion(
                overlap_mask,
                structure=np.ones((3, 3), dtype=bool),
                iterations=overlap_erosion_cells,
                border_value=0,
            ),
            dtype=bool,
        )

    # Safety: never return an empty domain. If erosion removes everything,
    # fall back to pre-erosion overlap; if still empty, use nearest-mask overlap.
    if not np.any(overlap_mask):
        overlap_mask = pre_erosion_mask

    if not np.any(overlap_mask):
        nearest_overlap = np.ones_like(grid_x, dtype=bool)
        for raw in raw_surveys:
            valid_raw = np.isfinite(raw.z_raw).astype(float)
            x_axis = raw.x_raw[:, 0]
            y_axis = raw.y_raw[0, :]
            valid_interp = RegularGridInterpolator(
                (y_axis, x_axis),
                valid_raw.T,
                method="nearest",
                bounds_error=False,
                fill_value=0.0,
            )
            valid_on_grid = valid_interp(np.column_stack((grid_y.ravel(), grid_x.ravel()))).reshape(grid_x.shape)
            nearest_overlap &= valid_on_grid >= 0.5
        overlap_mask = nearest_overlap

    return overlap_mask


def _compute_single_survey_boundary_mask(
    raw: RawSurvey,
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    boundary_tightness: float,
) -> np.ndarray:
    """Compute a boundary-like footprint mask for one survey.

    This emulates MATLAB's `boundary(x,y)` behavior more closely than a convex hull
    by using an alpha-shape-like filter on Delaunay triangles.
    """
    valid = np.isfinite(raw.z_raw)
    if np.count_nonzero(valid) < 3:
        return np.zeros_like(grid_x, dtype=bool)

    # Use edge points to reduce point count and focus the boundary geometry.
    eroded = np.asarray(binary_erosion(valid, structure=np.ones((3, 3), dtype=bool), border_value=0), dtype=bool)
    edge = valid & (~eroded)
    if np.count_nonzero(edge) < 3:
        edge = valid

    points = np.column_stack((raw.x_raw[edge], raw.y_raw[edge]))
    if points.shape[0] < 3:
        return np.zeros_like(grid_x, dtype=bool)

    # Remove duplicates to keep triangulation stable.
    points = np.unique(points, axis=0)
    if points.shape[0] < 3:
        return np.zeros_like(grid_x, dtype=bool)

    if points.shape[0] < 4:
        # Too few points for a robust concave boundary; nearest mask fallback.
        valid_raw = valid.astype(float)
        x_axis = raw.x_raw[:, 0]
        y_axis = raw.y_raw[0, :]
        valid_interp = RegularGridInterpolator(
            (y_axis, x_axis),
            valid_raw.T,
            method="nearest",
            bounds_error=False,
            fill_value=0.0,
        )
        valid_on_grid = valid_interp(np.column_stack((grid_y.ravel(), grid_x.ravel()))).reshape(grid_x.shape)
        return valid_on_grid >= 0.5

    # Derive an adaptive circumradius threshold from local point spacing,
    # similar in spirit to MATLAB boundary shrink behavior.
    tree = KDTree(points)
    nn = tree.query(points, k=2)[0][:, 1]
    nn = nn[np.isfinite(nn)]
    if nn.size == 0:
        valid_raw = valid.astype(float)
        x_axis = raw.x_raw[:, 0]
        y_axis = raw.y_raw[0, :]
        valid_interp = RegularGridInterpolator(
            (y_axis, x_axis),
            valid_raw.T,
            method="nearest",
            bounds_error=False,
            fill_value=0.0,
        )
        valid_on_grid = valid_interp(np.column_stack((grid_y.ravel(), grid_x.ravel()))).reshape(grid_x.shape)
        return valid_on_grid >= 0.5

    median_spacing = float(np.median(nn))
    max_circumradius = max(boundary_tightness, 0.5) * median_spacing

    tri = Delaunay(points)
    triangles = tri.simplices
    p1 = points[triangles[:, 0]]
    p2 = points[triangles[:, 1]]
    p3 = points[triangles[:, 2]]

    a = np.linalg.norm(p2 - p3, axis=1)
    b = np.linalg.norm(p1 - p3, axis=1)
    c = np.linalg.norm(p1 - p2, axis=1)
    s = 0.5 * (a + b + c)
    area_sq = s * (s - a) * (s - b) * (s - c)
    area_sq = np.maximum(area_sq, 0.0)
    area = np.sqrt(area_sq)

    circumradius = np.full_like(area, np.inf)
    nondeg = area > 0
    circumradius[nondeg] = (a[nondeg] * b[nondeg] * c[nondeg]) / (4.0 * area[nondeg])

    keep = np.isfinite(circumradius) & (circumradius <= max_circumradius)
    if not np.any(keep):
        # Fallback if threshold is too strict.
        keep = np.isfinite(circumradius)

    triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles=triangles[keep])
    trifinder = triang.get_trifinder()
    tri_ids = np.asarray(trifinder(grid_x.ravel(), grid_y.ravel()))
    inside = tri_ids != -1
    return np.asarray(inside, dtype=bool).reshape(grid_x.shape)


def save_bathy_grid_mat(bathy: BathyGrid, out_path: str | Path) -> None:
    """Save regridded bathymetry in MATLAB-struct style."""
    savemat(
        str(out_path),
        {
            "bathy": {
                "location": bathy.location,
                "t": bathy.t,
                "x": bathy.x,
                "y": bathy.y,
                "z": bathy.z,
            }
        },
    )


def _coerce_mat_struct(value: Any) -> Any:
    """Extract a scipy.io-loaded MATLAB struct from possible ndarray wrappers."""
    if isinstance(value, np.ndarray) and value.dtype == object:
        return value.flat[0]
    return value


def _to_time_column_vector(t_values: np.ndarray) -> np.ndarray:
    """Normalize time values to an (nt, 1) MATLAB-like column vector."""
    t_flat = np.asarray(t_values, dtype=float).reshape(-1)
    return t_flat.reshape(-1, 1)


def load_bathy_grid_mat(in_path: str | Path) -> BathyGrid:
    """Load regridded bathymetry from MATLAB .mat output produced by this workflow."""
    mat = loadmat(str(in_path), squeeze_me=True, struct_as_record=False)
    if "bathy" not in mat:
        raise ValueError(f"MAT file does not contain 'bathy' struct: {in_path}")

    bathy_struct = _coerce_mat_struct(mat["bathy"])
    location = str(getattr(bathy_struct, "location"))
    t = _to_time_column_vector(np.asarray(getattr(bathy_struct, "t"), dtype=float))
    x = np.asarray(getattr(bathy_struct, "x"), dtype=float)
    y = np.asarray(getattr(bathy_struct, "y"), dtype=float)
    z = np.asarray(getattr(bathy_struct, "z"), dtype=float)

    return BathyGrid(location=location, t=t, x=x, y=y, z=z)


def load_bathy_grid_netcdf(in_path: str | Path) -> BathyGrid:
    """Load regridded bathymetry from CF-style netCDF output produced by this workflow."""
    with Dataset(str(in_path), "r") as nc:
        if "bathymetry" not in nc.variables:
            raise ValueError(f"NetCDF file missing 'bathymetry' variable: {in_path}")
        if "time" not in nc.variables or "x" not in nc.variables or "y" not in nc.variables:
            raise ValueError(f"NetCDF file missing one of required coordinates: time, x, y in {in_path}")

        z_tyx = np.asarray(nc.variables["bathymetry"][:], dtype=float)
        x_1d = np.asarray(nc.variables["x"][:], dtype=float)
        y_1d = np.asarray(nc.variables["y"][:], dtype=float)
        t_cf = np.asarray(nc.variables["time"][:], dtype=float).reshape(-1)

        # Convert CF time (days since 1970-01-01) back to MATLAB datenum.
        cf_epoch_datenum = 719163.0
        t_matlab = _to_time_column_vector(t_cf + cf_epoch_datenum)

        x_grid, y_grid = np.meshgrid(x_1d, y_1d)
        z_yxt = np.moveaxis(z_tyx, 0, -1)

        location = str(getattr(nc, "location", Path(in_path).stem))

    return BathyGrid(location=location, t=t_matlab, x=x_grid, y=y_grid, z=z_yxt)


def load_bathy_grid(in_path: str | Path) -> BathyGrid:
    """Load regridded bathymetry from .mat or .nc file path."""
    path = Path(in_path)
    suffix = path.suffix.lower()
    if suffix == ".mat":
        return load_bathy_grid_mat(path)
    if suffix == ".nc":
        return load_bathy_grid_netcdf(path)
    raise ValueError(f"Unsupported bathy file extension: {path.suffix}. Use .mat or .nc")


def load_for_analysis(source: BathyGrid | BathyProcessResult | str | Path) -> BathyGrid:
    """Return a BathyGrid from either an in-memory object or a saved .mat/.nc file."""
    if isinstance(source, BathyGrid):
        return source
    if isinstance(source, BathyProcessResult):
        return source.bathy
    # Notebook-friendly fallback for module reloads where class identity changes.
    if hasattr(source, "bathy"):
        bathy_obj = getattr(source, "bathy")
        if hasattr(bathy_obj, "x") and hasattr(bathy_obj, "y") and hasattr(bathy_obj, "z") and hasattr(bathy_obj, "t"):
            return bathy_obj
    if hasattr(source, "x") and hasattr(source, "y") and hasattr(source, "z") and hasattr(source, "t"):
        return source
    return load_bathy_grid(source)


def _as_bathy_process_result(process_result: BathyProcessResult | str | Path) -> tuple[BathyProcessResult, bool]:
    """Return a BathyProcessResult and whether it was loaded from file path."""
    if isinstance(process_result, BathyProcessResult):
        return process_result, False

    bathy = load_bathy_grid(process_result)
    # Minimal shell result for plotting from persisted data only.
    dummy = BathyProcessResult(
        surveys_all=[],
        surveys_processed=[],
        bathy=bathy,
        x_lims=np.array([np.nan, np.nan], dtype=float),
        y_lims=np.array([np.nan, np.nan], dtype=float),
        x_min=np.array([np.nan, np.nan], dtype=float),
        y_min=np.array([np.nan, np.nan], dtype=float),
        output_base=bathy.location,
        min_year_processed=0,
        max_year_processed=0,
    )
    return dummy, True


def save_bathy_grid_netcdf(bathy: BathyGrid, out_path: str | Path) -> None:
    """Save regridded bathymetry as CF-compliant netCDF file.
    
    Creates a netCDF file with:
    - Dimensions: time, y, x
    - Coordinate variables: time (days since 1970-01-01), y, x (meters UTM)
    - Data variable: bathymetry (time, y, x) with bathymetric depth
    - Grid mapping: UTM 18N projection
    - Standard CF metadata attributes
    """
    from datetime import datetime
    
    out_path = Path(out_path)
    
    # Ensure bathy.t is a 1D numpy array
    time_vals = np.atleast_1d(np.asarray(bathy.t).flatten())
    
    # Extract 1D coordinates from bathy.x and bathy.y
    # If they're 2D (meshgrid format), extract the unique values; else use as-is
    if bathy.x.ndim == 2:
        x_1d = bathy.x[0, :]  # Extract first row (x doesn't vary in y direction)
    else:
        x_1d = bathy.x
    
    if bathy.y.ndim == 2:
        y_1d = bathy.y[:, 0]  # Extract first column (y doesn't vary in x direction)
    else:
        y_1d = bathy.y
    
    # Convert MATLAB datenum to CF time (days since 1970-01-01)
    # MATLAB datenum epoch: 0000-01-01 (year 0)
    # CF epoch: 1970-01-01
    # Offset: datenum(1970,1,1) = 719163
    cf_epoch_datenum = 719163.0
    time_cf = time_vals - cf_epoch_datenum
    
    # Create netCDF dataset
    with Dataset(str(out_path), 'w', format='NETCDF4') as nc:
        # Create dimensions
        nt = len(time_vals)
        ny = len(y_1d)
        nx = len(x_1d)
        nc.createDimension('time', nt)
        nc.createDimension('y', ny)
        nc.createDimension('x', nx)
        
        # Create coordinate variables
        time_var = nc.createVariable('time', 'f8', ('time',))
        time_var[:] = time_cf
        time_var.units = 'days since 1970-01-01 00:00:00'
        time_var.calendar = 'standard'
        time_var.long_name = 'time'
        time_var.standard_name = 'time'
        
        y_var = nc.createVariable('y', 'f8', ('y',))
        y_var[:] = y_1d
        y_var.units = 'meters'
        y_var.long_name = 'Northing coordinate in UTM zone 18N'
        y_var.standard_name = 'projection_y_coordinate'
        y_var.axis = 'Y'
        
        x_var = nc.createVariable('x', 'f8', ('x',))
        x_var[:] = x_1d
        x_var.units = 'meters'
        x_var.long_name = 'Easting coordinate in UTM zone 18N'
        x_var.standard_name = 'projection_x_coordinate'
        x_var.axis = 'X'
        
        # Create bathymetry data variable
        z_var = nc.createVariable('bathymetry', 'f4', ('time', 'y', 'x'),
                                   fill_value=np.nan)
        # Note: bathy.z is stored as (y, x, time), but CF netCDF expects (time, y, x)
        z_var[:] = np.moveaxis(bathy.z, -1, 0)  # Move time axis from -1 to 0: (y,x,t) -> (t,y,x)
        z_var.units = 'meters'
        z_var.long_name = 'bathymetric depth (positive down from NAVD88)'
        z_var.standard_name = 'sea_floor_depth'
        z_var.grid_mapping = 'crs'
        z_var.coordinates = 'time y x'
        
        # Create grid mapping variable for UTM 18N
        crs_var = nc.createVariable('crs', 'i4')
        crs_var.grid_mapping_name = 'transverse_mercator'
        crs_var.longitude_of_central_meridian = -75.0
        crs_var.latitude_of_projection_origin = 0.0
        crs_var.scale_factor_at_central_meridian = 0.9996
        crs_var.false_easting = 500000.0
        crs_var.false_northing = 0.0
        crs_var.spatial_ref = 'PROJCS["WGS 84 / UTM zone 18N",GEOGCS["WGS 84",' \
                              'DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],' \
                              'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],' \
                              'PROJECTION["Transverse_Mercator"],' \
                              'PARAMETER["latitude_of_origin",0],' \
                              'PARAMETER["central_meridian",-75],' \
                              'PARAMETER["scale_factor",0.9996],' \
                              'PARAMETER["false_easting",500000],' \
                              'PARAMETER["false_northing",0],' \
                              'UNIT["metre",1]]'
        crs_var.GeoTransform = ''
        
        # Global attributes
        nc.title = f'Regridded Bathymetry: {bathy.location}'
        nc.institution = 'University of North Carolina'
        nc.source = 'bathyFormatter.py - MATLAB to Python conversion'
        nc.history = f'Created {datetime.now().isoformat()}'
        nc.Conventions = 'CF-1.8'
        nc.references = 'https://www.unidata.ucar.edu/software/netcdf/'
        nc.location = bathy.location
        nc.vertical_datum = 'NAVD88'
        nc.horizontal_datum = 'WGS84 / UTM zone 18N'


def plot_regridded_surveys(
    bathy: BathyGrid,
    raw_surveys: list[RawSurvey] | None,
    out_dir: str | Path,
    bathy_cmap,
    extent_boundary_method: str = "mask",
    overlap_mask: np.ndarray | None = None,
    mask_nan: bool = True,
) -> None:
    """Plot each regridded survey.
    
    Parameters
    ----------
    mask_nan : bool
        If True, mask out NaN values so they appear white instead of filled.
    """
    out_path = Path(out_dir)
    _ensure_dir(out_path)

    method_key = extent_boundary_method.lower()

    t_vals = np.asarray(bathy.t, dtype=float).reshape(-1)
    if bathy.z.shape[2] != t_vals.size:
        raise ValueError("Bathy time dimension mismatch: expected z.shape[2] == len(t)")

    if raw_surveys is None:
        raw_seq: list[RawSurvey] = []
    else:
        raw_seq = raw_surveys
    have_raw = bool(raw_seq) and len(raw_seq) == t_vals.size

    for i in range(t_vals.size):
        fig, ax = plt.subplots(figsize=(11.4, 7.9))

        levels = np.arange(-30, 5.2, 0.2)
        z_data = bathy.z[:, :, i]

        if method_key == "mask":
            if overlap_mask is not None:
                outside_mask = ~overlap_mask
            elif not have_raw:
                outside_mask = np.zeros_like(z_data, dtype=bool)
            else:
                # Backward-compatible fallback: use per-survey footprint.
                raw = raw_seq[i]
                valid_raw = np.isfinite(raw.z_raw).astype(float)
                x_axis = raw.x_raw[:, 0]
                y_axis = raw.y_raw[0, :]
                valid_interp = RegularGridInterpolator(
                    (y_axis, x_axis),
                    valid_raw.T,
                    method="nearest",
                    bounds_error=False,
                    fill_value=0.0,
                )
                valid_on_grid = valid_interp(np.column_stack((bathy.y.ravel(), bathy.x.ravel()))).reshape(bathy.x.shape)
                outside_mask = valid_on_grid < 0.5
        elif method_key == "convex":
            outside_mask = np.zeros_like(z_data, dtype=bool)
        else:
            raise ValueError("extent_boundary_method must be one of: mask, convex")

        if mask_nan:
            z_plot = np.ma.masked_where(np.isnan(z_data) | outside_mask, z_data)
        else:
            z_plot = np.ma.masked_where(outside_mask, z_data)

        cf = ax.contourf(bathy.x / 1000.0, bathy.y / 1000.0, z_plot, levels=levels, cmap=bathy_cmap)
        z_contour = np.ma.masked_where(outside_mask | np.isnan(z_data), z_data)
        ax.contour(bathy.x / 1000.0, bathy.y / 1000.0, z_contour, levels=[MLW], colors=[(0.5, 0.5, 0.5)], linewidths=1.0)
        ax.contour(bathy.x / 1000.0, bathy.y / 1000.0, z_contour, levels=[-6], colors="k", linestyles=":", linewidths=0.5)

        y, m = _datenum_to_year_month(float(t_vals[i]))
        ax.text(307.9, 3835.9, f"{y:04d}-{m:02d}", fontsize=16, weight="bold", style="italic")

        ax.set_xlim((307.0, 309.0))
        ax.set_ylim((3834.3, 3836.1))
        ax.set_xticks(np.arange(307, 309.01, 0.2))
        ax.set_yticks(np.arange(3834.3, 3836.11, 0.2))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
        cf.set_clim(-10, 5)

        fig.tight_layout()
        if have_raw:
            y_raw, m_raw = _datenum_to_year_month(raw_seq[i].datenum)
            location = raw_seq[i].location
        else:
            y_raw, m_raw = y, m
            location = bathy.location
        fig.savefig(out_path / f"{location}_ReGridBathy_{y_raw:04d}_{m_raw:02d}_py.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def run_bathy_formatter(
    nc_dir: str | Path,
    out_dir: str | Path,
    plot_dir: str | Path,
    dx: float = 20.0,
    drop_year: int | None = 2008,
    interp_method: str = "linear",
    boundary_tightness: float = 3.0,
    overlap_erosion_cells: int = 0,
    output_format: str = "netcdf",
    bathy_cmap_name: str = "kg2",
    bathy_cmap_file: str | Path | None = None,
    extents_cmap_name: str = "viridis",
    extents_cmap_file: str | Path | None = None,
    extent_boundary_method: str = "mask",
    mask_nan_plots: bool = True,
) -> BathyGrid:
    """Run the full bathy formatter workflow (processing + plotting).
    
    Parameters
    ----------
    output_format : str
        Format for regridded output: "netcdf" (CF-compliant, default),
        "mat" (MATLAB struct), or "both".
    """
    process_result = process_bathy_formatter(
        nc_dir=nc_dir,
        out_dir=out_dir,
        dx=dx,
        drop_year=drop_year,
        interp_method=interp_method,
        boundary_tightness=boundary_tightness,
        overlap_erosion_cells=overlap_erosion_cells,
        output_format=output_format,
    )

    plot_bathy_formatter_outputs(
        process_result=process_result,
        plot_dir=plot_dir,
        bathy_cmap_name=bathy_cmap_name,
        bathy_cmap_file=bathy_cmap_file,
        extents_cmap_name=extents_cmap_name,
        extents_cmap_file=extents_cmap_file,
        extent_boundary_method=extent_boundary_method,
        boundary_tightness=boundary_tightness,
        overlap_erosion_cells=overlap_erosion_cells,
        mask_nan_plots=mask_nan_plots,
    )
    return process_result.bathy


def process_bathy_formatter(
    nc_dir: str | Path,
    out_dir: str | Path,
    dx: float = 20.0,
    drop_year: int | None = 2008,
    interp_method: str = "linear",
    boundary_tightness: float = 3.0,
    overlap_erosion_cells: int = 0,
    output_format: str = "netcdf",
) -> BathyProcessResult:
    """Run non-plotting bathy formatter steps and save outputs.
    
    Parameters
    ----------
    output_format : str
        Format for regridded output: "netcdf" (CF-compliant, default),
        "mat" (MATLAB struct), or "both".
    """
    out_path = Path(out_dir)
    _ensure_dir(out_path)

    surveys = load_raw_surveys(nc_dir)

    # Save full raw set (with all years).
    min_year, _ = _datenum_to_year_month(min(s.datenum for s in surveys))
    max_year, _ = _datenum_to_year_month(max(s.datenum for s in surveys))
    base = surveys[0].location
    save_raw_surveys_mat(surveys, out_path / f"{base}_{min_year}-{max_year}_bathyRaw_with2008_2023s_py.mat")

    # Optionally drop a year to improve overlap, mirroring MATLAB behavior.
    surveys_proc = surveys
    if drop_year is not None:
        surveys_proc = remove_year(surveys_proc, drop_year)

    min_year2, _ = _datenum_to_year_month(min(s.datenum for s in surveys_proc))
    max_year2, _ = _datenum_to_year_month(max(s.datenum for s in surveys_proc))
    save_raw_surveys_mat(surveys_proc, out_path / f"{base}_{min_year2}-{max_year2}_bathyRaw_2023s_py.mat")

    x_lims, y_lims, x_min, y_min = compute_domain_extents(surveys_proc)

    bathy = regrid_surveys(
        surveys_proc,
        dx=dx,
        x_min=x_min,
        y_min=y_min,
        interp_method=interp_method,
        boundary_tightness=boundary_tightness,
        overlap_erosion_cells=overlap_erosion_cells,
    )
    
    # Save regridded bathymetry in requested format(s)
    if output_format in ("mat", "both"):
        save_bathy_grid_mat(
            bathy,
            out_path / f"{base}_{min_year2}-{max_year2}_bathyReGrid_2023s_dx{int(dx)}m_py.mat",
        )
    if output_format in ("netcdf", "both"):
        save_bathy_grid_netcdf(
            bathy,
            out_path / f"{base}_{min_year2}-{max_year2}_bathyReGrid_2023s_dx{int(dx)}m_py.nc",
        )

    return BathyProcessResult(
        surveys_all=surveys,
        surveys_processed=surveys_proc,
        bathy=bathy,
        x_lims=x_lims,
        y_lims=y_lims,
        x_min=x_min,
        y_min=y_min,
        output_base=base,
        min_year_processed=min_year2,
        max_year_processed=max_year2,
    )

def plot_bathy_formatter_outputs(
    process_result: BathyProcessResult | str | Path,
    plot_dir: str | Path,
    bathy_cmap_name: str = "kg2",
    bathy_cmap_file: str | Path | None = None,
    extents_cmap_name: str = "viridis",
    extents_cmap_file: str | Path | None = None,
    extent_boundary_method: str = "mask",
    boundary_tightness: float = 3.0,
    overlap_erosion_cells: int = 0,
    mask_nan_plots: bool = True,
) -> None:
    """Generate bathy formatter QC plots from a process result or saved .mat/.nc file path.

    If a file path is provided, only regridded plots are generated because raw-survey
    extents and raw-survey plots require per-survey source data.
    """
    process_result, loaded_from_file = _as_bathy_process_result(process_result)

    plot_path = Path(plot_dir)
    _ensure_dir(plot_path)
    _ensure_dir(plot_path / "raw")
    _ensure_dir(plot_path / "regridded")

    extents_cmap = get_named_colormap(extents_cmap_name, extents_cmap_file)
    bathy_cmap = get_named_colormap(bathy_cmap_name, bathy_cmap_file)

    overlap_mask = None
    if extent_boundary_method.lower() == "mask" and not loaded_from_file:
        overlap_mask = compute_minimum_overlap_mask(
            process_result.surveys_processed,
            process_result.bathy.x,
            process_result.bathy.y,
            boundary_tightness=boundary_tightness,
            overlap_erosion_cells=overlap_erosion_cells,
        )

    if not loaded_from_file:
        plot_extents(
            process_result.surveys_processed,
            plot_path / f"{process_result.output_base}_bathyExtents_{process_result.min_year_processed}-{process_result.max_year_processed}_py.png",
            extents_cmap,
            extent_boundary_method=extent_boundary_method,
        )

        plot_raw_surveys(
            process_result.surveys_processed,
            process_result.x_lims,
            process_result.y_lims,
            plot_path / "raw",
            bathy_cmap,
            mask_nan=mask_nan_plots,
        )
    else:
        print("Loaded from file path: skipping raw extents/raw survey plots (source surveys unavailable).")

    plot_regridded_surveys(
        process_result.bathy,
        process_result.surveys_processed,
        plot_path / "regridded",
        bathy_cmap,
        extent_boundary_method=extent_boundary_method,
        overlap_mask=overlap_mask,
        mask_nan=mask_nan_plots,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert and regrid Bogue Inlet bathymetry .nc files")
    parser.add_argument("--nc-dir", required=True, help="Directory containing input .nc files")
    parser.add_argument("--out-dir", required=True, help="Output directory for .mat files")
    parser.add_argument("--plot-dir", required=True, help="Output directory for QC plots")
    parser.add_argument("--dx", type=float, default=20.0, help="Grid spacing in meters (default: 20)")
    parser.add_argument(
        "--drop-year",
        type=int,
        default=2008,
        help="Survey year to remove before regridding (default: 2008). Use -1 to keep all years.",
    )
    parser.add_argument(
        "--interp-method",
        choices=["linear", "nearest", "cubic"],
        default="linear",
        help="Interpolation method for regridding (default: linear, MATLAB-like).",
    )
    parser.add_argument(
        "--boundary-tightness",
        type=float,
        default=3.0,
        help="Boundary tightness for overlap domain; lower is tighter (default: 3.0).",
    )
    parser.add_argument(
        "--overlap-erosion-cells",
        type=int,
        default=0,
        help="Additional binary-erosion steps on overlap mask (default: 0).",
    )
    parser.add_argument(
        "--bathy-cmap",
        choices=["kg2", "vintage", "turbo", "viridis"],
        default="kg2",
        help="Colormap for raw/regridded bathy plots (default: kg2).",
    )
    parser.add_argument(
        "--bathy-cmap-file",
        default=None,
        help="Optional path to custom .clrmap file for bathy plots.",
    )
    parser.add_argument(
        "--extents-cmap",
        choices=["kg2", "vintage", "turbo", "viridis"],
        default="viridis",
        help="Colormap for extents plot and year colorbar (default: viridis).",
    )
    parser.add_argument(
        "--extents-cmap-file",
        default=None,
        help="Optional path to custom .clrmap file for extents plot.",
    )
    parser.add_argument(
        "--extent-boundary-method",
        choices=["mask", "convex"],
        default="mask",
        help="Extents-outline method: mask (tighter, default) or convex.",
    )
    parser.add_argument(
        "--no-mask-nan",
        action="store_true",
        help="If set, plot NaN areas instead of masking them (showing filled contours everywhere).",
    )
    parser.add_argument(
        "--output-format",
        choices=["netcdf", "mat", "both"],
        default="netcdf",
        help="Output format for regridded bathymetry: netcdf (CF-compliant, default), mat (MATLAB), or both.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    drop_year = None if args.drop_year == -1 else args.drop_year
    mask_nan = not args.no_mask_nan

    run_bathy_formatter(
        nc_dir=args.nc_dir,
        out_dir=args.out_dir,
        plot_dir=args.plot_dir,
        dx=args.dx,
        drop_year=drop_year,
        interp_method=args.interp_method,
        boundary_tightness=args.boundary_tightness,
        overlap_erosion_cells=args.overlap_erosion_cells,
        output_format=args.output_format,
        bathy_cmap_name=args.bathy_cmap,
        bathy_cmap_file=args.bathy_cmap_file,
        extents_cmap_name=args.extents_cmap,
        extents_cmap_file=args.extents_cmap_file,
        extent_boundary_method=args.extent_boundary_method,
        mask_nan_plots=mask_nan,
    )
    print("bathy_formatter workflow complete")


if __name__ == "__main__":
    main()
