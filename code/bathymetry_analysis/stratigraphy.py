from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Colormap
from pyproj import CRS, Transformer
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat


MHW = 0.358
MSL = -0.112
MLW = -0.590


@dataclass
class BathyCube:
    location: str
    t: np.ndarray  # MATLAB datenum, shape (nt,)
    x: np.ndarray  # 2D grid, meters
    y: np.ndarray  # 2D grid, meters
    z: np.ndarray  # 3D cube, shape (ny, nx, nt) or (nx, ny, nt)


@dataclass
class Transect:
    name: str
    x: np.ndarray  # meters
    y: np.ndarray  # meters


@dataclass
class StratigraphyConfig:
    initial_index: int = 0
    dx: float = 20.0
    target_crs: str = "EPSG:32618"
    dot_spacing_km: float = 0.1
    end_dot_size: float = 40.0
    mid_dot_size: float = 15.0


@dataclass
class StratigraphyResult:
    t: np.ndarray
    t_full: np.ndarray
    min_surf: np.ndarray
    max_surf: np.ndarray
    deposit_elev: np.ndarray
    deposit_thk_full: np.ndarray
    deposit_per_year: np.ndarray
    total_sed_vol_per_year: np.ndarray
    theseus_ratio: np.ndarray


def datenum_to_datetime64(datenum: np.ndarray) -> np.ndarray:
    import datetime as dt

    out = []
    for value in np.asarray(datenum, dtype=float).reshape(-1):
        ordinal = int(value)
        frac = value - ordinal
        py_dt = dt.datetime.fromordinal(ordinal) + dt.timedelta(days=frac) - dt.timedelta(days=366)
        out.append(np.datetime64(py_dt))
    return np.asarray(out)


def datetime64_to_datenum(values: Iterable[np.datetime64]) -> np.ndarray:
    import datetime as dt

    out = []
    for value in values:
        py_dt = pd.Timestamp(value).to_pydatetime()
        out.append(py_dt.toordinal() + 366 + (py_dt.hour / 24.0) + (py_dt.minute / 1440.0) + (py_dt.second / 86400.0))
    return np.asarray(out, dtype=float)


def _coerce_mat_struct(value):
    if isinstance(value, np.ndarray) and value.dtype == object:
        return value.flat[0]
    return value


def load_bathy_mat_known_structure(path: str | Path) -> BathyCube:
    mat = loadmat(str(path), squeeze_me=True, struct_as_record=False)
    if "bathy" not in mat:
        raise ValueError(f"Expected a 'bathy' struct in {path}")

    b = _coerce_mat_struct(mat["bathy"])
    location = str(getattr(b, "location", Path(path).stem))
    t = np.asarray(getattr(b, "t"), dtype=float).reshape(-1)
    x = np.asarray(getattr(b, "x"), dtype=float)
    y = np.asarray(getattr(b, "y"), dtype=float)
    z = np.asarray(getattr(b, "z"), dtype=float)

    if x.shape != y.shape:
        raise ValueError(f"x and y shapes differ: {x.shape} vs {y.shape}")
    if z.ndim != 3:
        raise ValueError(f"Expected 3D z array. Got shape: {z.shape}")

    # Align z first two dimensions to x/y if needed.
    if z.shape[:2] != x.shape and z.shape[:2] == x.T.shape:
        z = np.transpose(z, (1, 0, 2))
    if z.shape[:2] != x.shape:
        raise ValueError(f"z spatial shape {z.shape[:2]} does not match x/y shape {x.shape}")

    return BathyCube(location=location, t=t, x=x, y=y, z=z)


def _derive_xy_axes(x2d: np.ndarray, y2d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_axis = np.asarray(x2d[0, :], dtype=float)
    y_axis = np.asarray(y2d[:, 0], dtype=float)

    if not (np.all(np.diff(x_axis) > 0) and np.all(np.diff(y_axis) > 0)):
        raise ValueError("Only monotonic regular grids are currently supported for transect interpolation")

    return x_axis, y_axis


def _nanminmax_surface(stack: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    all_nan = np.all(np.isnan(stack), axis=axis)
    stack_safe = np.where(all_nan[..., None], np.inf, stack)
    min_surf = np.nanmin(stack_safe, axis=axis)
    min_surf[all_nan] = np.nan

    stack_safe = np.where(all_nan[..., None], -np.inf, stack)
    max_surf = np.nanmax(stack_safe, axis=axis)
    max_surf[all_nan] = np.nan
    return min_surf, max_surf


def compute_stratigraphy(bathy: BathyCube, config: StratigraphyConfig) -> StratigraphyResult:
    initial_idx = config.initial_index
    x = bathy.x
    y = bathy.y
    z = np.asarray(bathy.z, dtype=float).copy()
    t = np.asarray(bathy.t, dtype=float).reshape(-1)

    if initial_idx < 0 or initial_idx >= z.shape[2]:
        raise ValueError(f"initial_index out of bounds: {initial_idx}")

    # If any survey at a cell is NaN, keep it NaN across all surveys.
    any_nan = np.any(np.isnan(z), axis=2)
    z[any_nan] = np.nan

    min_surf, max_surf = _nanminmax_surface(z[:, :, initial_idx:], axis=2)

    deposit_elev = np.full_like(z, np.nan)
    deposit_elev[:, :, initial_idx] = z[:, :, initial_idx]

    nt = z.shape[2]
    deposit_per_year = np.zeros((nt, nt), dtype=float)
    t_full_global = np.array([], dtype=float)
    deposit_thk_full_global = np.zeros((z.shape[0], z.shape[1], 0), dtype=float)

    for t_deposit_year in range(initial_idx + 1, nt):
        for tt in range(initial_idx + 1, t_deposit_year + 1):
            dz = z[:, :, tt] - z[:, :, tt - 1]
            deposit_elev[:, :, tt] = z[:, :, tt]

            erosion_mask = dz < 0
            if np.any(erosion_mask):
                cur = deposit_elev[:, :, tt]
                for k in range(initial_idx, tt):
                    prev = deposit_elev[:, :, k]
                    update = erosion_mask & ~np.isnan(prev) & (cur < prev)
                    prev[update] = cur[update]
                    deposit_elev[:, :, k] = prev

        deposit_thk = np.diff(deposit_elev, axis=2)
        if initial_idx > 0:
            deposit_thk = deposit_thk[:, :, initial_idx - 1 :]
        deposit_thk = np.concatenate([min_surf[:, :, None], deposit_thk], axis=2)

        survey_datetimes = datenum_to_datetime64(t)
        start_year = pd.Timestamp(survey_datetimes[0]).year
        end_year = pd.Timestamp(survey_datetimes[-1]).year

        t_full = []
        for year in range(start_year, end_year + 1):
            for month in range(1, 13):
                t_full.append(np.datetime64(f"{year:04d}-{month:02d}-01"))
        t_full = np.asarray(t_full)
        t_full_dn = datetime64_to_datenum(t_full)

        deposit_thk_full = np.zeros((x.shape[0], x.shape[1], len(t_full_dn)), dtype=float)

        survey_count = 0
        for tt in range(len(t_full_dn) - 1):
            if survey_count >= len(t):
                break
            if t_full_dn[tt] <= t[survey_count] <= t_full_dn[tt + 1]:
                if survey_count < deposit_thk.shape[2]:
                    deposit_thk_full[:, :, tt] = deposit_thk[:, :, survey_count]
                if survey_count >= len(t) - 1:
                    break
                survey_count += 1

        deposit_per_year[:, t_deposit_year] = np.nansum(deposit_thk * (config.dx**2), axis=(0, 1))

        t_full_global = t_full_dn
        deposit_thk_full_global = deposit_thk_full

    z_min = np.nanmin(z[:, :, initial_idx:])
    total_sed_vol_per_year = np.nansum((z - z_min) * (config.dx**2), axis=(0, 1))

    if nt > 0:
        deposit_per_year[0, 0] = total_sed_vol_per_year[0]

    theseus_ratio = np.full((nt, nt), np.nan, dtype=float)
    for tt in range(nt):
        deposit_per_year[0, tt] = total_sed_vol_per_year[tt] - np.nansum(deposit_per_year[1:, tt])
        denom = deposit_per_year[tt, tt]
        if np.isfinite(denom) and denom != 0:
            theseus_ratio[tt, :] = deposit_per_year[tt, :] / denom

    return StratigraphyResult(
        t=t,
        t_full=t_full_global,
        min_surf=min_surf,
        max_surf=max_surf,
        deposit_elev=deposit_elev,
        deposit_thk_full=deposit_thk_full_global,
        deposit_per_year=deposit_per_year,
        total_sed_vol_per_year=total_sed_vol_per_year,
        theseus_ratio=theseus_ratio,
    )


def _extract_line_xy(gdf) -> tuple[np.ndarray, np.ndarray] | None:
    for geom in gdf.geometry:
        if geom is None:
            continue
        if geom.geom_type == "LineString":
            coords = np.asarray(geom.coords)
            return coords[:, 0], coords[:, 1]
        if geom.geom_type == "MultiLineString":
            all_coords = []
            for line in geom.geoms:
                all_coords.extend(list(line.coords))
            coords = np.asarray(all_coords)
            return coords[:, 0], coords[:, 1]
    return None


def _extract_xy_columns(gdf) -> tuple[np.ndarray, np.ndarray] | None:
    lower_to_col = {c.lower(): c for c in gdf.columns}

    x_keys = ["lon1", "x", "easting", "lon", "longitude"]
    y_keys = ["lat1", "y", "northing", "lat", "latitude"]

    x_col = next((lower_to_col[k] for k in x_keys if k in lower_to_col), None)
    y_col = next((lower_to_col[k] for k in y_keys if k in lower_to_col), None)

    if x_col is None or y_col is None:
        return None

    x = pd.to_numeric(gdf[x_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(gdf[y_col], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    return x[valid], y[valid]


def _resolve_source_crs(gdf, source_crs: str | None, prompt_if_missing: bool, shp_name: str) -> CRS:
    if source_crs is not None:
        return CRS.from_user_input(source_crs)
    if gdf.crs is not None:
        return CRS.from_user_input(gdf.crs)
    if not prompt_if_missing:
        raise ValueError(f"Shapefile {shp_name} has no CRS and source_crs was not provided")

    response = input(
        f"Shapefile {shp_name} has no CRS metadata. Enter source CRS (e.g., EPSG:4326 or EPSG:32618): "
    ).strip()
    if not response:
        raise ValueError(f"No CRS provided for shapefile {shp_name}")
    return CRS.from_user_input(response)


def load_transects_from_shapefiles(
    shp_dir: str | Path,
    target_crs: str,
    source_crs: str | None = None,
    prompt_if_missing_crs: bool = True,
) -> list[Transect]:
    try:
        import geopandas as gpd
    except ImportError as exc:
        raise ImportError("geopandas is required for shapefile transect loading") from exc

    shp_paths = sorted(Path(shp_dir).glob("*.shp"))
    transects: list[Transect] = []

    if not shp_paths:
        raise FileNotFoundError(f"No shapefiles found in {shp_dir}")

    dst_crs = CRS.from_user_input(target_crs)

    for shp in shp_paths:
        gdf = gpd.read_file(shp)

        xy = _extract_line_xy(gdf)
        if xy is None:
            xy = _extract_xy_columns(gdf)
        if xy is None:
            raise ValueError(
                f"Could not derive transect coordinates from {shp.name}. "
                "Expected LineString geometry or columns like Lon1/Lat1 or x/y."
            )

        x_raw, y_raw = xy
        if len(x_raw) < 2:
            continue

        src_crs_obj = _resolve_source_crs(gdf, source_crs, prompt_if_missing_crs, shp.name)
        transformer = Transformer.from_crs(src_crs_obj, dst_crs, always_xy=True)
        x_m, y_m = transformer.transform(x_raw, y_raw)

        transects.append(Transect(name=shp.stem, x=np.asarray(x_m, dtype=float), y=np.asarray(y_m, dtype=float)))

    return transects


def build_manual_transects(transect_points: list[np.ndarray], names: list[str] | None = None) -> list[Transect]:
    out: list[Transect] = []
    for i, points in enumerate(transect_points):
        arr = np.asarray(points, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("Each transect in transect_points must be an Nx2 array")
        name = names[i] if names is not None and i < len(names) else f"T{i + 1}"
        out.append(Transect(name=name, x=arr[:, 0], y=arr[:, 1]))
    return out


def pick_interactive_transect(ax, name: str = "GUI") -> Transect:
    print("Left-click to add points, press Enter to finish transect")
    pts = np.asarray(plt.ginput(n=-1, timeout=0), dtype=float)
    if pts.shape[0] < 2:
        raise ValueError("At least two points are required for a transect")
    ax.plot(pts[:, 0], pts[:, 1], "xk")
    return Transect(name=name, x=pts[:, 0], y=pts[:, 1])


def _resample_polyline(x: np.ndarray, y: np.ndarray, n_points: int = 1000) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    keep = np.ones_like(x, dtype=bool)
    keep[1:] = (np.diff(x) != 0) | (np.diff(y) != 0)
    x = x[keep]
    y = y[keep]

    dx = np.diff(x)
    dy = np.diff(y)
    seg_len = np.hypot(dx, dy)
    s = np.concatenate([[0.0], np.cumsum(seg_len)])

    if s[-1] <= 0:
        raise ValueError("Transect has zero length")

    s_q = np.linspace(0.0, s[-1], n_points)
    x_q = np.interp(s_q, s, x)
    y_q = np.interp(s_q, s, y)
    return x_q, y_q, s_q


def _interp_stack_along_transect(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    stack: np.ndarray,
    x_q: np.ndarray,
    y_q: np.ndarray,
) -> np.ndarray:
    out = np.full((stack.shape[2], len(x_q)), np.nan, dtype=float)
    pts = np.column_stack([y_q, x_q])

    for k in range(stack.shape[2]):
        interp = RegularGridInterpolator(
            (y_axis, x_axis),
            stack[:, :, k],
            bounds_error=False,
            fill_value=np.nan,
        )
        out[k, :] = interp(pts)

    return out


def _interp_surface_along_transect(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    surface: np.ndarray,
    x_q: np.ndarray,
    y_q: np.ndarray,
) -> np.ndarray:
    interp = RegularGridInterpolator((y_axis, x_axis), surface, bounds_error=False, fill_value=np.nan)
    pts = np.column_stack([y_q, x_q])
    return interp(pts)


def _default_cmap(n: int) -> Colormap:
    return plt.get_cmap("viridis", max(n, 2))


def plot_transect_location_map(
    bathy: BathyCube,
    transects: list[Transect],
    out_path: str | Path,
    config: StratigraphyConfig,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 8), dpi=150)
    x_km = bathy.x / 1000.0
    y_km = bathy.y / 1000.0

    cont = ax.contourf(x_km, y_km, bathy.z[:, :, -1], levels=np.arange(-30, 5.2, 0.2), cmap="viridis")
    ax.contour(x_km, y_km, bathy.z[:, :, -1], levels=[MLW], colors=[(0.5, 0.5, 0.5)], linewidths=1.0)
    ax.contour(x_km, y_km, bathy.z[:, :, -1], levels=[-6], colors="k", linestyles=":", linewidths=0.5)

    for tr in transects:
        xq, yq, dq = _resample_polyline(tr.x, tr.y, n_points=1000)
        xq_km = xq / 1000.0
        yq_km = yq / 1000.0
        dq_km = dq / 1000.0

        ax.plot(xq_km, yq_km, "-k", lw=1.0)
        ax.scatter([xq_km[0]], [yq_km[0]], s=config.end_dot_size, c="w", edgecolors="k", zorder=3)
        ax.scatter([xq_km[-1]], [yq_km[-1]], s=config.end_dot_size, c="w", edgecolors="k", zorder=3)

        prev_bin = -1
        for i in range(len(dq_km)):
            cur_bin = int(np.floor(dq_km[i] / config.dot_spacing_km))
            if cur_bin > prev_bin:
                ax.scatter([xq_km[i]], [yq_km[i]], s=config.mid_dot_size, c="w", edgecolors="k", zorder=3)
                prev_bin = cur_bin

        ax.text(xq_km[0] + 0.03, yq_km[0] - 0.03, tr.name, fontsize=10, fontweight="bold")
        ax.text(xq_km[-1] - 0.08, yq_km[-1] + 0.03, f"{tr.name}'", fontsize=10, fontweight="bold")

    cb = fig.colorbar(cont, ax=ax)
    cb.set_label("Elevation [m NAVD88]")

    survey_year = pd.Timestamp(datenum_to_datetime64(np.asarray([bathy.t[-1]]))[0]).year
    ax.set_title(f"Transect Locations ({survey_year} Bathymetry)")
    ax.set_xlabel("Easting [km]")
    ax.set_ylabel("Northing [km]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_cross_sections(
    bathy: BathyCube,
    result: StratigraphyResult,
    transects: list[Transect],
    out_dir: str | Path,
) -> None:
    x_axis, y_axis = _derive_xy_axes(bathy.x, bathy.y)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmap = _default_cmap(result.deposit_thk_full.shape[2])
    colors = cmap(np.linspace(0, 1, max(result.deposit_thk_full.shape[2], 2)))

    for tr in transects:
        xq, yq, dq_m = _resample_polyline(tr.x, tr.y, n_points=1000)
        dq_km = dq_m / 1000.0

        z_min = _interp_surface_along_transect(x_axis, y_axis, result.min_surf, xq, yq)
        z_max = _interp_surface_along_transect(x_axis, y_axis, result.max_surf, xq, yq)
        z_now = _interp_surface_along_transect(x_axis, y_axis, bathy.z[:, :, -1], xq, yq)
        dep_xs = _interp_stack_along_transect(x_axis, y_axis, result.deposit_thk_full, xq, yq)

        fig, ax = plt.subplots(figsize=(14, 4), dpi=150)

        baseline = -35.0 * np.ones_like(dq_km)
        running = baseline.copy()
        for k in range(dep_xs.shape[0]):
            layer = np.nan_to_num(dep_xs[k, :], nan=0.0)
            next_running = running + layer
            ax.fill_between(dq_km, running, next_running, color=colors[k], linewidth=0.0)
            running = next_running

        ax.plot(dq_km, z_min, "-k", lw=1.0)
        ax.plot(dq_km, z_max, ":k", lw=1.0)
        ax.plot(dq_km, z_now, "-k", lw=1.5)

        ax.axhline(MLW, linestyle="--", color=(0.1, 0.2, 0.5), linewidth=0.6)
        ax.axhline(MHW, linestyle="--", color=(0.1, 0.2, 0.5), linewidth=0.6)
        ax.text(dq_km[0] + 0.02, MHW + 0.3, "MHW", color=(0.1, 0.2, 0.5), fontsize=9)
        ax.text(dq_km[0] + 0.02, MLW + 0.3, "MLW", color=(0.1, 0.2, 0.5), fontsize=9)

        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("Elevation [m NAVD88]")
        ax.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)
        ax.set_xlim(0.0, np.nanmax(dq_km))

        y_min = min(np.nanmin(z_min) - 1.0, np.nanmin(baseline) - 1.0)
        y_max = np.nanmax([np.nanmax(z_now) + 1.0, 3.0])
        if np.isfinite(y_min) and np.isfinite(y_max) and y_max > y_min:
            ax.set_ylim(y_min, y_max)

        fig.tight_layout()
        fig.savefig(out_dir / f"Cross-section {tr.name}-{tr.name}'.png")
        plt.close(fig)


def plot_theseus_ratio(result: StratigraphyResult, out_path: str | Path) -> None:
    t_dt = datenum_to_datetime64(result.t)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=150)

    for tt in range(result.theseus_ratio.shape[0]):
        y = result.theseus_ratio[tt, tt:]
        x = t_dt[tt:]
        valid = np.isfinite(y)
        if np.any(valid):
            ax1.plot(x[valid], y[valid], linewidth=1.0)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Theta (Fraction Preserved) [-]")
    ax1.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)

    for tt in range(result.theseus_ratio.shape[0]):
        y = result.theseus_ratio[tt, tt:]
        x = (result.t[tt:] - result.t[tt]) / 10.0
        valid = np.isfinite(y) & np.isfinite(x) & (x > 0) & (y > 0)
        if np.any(valid):
            ax2.plot(x[valid], y[valid], linewidth=1.0)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Theta (Fraction Preserved) [-]")
    ax2.grid(True, color=(0.5, 0.5, 0.5), alpha=0.4)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def save_metrics(result: StratigraphyResult, out_dir: str | Path) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(result.deposit_per_year).to_csv(out_dir / "deposit_per_year.csv", index=False)
    pd.DataFrame(result.theseus_ratio).to_csv(out_dir / "theseus_ratio.csv", index=False)
    pd.DataFrame({"total_sed_vol_per_year": result.total_sed_vol_per_year}).to_csv(
        out_dir / "total_sed_vol_per_year.csv", index=False
    )
    pd.DataFrame({"survey_datenum": result.t}).to_csv(out_dir / "survey_datenum.csv", index=False)
    pd.DataFrame({"monthly_datenum": result.t_full}).to_csv(out_dir / "monthly_datenum.csv", index=False)

    np.savez_compressed(out_dir / "deposit_thk_full.npz", deposit_thk_full=result.deposit_thk_full)
    np.savez_compressed(out_dir / "min_surf.npz", min_surf=result.min_surf)
    np.savez_compressed(out_dir / "max_surf.npz", max_surf=result.max_surf)


def run_stratigraphy_workflow(
    bathy_mat_path: str | Path,
    output_root: str | Path,
    transect_mode: str,
    shp_dir: str | Path | None = None,
    manual_transects: list[np.ndarray] | None = None,
    source_crs: str | None = None,
    config: StratigraphyConfig | None = None,
) -> StratigraphyResult:
    config = config or StratigraphyConfig()
    bathy = load_bathy_mat_known_structure(bathy_mat_path)
    result = compute_stratigraphy(bathy, config)

    output_root = Path(output_root)
    strat_dir = output_root / "plots" / "python" / "stratigraphy"
    metrics_dir = output_root / "metrics" / "stratigraphy"

    transects: list[Transect] = []
    if transect_mode == "shapefile":
        if shp_dir is None:
            raise ValueError("shp_dir is required when transect_mode='shapefile'")
        transects = load_transects_from_shapefiles(
            shp_dir=shp_dir,
            target_crs=config.target_crs,
            source_crs=source_crs,
            prompt_if_missing_crs=True,
        )
    elif transect_mode == "manual":
        if not manual_transects:
            raise ValueError("manual_transects is required when transect_mode='manual'")
        transects = build_manual_transects(manual_transects)
    elif transect_mode == "gui":
        fig, ax = plt.subplots(figsize=(12, 8), dpi=150)
        ax.contourf(bathy.x / 1000.0, bathy.y / 1000.0, bathy.z[:, :, -1], levels=np.arange(-30, 5.2, 0.2))
        ax.set_title("Click transect points, then press Enter")
        tr = pick_interactive_transect(ax=ax, name="GUI")
        plt.close(fig)
        transects = [tr]
    else:
        raise ValueError("transect_mode must be one of: shapefile, manual, gui")

    if transects:
        plot_transect_location_map(
            bathy,
            transects,
            strat_dir / "Stratigraphic Transect Map (Python).png",
            config,
        )
        plot_cross_sections(bathy, result, transects, strat_dir)

    plot_theseus_ratio(result, strat_dir / "TheseusRatio.png")
    save_metrics(result, metrics_dir)

    return result


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute stratigraphy from MATLAB bathy input")
    parser.add_argument("--bathy-mat", required=True, help="Path to bathy MAT file with known structure")
    parser.add_argument("--output-root", default=".", help="Root output folder")
    parser.add_argument(
        "--transect-mode",
        required=True,
        choices=["shapefile", "manual", "gui"],
        help="Transect source mode",
    )
    parser.add_argument("--shp-dir", default=None, help="Folder of shapefiles when --transect-mode shapefile")
    parser.add_argument("--source-crs", default=None, help="Optional source CRS override, e.g. EPSG:4326")
    parser.add_argument("--target-crs", default="EPSG:32618", help="Target projected CRS")
    parser.add_argument("--initial-index", type=int, default=0, help="Initial survey index used as baseline")
    parser.add_argument("--dx", type=float, default=20.0, help="Grid cell size [m]")
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    config = StratigraphyConfig(
        initial_index=args.initial_index,
        dx=args.dx,
        target_crs=args.target_crs,
    )

    run_stratigraphy_workflow(
        bathy_mat_path=args.bathy_mat,
        output_root=args.output_root,
        transect_mode=args.transect_mode,
        shp_dir=args.shp_dir,
        source_crs=args.source_crs,
        config=config,
    )


if __name__ == "__main__":
    main()
