"""Synthetic stratigraphy utilities for scenario testing and plotting."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors

from .stratigraphy import (
    BathyCube,
    StratigraphyConfig,
    Transect,
    compute_stratigraphy,
    plot_cross_sections,
    plot_stratigraphy_stack,
)


def _matlab_datenum(dt: datetime) -> float:
    midnight = datetime(dt.year, dt.month, dt.day)
    return dt.toordinal() + 366 + (dt - midnight).total_seconds() / 86400.0


def _datenum_to_years(t_vals: np.ndarray) -> np.ndarray:
    years = []
    for value in np.asarray(t_vals, dtype=float).reshape(-1):
        ordinal = int(value)
        frac = value - ordinal
        py_dt = datetime.fromordinal(ordinal) + timedelta(days=frac) - timedelta(days=366)
        years.append(py_dt.year + (py_dt.timetuple().tm_yday - 1) / 365.0)
    return np.asarray(years, dtype=float)


def _make_grid(nx: int = 101, ny: int = 1, dx: float = 5.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # 1D transect grid in meters (ny=1 gives a single cross-shore line).
    x_axis = np.arange(nx, dtype=float) * dx
    y_axis = np.zeros(ny, dtype=float)
    x, y = np.meshgrid(x_axis, y_axis)
    return x, y, x_axis, y_axis


def _gaussian(x: np.ndarray, x0: float, sigma: float, amp: float) -> np.ndarray:
    # 1D Gaussian profile used for hump/channel scenarios.
    r2 = (x - x0) ** 2
    return amp * np.exp(-0.5 * r2 / (sigma**2))


@dataclass(frozen=True)
class SyntheticScenario:
    """Configuration for a synthetic stratigraphy scenario."""

    name: str
    nsteps: int
    translation_m_per_year: float


@dataclass(frozen=True)
class SyntheticScenarioSuite:
    """Collection of scenario configs with shared defaults."""

    scenarios: tuple[SyntheticScenario, ...] = field(default_factory=lambda: (
        SyntheticScenario("gaussian_hump_migration", nsteps=21, translation_m_per_year=10.0),
        SyntheticScenario("gaussian_channel_migration", nsteps=21, translation_m_per_year=10.0),
        SyntheticScenario("gaussian_channel_fill", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("gaussian_channel_growth", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("gaussian_hump_growth", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("gaussian_hump_erosion", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("gilbert_delta_ideal", nsteps=21, translation_m_per_year=10.0),
        SyntheticScenario("accreting_flat", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("eroding_flat", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("static_flat", nsteps=21, translation_m_per_year=0.0),
        SyntheticScenario("bruun_slr", nsteps=20, translation_m_per_year=0.0),
    ))


def build_synthetic_suite(
    include: Iterable[str] | None = None,
    overrides: dict[str, dict[str, float | int]] | None = None,
) -> SyntheticScenarioSuite:
    """Build a scenario suite with optional filtering and per-scenario overrides."""
    base = SyntheticScenarioSuite().scenarios
    if include is not None:
        include_set = {name for name in include}
        base = tuple(s for s in base if s.name in include_set)

    override_map: dict[str, dict[str, float | int]] = overrides or {}
    updated: list[SyntheticScenario] = []
    for scenario in base:
        override = override_map.get(scenario.name, {})
        updated.append(
            SyntheticScenario(
                scenario.name,
                nsteps=int(override.get("nsteps", scenario.nsteps)),
                translation_m_per_year=float(override.get("translation_m_per_year", scenario.translation_m_per_year)),
            )
        )
    return SyntheticScenarioSuite(tuple(updated))


def list_synthetic_scenarios() -> list[SyntheticScenario]:
    """Return the default synthetic scenario configurations."""
    return list(SyntheticScenarioSuite().scenarios)


def build_synthetic_bathy(
    scenario: str,
    nsteps: int = 2,
    dx: float = 5.0,
    domain_m: float = 500.0,
    width_m: float = 100.0,
    amp_m: float = 10.0,
    translation_m_per_year: float = 1.0,
    years_per_step: float = 1.0,
    bruun_years_per_step: float = 10.0,
    bruun_A: float = 0.05,
    bruun_slr_rate: float = 0.01,
    bruun_retreat_factor: float = 100.0,
    bruun_m: float = 2.0 / 3.0,
) -> BathyCube:
    scenario_aliases = {
        "gaussian_hump": "gaussian_hump_migration",
        "gaussian_channel": "gaussian_channel_migration",
    }
    scenario = scenario_aliases.get(scenario, scenario)

    if scenario == "bruun_slr":
        x_axis = np.arange(-300.0, 800.0 + dx, dx)
        y_axis = np.zeros(1, dtype=float)
        x, y = np.meshgrid(x_axis, y_axis)
        base = -5.0
        z = np.full((y.shape[0], x.shape[1], nsteps), base, dtype=float)
    else:
        nx = int(domain_m / dx) + 1
        x, y, x_axis, _y_axis = _make_grid(nx=nx, ny=1, dx=dx)
        base = -5.0
        z = np.full((y.shape[0], x.shape[1], nsteps), base, dtype=float)

    # Annual time stamps in MATLAB datenum format.
    t0 = datetime(2020, 1, 1)
    t_vals = np.array([_matlab_datenum(t0 + timedelta(days=365 * years_per_step * i)) for i in range(nsteps)], dtype=float)

    # Start near the left edge so the feature stays on-screen longer.
    x_center = max(0.5 * width_m, 0.3 * domain_m)
    sigma = max(width_m / 2.355, dx)
    shift_per_step = translation_m_per_year

    step_denom = max(nsteps - 1, 1)

    if scenario == "gaussian_hump_migration":
        for i in range(nsteps):
            x0 = x_center + i * shift_per_step
            z[0, :, i] = base + _gaussian(x_axis, x0, sigma, amp=amp_m)
    elif scenario == "gaussian_channel_migration":
        for i in range(nsteps):
            x0 = x_center + i * shift_per_step
            z[0, :, i] = base - _gaussian(x_axis, x0, sigma, amp=amp_m)
    elif scenario == "gaussian_channel_fill":
        for i in range(nsteps):
            scale = 1.0 - (i / step_denom)
            z[0, :, i] = base - _gaussian(x_axis, x_center, sigma, amp=amp_m * scale)
    elif scenario == "gaussian_channel_growth":
        for i in range(nsteps):
            scale = i / step_denom
            z[0, :, i] = base - _gaussian(x_axis, x_center, sigma, amp=amp_m * scale)
    elif scenario == "gaussian_hump_growth":
        for i in range(nsteps):
            scale = i / step_denom
            z[0, :, i] = base + _gaussian(x_axis, x_center, sigma, amp=amp_m * scale)
    elif scenario == "gaussian_hump_erosion":
        for i in range(nsteps):
            scale = 1.0 - (i / step_denom)
            z[0, :, i] = base + _gaussian(x_axis, x_center, sigma, amp=amp_m * scale)
    elif scenario == "gilbert_delta_ideal":
        top_elev = 0.0
        foreset_angle_deg = 25.0
        foreset_slope = -np.tan(np.deg2rad(foreset_angle_deg))
        foreset_len = 80.0
        toe_len = 30.0
        toe_slope = -0.02
        bottom_slope = -0.005
        x_start = 50.0
        prograde_per_step = max(translation_m_per_year, dx)
        for i in range(nsteps):
            x_shore = x_start + i * prograde_per_step
            x_foreset_end = x_shore + foreset_len
            x_toe_end = x_foreset_end + toe_len
            depth_foreset_end = top_elev + foreset_slope * foreset_len
            depth_toe_end = depth_foreset_end + toe_slope * toe_len

            profile = np.full_like(x_axis, top_elev)
            mask_fore = (x_axis >= x_shore) & (x_axis <= x_foreset_end)
            profile[mask_fore] = top_elev + foreset_slope * (x_axis[mask_fore] - x_shore)
            mask_toe = (x_axis > x_foreset_end) & (x_axis <= x_toe_end)
            profile[mask_toe] = depth_foreset_end + toe_slope * (x_axis[mask_toe] - x_foreset_end)
            mask_bottom = x_axis > x_toe_end
            profile[mask_bottom] = depth_toe_end + bottom_slope * (x_axis[mask_bottom] - x_toe_end)
            z[0, :, i] = profile
    elif scenario == "accreting_flat":
        rate = 0.2
        for i in range(nsteps):
            z[:, :, i] = base + rate * i
    elif scenario == "eroding_flat":
        rate = 0.2
        for i in range(nsteps):
            z[:, :, i] = base - rate * i
    elif scenario == "static_flat":
        pass
    elif scenario == "bruun_slr":
        # Dean profile with Bruun-rule retreat under SLR.
        A = bruun_A
        slr_rate = bruun_slr_rate
        retreat_factor = bruun_retreat_factor
        m_bruun = bruun_m
        years = np.arange(nsteps, dtype=float) * bruun_years_per_step
        t_vals = np.array([_matlab_datenum(t0 + timedelta(days=365 * y)) for y in years], dtype=float)
        for i, years_elapsed in enumerate(years):
            slr = slr_rate * years_elapsed
            retreat = slr * retreat_factor
            profile = -A * np.power(np.maximum(x_axis + retreat, 0.0), m_bruun)
            z[0, :, i] = profile + slr
    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    return BathyCube(location=f"Synthetic_{scenario}", t=t_vals, x=x, y=y, z=z)


def plot_time_colorbar(t_vals: np.ndarray, out_path: Path) -> None:
    # Standalone colorbar for time in years.
    years = _datenum_to_years(t_vals)
    vmin = np.nanmin(years)
    vmax = np.nanmax(years)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])

    fig, ax = plt.subplots(figsize=(6.5, 1.2))
    fig.subplots_adjust(bottom=0.35)
    cbar = fig.colorbar(sm, cax=ax, orientation="horizontal")
    cbar.set_label("Year")
    cbar.set_ticks(np.linspace(vmin, vmax, 5))
    cbar.ax.xaxis.set_major_formatter(lambda x, pos: f"{x:.1f}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def run_synthetic_scenario(
    scenario: str,
    cfg: StratigraphyConfig,
    output_root: Path,
    nsteps: int,
    translation_m_per_year: float,
    transects: Iterable[Transect] | None = None,
    start_year_at_zero: bool = True,
    years_per_step: float = 1.0,
    bruun_years_per_step: float = 20.0,
    bruun_A: float = 0.05,
    bruun_slr_rate: float = 0.01,
    bruun_retreat_factor: float = 100.0,
    bruun_m: float = 2.0 / 3.0,
    plot_xlim: tuple[float, float] | None = None,
    plot_ylim: tuple[float, float] | None = None,
) -> None:
    cube = build_synthetic_bathy(
        scenario,
        nsteps=nsteps,
        dx=cfg.dx,
        domain_m=500.0,
        width_m=100.0,
        amp_m=10.0,
        translation_m_per_year=translation_m_per_year,
        years_per_step=years_per_step,
        bruun_years_per_step=bruun_years_per_step,
        bruun_A=bruun_A,
        bruun_slr_rate=bruun_slr_rate,
        bruun_retreat_factor=bruun_retreat_factor,
        bruun_m=bruun_m,
    )
    result = compute_stratigraphy(cube, cfg)

    scenario_dir = output_root / scenario
    plot_stratigraphy_stack(
        cube,
        result,
        scenario_dir / "stratigraphy_stack.png",
        scenario,
        time_vals=cube.t,
        start_year_at_zero=start_year_at_zero,
        plot_xlim=plot_xlim,
        plot_ylim=plot_ylim,
    )
    plot_time_colorbar(cube.t, scenario_dir / "stratigraphy_time_colorbar.png")

    if transects:
        plot_cross_sections(cube, result, list(transects), scenario_dir)


def run_synthetic_scenario_by_name(
    name: str,
    cfg: StratigraphyConfig,
    output_root: Path,
    transects: Iterable[Transect] | None = None,
    nsteps: int | None = None,
    translation_m_per_year: float | None = None,
    start_year_at_zero: bool = True,
    years_total: float | None = None,
    years_step: float | None = None,
    plot_xlim: tuple[float, float] | None = None,
    plot_ylim: tuple[float, float] | None = None,
    bruun_A: float = 0.05,
    bruun_slr_rate: float = 0.01,
    bruun_retreat_factor: float = 100.0,
    bruun_m: float = 2.0 / 3.0,
) -> None:
    """Run a single synthetic scenario by name with optional overrides."""
    nsteps_local = nsteps
    years_step_local = years_step
    if name in {
        "gaussian_channel_fill",
        "gaussian_channel_growth",
        "gaussian_hump_growth",
        "gaussian_hump_erosion",
    }:
        if plot_xlim is None:
            plot_xlim = (0.0, 300.0)

    if name == "bruun_slr":
        if years_total is None:
            years_total = 200.0
        if years_step_local is None:
            years_step_local = 20.0
        if plot_xlim is None:
            plot_xlim = (-300.0, 800.0)
        if plot_ylim is None:
            plot_ylim = (-5.0, 3.0)
    if years_total is not None and years_step is not None:
        if years_step <= 0:
            raise ValueError("years_step must be positive")
        nsteps_local = int(round(years_total / years_step)) + 1

    overrides: dict[str, dict[str, float | int]] = {}
    if nsteps_local is not None or translation_m_per_year is not None:
        overrides[name] = {}
        if nsteps_local is not None:
            overrides[name]["nsteps"] = int(nsteps_local)
        if translation_m_per_year is not None:
            overrides[name]["translation_m_per_year"] = float(translation_m_per_year)

    suite = build_synthetic_suite(include=[name], overrides=overrides or None)
    if not suite.scenarios:
        raise ValueError(f"Unknown scenario: {name}")

    scenario = suite.scenarios[0]
    run_synthetic_scenario(
        scenario.name,
        cfg=cfg,
        output_root=output_root,
        nsteps=scenario.nsteps,
        translation_m_per_year=scenario.translation_m_per_year,
        transects=transects,
        start_year_at_zero=start_year_at_zero,
        years_per_step=years_step_local or 1.0,
        bruun_years_per_step=years_step_local or 10.0,
        plot_xlim=plot_xlim,
        plot_ylim=plot_ylim,
        bruun_A=bruun_A,
        bruun_slr_rate=bruun_slr_rate,
        bruun_retreat_factor=bruun_retreat_factor,
        bruun_m=bruun_m,
    )


def run_synthetic_scenarios(
    cfg: StratigraphyConfig,
    output_root: Path,
    transects: Iterable[Transect] | None = None,
    start_year_at_zero: bool = True,
    suite: SyntheticScenarioSuite | None = None,
) -> None:
    scenarios = suite.scenarios if suite is not None else SyntheticScenarioSuite().scenarios
    for scenario in scenarios:
        run_synthetic_scenario(
            scenario.name,
            cfg=cfg,
            output_root=output_root,
            nsteps=scenario.nsteps,
            translation_m_per_year=scenario.translation_m_per_year,
            transects=transects,
            start_year_at_zero=start_year_at_zero,
        )
