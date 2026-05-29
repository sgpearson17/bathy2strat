"""Tests for synthetic stratigraphy scenarios."""

from __future__ import annotations

import numpy as np

from bathymetry_analysis.synthetic_stratigraphy import (
    SyntheticScenarioSuite,
    build_synthetic_bathy,
    compute_deposit_elev_1d,
    compute_remaining_volumes,
)


def _scenario_metrics(scenario: str, nsteps: int = 6, translation: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    cube = build_synthetic_bathy(
        scenario,
        nsteps=nsteps,
        dx=5.0,
        domain_m=500.0,
        width_m=100.0,
        amp_m=10.0,
        translation_m_per_year=translation,
    )
    z_stack = cube.z[0, :, :].T
    deposit_elev = compute_deposit_elev_1d(z_stack)
    min_elev = np.nanmin(deposit_elev)
    base = min_elev - 0.1 * abs(min_elev)
    dx = float(cube.x[0, 1] - cube.x[0, 0]) if cube.x.shape[1] > 1 else 1.0
    remaining, initial_total = compute_remaining_volumes(z_stack, dx=dx, base=base)
    ratio = remaining[:, 0] / initial_total
    total_preserved = np.nansum(remaining, axis=1)
    return ratio, total_preserved


def test_static_flat_preservation_constant() -> None:
    ratio, total_preserved = _scenario_metrics("static_flat")
    assert np.allclose(ratio, 1.0, atol=1e-6)
    assert np.allclose(total_preserved, total_preserved[0], atol=1e-6)


def test_accreting_flat_increases_total_volume() -> None:
    ratio, total_preserved = _scenario_metrics("accreting_flat")
    assert total_preserved[-1] > total_preserved[0]
    assert ratio[-1] >= 1.0 - 1e-6


def test_eroding_flat_decreases_total_volume() -> None:
    ratio, total_preserved = _scenario_metrics("eroding_flat")
    assert total_preserved[-1] < total_preserved[0]
    assert ratio[-1] < 1.0


def test_gaussian_hump_preservation_nonincreasing() -> None:
    suite = SyntheticScenarioSuite()
    hump = next(s for s in suite.scenarios if s.name == "gaussian_hump")
    ratio, _total_preserved = _scenario_metrics(
        "gaussian_hump",
        nsteps=hump.nsteps,
        translation=hump.translation_m_per_year,
    )
    diffs = np.diff(ratio)
    assert np.all(diffs <= 1e-6)


def test_bruun_slr_profiles_deepen() -> None:
    suite = SyntheticScenarioSuite()
    bruun = next(s for s in suite.scenarios if s.name == "bruun_slr")
    cube = build_synthetic_bathy("bruun_slr", nsteps=bruun.nsteps, dx=5.0)
    z_stack = cube.z[0, :, :].T
    mean_depth = np.nanmean(z_stack, axis=1)
    diffs = np.diff(mean_depth)
    assert np.all(diffs <= 1e-6)
