"""Benchmark script for bathy formatter regridding performance."""

from __future__ import annotations

from time import perf_counter
import argparse

from pathlib import Path
import sys

import numpy as np


def make_synthetic_surveys(n_surveys: int = 5):
    """Create a deterministic synthetic dataset for benchmarking.

    Args:
        n_surveys: Number of synthetic surveys to generate.

    Returns:
        List of :class:`bathy_formatter.RawSurvey` instances.
    """
    from bathy_formatter import RawSurvey, _matlab_datenum

    x_vals = np.arange(0.0, 200.0, 20.0)
    y_vals = np.arange(0.0, 200.0, 20.0)
    x_raw, y_raw = np.meshgrid(x_vals, y_vals)

    surveys: list[RawSurvey] = []
    for i in range(n_surveys):
        z_base = -1.0 + 0.05 * i
        z = z_base * np.ones_like(x_raw)
        z[0, 0] = np.nan  # consistent NaN to keep masks exercised
        surveys.append(
            RawSurvey(
                location="BenchUnit",
                datenum=_matlab_datenum(2000 + i, 1, 1),
                vertical_datum="NAVD88",
                horizontal_datum=None,
                ncfile=f"BenchUnit_{2000 + i:04d}_01_NAVD88.nc",
                x_raw=x_raw,
                y_raw=y_raw,
                z_raw=z,
            )
        )
    return surveys


def run_benchmark(run_real: bool) -> None:
    """Run synthetic and optional real-data benchmarks.

    Args:
        run_real: If True, include the notebook-sized dataset benchmark.
    """
    code_dir = Path(__file__).resolve().parents[1]
    if str(code_dir) not in sys.path:
        sys.path.insert(0, str(code_dir))

    from bathy_formatter import compute_domain_extents, load_raw_surveys, regrid_surveys

    def _bench_case(label: str, surveys: list, dx: float) -> None:
        """Time a single regridding scenario and print summary stats."""
        _, _, x_min, y_min = compute_domain_extents(surveys)

        start = perf_counter()
        bathy, overlap_mask = regrid_surveys(
            surveys,
            dx=dx,
            x_min=x_min,
            y_min=y_min,
            interp_method="linear",
            boundary_tightness=2.0,
            overlap_erosion_cells=1,
            parallel=False,
        )
        elapsed = perf_counter() - start

        print(f"{label} (serial)")
        print(f"surveys: {len(surveys)}")
        print(f"grid shape (y, x): {bathy.x.shape}")
        print(f"overlap mask true count: {int(np.sum(overlap_mask))}")
        print(f"elapsed seconds: {elapsed:.4f}")

        try:
            import joblib  # noqa: F401
        except ModuleNotFoundError:
            print(f"{label} (parallel) skipped: joblib not installed")
            return

        start = perf_counter()
        bathy, overlap_mask = regrid_surveys(
            surveys,
            dx=dx,
            x_min=x_min,
            y_min=y_min,
            interp_method="linear",
            boundary_tightness=2.0,
            overlap_erosion_cells=1,
            parallel=True,
        )
        elapsed = perf_counter() - start

        print(f"{label} (parallel)")
        print(f"surveys: {len(surveys)}")
        print(f"grid shape (y, x): {bathy.x.shape}")
        print(f"overlap mask true count: {int(np.sum(overlap_mask))}")
        print(f"elapsed seconds: {elapsed:.4f}")

    _bench_case("Synthetic benchmark", make_synthetic_surveys(), dx=20.0)

    if run_real:
        nc_dir = Path(r"C:\surf\300_Data\320_BogueInlet\netcdf\UTM_renamed")
        if nc_dir.exists():
            real_surveys = load_raw_surveys(nc_dir)
            _bench_case("Notebook benchmark", real_surveys, dx=20.0)
        else:
            print("Notebook benchmark skipped: nc_dir not found")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark bathy_formatter regridding")
    parser.add_argument(
        "--real",
        action="store_true",
        help="Include the notebook (real-data) benchmark.",
    )
    args = parser.parse_args()
    run_benchmark(run_real=args.real)
