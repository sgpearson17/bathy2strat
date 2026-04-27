"""Carolina Inlets analysis package."""

from .bathy import (
    BathyProcessResult,
    BathyGrid,
    RawSurvey,
    load_clrmap_file,
    load_bathy_grid,
    load_for_analysis,
    load_raw_surveys,
    plot_bathy_formatter_outputs,
    process_bathy_formatter,
    run_bathy_formatter,
    save_bathy_grid_netcdf,
)
from .stratigraphy import (
    BathyCube,
    StratigraphyConfig,
    StratigraphyResult,
    Transect,
    compute_stratigraphy,
    load_bathy_mat_known_structure,
    load_transects_from_shapefiles,
    run_stratigraphy_workflow,
)

__all__ = [
    "BathyProcessResult",
    "BathyGrid",
    "RawSurvey",
    "load_clrmap_file",
    "load_bathy_grid",
    "load_for_analysis",
    "load_raw_surveys",
    "plot_bathy_formatter_outputs",
    "process_bathy_formatter",
    "run_bathy_formatter",
    "save_bathy_grid_netcdf",
    "BathyCube",
    "StratigraphyConfig",
    "StratigraphyResult",
    "Transect",
    "compute_stratigraphy",
    "load_bathy_mat_known_structure",
    "load_transects_from_shapefiles",
    "run_stratigraphy_workflow",
]
