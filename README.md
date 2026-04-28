# bathy2strat
Code for computing stratigraphy from bathymetric data and other morphodynamic analyses, as per Pearson et al (2022).

## Python environment setup

Use a local virtual environment and install pinned dependencies so everyone runs the same stack.

### 1) Create and activate a virtual environment

Windows (PowerShell):

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

macOS/Linux (bash/zsh):

```bash
python -m venv .venv
source .venv/bin/activate
```

### 2) Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For strict reproducibility (CI or exact team sync), install from the lock file instead:

```bash
python -m pip install -r requirements-lock.txt
```

### One-command bootstrap scripts

From the repo root:

Windows (PowerShell):

```powershell
.\setup.ps1
```

Use exact lock versions:

```powershell
.\setup.ps1 -UseLockFile
```

macOS/Linux (bash/zsh):

```bash
bash setup.sh
```

Use exact lock versions:

```bash
bash setup.sh --lock
```

### 3) Register and select the notebook kernel

```bash
python -m ipykernel install --user --name bathy2strat --display-name "Python (.venv bathy2strat)"
```

In VS Code, open the notebook and select the kernel named `Python (.venv bathy2strat)`.

### 4) Verify core imports

```bash
python -c "import numpy, pandas, scipy, matplotlib, netCDF4, pyproj, geopandas; print('Environment OK')"
```

## Code

### Python stratigraphy workflow
Initial Python conversion of the stratigraphy workflow is available in `code/bathymetry_analysis/stratigraphy.py`.

Supported first-pass features:
- Known MATLAB bathy `.mat` input structure
- Transect sources: shapefile, manual points, and interactive GUI picking
- Configurable CRS transformation for transects
- Plot outputs to `plots/python/stratigraphy`
- Separate metric outputs under `metrics/stratigraphy`

Example run (from `code/`):

```bash
python -m bathymetry_analysis.stratigraphy \
	--bathy-mat ..\\data\\BogueInlet_2005-2023_bathyReGrid_2023s_dx20m.mat \
	--output-root .. \
	--transect-mode shapefile \
	--shp-dir C:\\surfdrive\\500_Analysis\\522_CarolinaInlets\\data\\SHP \
	--target-crs EPSG:32618
```

### depositAgeCalculator_v016_2022_Natascia.m
Code for computing stratigraphy from bathymetry as per Figure 5 of Pearson et al (2022)

### depositAgeCalculator_v019_BogueInlet_GUI_Ameland.m
Code for computing stratigraphy from bathymetry with GUI as per Pearson et al (2023)

### polarETDanalysis_v014_resubmission_final_newSensitivity_plot_Ale.m
Code to compute polar analysis of ebb-tidal delta changes as per Figure 4 of Pearson et al (2022)

### bathyChangeEnvelope_Weighted_v008.m
Code to compute weighted bathymetric change envelope as per Figure 3 of Pearson et al (2022)

### kg2Colormap.m
Custom colourmap developed by Stuart during the Kustgenese2.0 project. See here for more details: https://coastallycurious.com/2022/02/20/custom-colourmapping/

### vintageColormap.m
Custom colourmap developed by Stuart based on old Rijkswaterstaat maps. See here for more details: https://coastallycurious.com/2022/02/20/custom-colourmapping/

## References

Pearson, S.G., Elias, E.P.L, van Prooijen, B.C., van der Vegt, H., van der Spek, A. & Wang, Z.B. (2022). _A Novel Approach to Mapping Ebb-Tidal Delta Morphodynamics and Stratigraphy_. Geomorphology. https://doi.org/10.1016/j.geomorph.2022.108185.

Pearson, S.G., Mallinson, D., Brown, C., Mulligan, R. (2023). _Bathymetry-derived Stratigraphic Mapping of Bogue Inlet, NC_. Poster. Presented by Dave Mallinson. Geological Society of America Joint 72nd Annual Southeastern/58th Annual Northeastern Section Meeting - 2023.
8/25. https://doi.org/10.1130/abs/2023SE-385906.
