import os
import numpy as np
import scipy.io as sio
import datetime
from scipy.interpolate import interp2d
from scipy.interpolate import splprep, splev
import matplotlib.pyplot as plt
import geopandas as gpd
from pyproj import Transformer

# Settings
data_dir = os.path.join('..', 'data')
plot_dir = os.path.join('..', 'plots', 'stratigraphy', 'shipTracks_v002')
ios.makedirs(plot_dir, exist_ok=True)

# Load bathymetric data (MAT-file)
mat = sio.loadmat(os.path.join(data_dir, 'BogueInlet_2005-2022_bathyReGrid_dx5m.mat'))
bathy = mat['bathy'].item()
X = bathy['x'] / 1000.0  # km
y = bathy['y'] / 1000.0
Z = bathy['z']           # m elevation
times = bathy['t'].flatten()  # MATLAB datenums
# Convert MATLAB datenum to Python datetime
T = [datetime.datetime.fromordinal(int(t)) + datetime.timedelta(days=t%1) - datetime.timedelta(days=366)
     for t in times]

# Tidal datums
MHW = 0.358
MSL = -0.112
MLW = -0.590

# Transect source: 2 = shapefile
transect_src = 2
ship_tracks = []
if transect_src == 2:
    shp_dir = r'C:\surfdrive\500_Analysis\522_CarolinaInlets\data\SHP'
    for shp in os.listdir(shp_dir):
        if shp.endswith('.shp'):
            gdf = gpd.read_file(os.path.join(shp_dir, shp))
            # Convert WGS84 to UTM zone 18N
            transformer = Transformer.from_crs(4326, 32618, always_xy=True)
            xs, ys = transformer.transform(gdf['geometry'].x.values, gdf['geometry'].y.values)
            ship_tracks.append({'name': shp[:-4], 'x': xs, 'y': ys})
    # filter short tracks
    ship_tracks = [st for st in ship_tracks if len(st['x']) > 1]

# Initialize parameters
dx = 20  # grid spacing (m?)
initial_idx = 1  # index start (Python 0-based adjust)

# Mask NaNs across all time steps
mask = np.any(np.isnan(Z), axis=2)
Z[mask[..., None]] = np.nan

# Compute min/max surfaces
min_surf = np.nanmin(Z[:, :, initial_idx:], axis=2)
max_surf = np.nanmax(Z[:, :, initial_idx:], axis=2)

# Compute deposit elevation time series
nz = Z.shape[2]
dep_elev = np.full_like(Z, np.nan)
dep_elev[:, :, initial_idx] = Z[:, :, initial_idx]

# March through time
for tt in range(initial_idx + 1, nz):
    dz = Z[:, :, tt] - Z[:, :, tt - 1]
    dep_elev[:, :, tt] = Z[:, :, tt]
    eros_idx = dz < 0
n    for (i, j), _ in np.ndenumerate(dz):
        if dz[i, j] < 0:
            for back in range(tt - initial_idx):
                prev = tt - back - 1
                if not np.isnan(dep_elev[i, j, prev]):
                    dep_elev[i, j, prev] = min(dep_elev[i, j, prev], dep_elev[i, j, tt])

# Convert elevations to thickness
dep_thk = np.diff(dep_elev, axis=2)
# prepend min surface
dep_thk = np.concatenate([min_surf[:, :, None], dep_thk[:, :, initial_idx:]], axis=2)

# Generate full monthly time series
start_year, end_year = T[0].year, T[-1].year
T_full = []
for yy in range(start_year, end_year + 1):
    for mm in range(1, 13):
        T_full.append(datetime.datetime(yy, mm, 1))
T_full = np.array(T_full)

# Map survey thickness onto full months
dep_thk_full = np.zeros((Z.shape[0], Z.shape[1], len(T_full))) * np.nan
survey_idx = 0
for idx, date in enumerate(T_full[:-1]):
    if date <= T[survey_idx] < T_full[idx + 1]:
        dep_thk_full[:, :, idx] = dep_thk[:, :, survey_idx]
        if survey_idx < len(T) - 1:
            survey_idx += 1

# Plotting example: final bathymetry and transects
plt.figure(figsize=(10, 8))
cs = plt.contourf(X, y, Z[:, :, -1], levels=np.arange(-30, 5.2, 0.2), cmap='viridis')
plt.colorbar(label='Elevation [m NAVD88]')
plt.contour(X, y, Z[:, :, -1], levels=[MLW], colors='gray', linewidths=1)
plt.contour(X, y, Z[:, :, -1], levels=[-6], linestyles=':', colors='k', linewidths=0.5)
plt.xlabel('Easting [km]')
plt.ylabel('Northing [km]')
plt.title(f"Transect Locations ({T[-1].year} Bathymetry)")
plt.gca().set_aspect('equal')
for st in ship_tracks:
    plt.plot(st['x'] / 1000, st['y'] / 1000, '-k')
plt.savefig(os.path.join(plot_dir, 'Stratigraphic_Transect_Map.png'), dpi=300)
plt.close()

# Note: cross-section plotting and interparc spline interpolation can be added similarly using
# shapely.geometry.LineString.interpolate and scipy.interpolate.splprep/splev.
