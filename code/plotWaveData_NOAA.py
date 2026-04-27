"""
Plot wave data from NOAA buoys 41159 and 41110
Reads meteorological data files and plots significant wave height time series

v001 - SGP/2026-04-22 - Converted from MATLAB
"""

import os
import glob
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime

# Set up matplotlib with Arial Bold Italic font for consistent styling
available_fonts = set(matplotlib.font_manager.get_font_names())
if 'Arial' in available_fonts:
    font_name = 'Arial'
else:
    font_name = 'DejaVu Sans'
    print("Note: Arial font not found. Using DejaVu Sans as fallback.")

plt.rcParams['font.family'] = font_name
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.style'] = 'italic'

# Define data directories
data_dirs = [
    r'C:\surf\500_Analysis\522_CarolinaInlets\data\meteo\NOAA_41159',
    r'C:\surf\500_Analysis\522_CarolinaInlets\data\meteo\NOAA_41110'
]

# Output directory for plots
plot_dir = r'C:\surf\500_Analysis\522_CarolinaInlets\plots'

# Create column names for the data
column_names = ['YY', 'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST', 'WVHT', 
                'DPD', 'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']

# Dictionary to store data for each buoy
buoy_data = {}

# Process each directory
for data_dir in data_dirs:
    # Extract buoy ID from directory name
    buoy_id = os.path.basename(data_dir)
    
    # Get all text files in the directory (skip README files)
    txt_files = sorted(glob.glob(os.path.join(data_dir, '*.txt')))
    txt_files = [f for f in txt_files if 'README' not in f.upper()]
    
    print(f"Processing {buoy_id}: Found {len(txt_files)} data files")
    
    # Initialize list to store all dataframes
    dfs = []
    
    # Read each file
    for file_path in txt_files:
        try:
            # Skip the header lines and read the data
            # NOAA files have comment lines starting with #
            df = pd.read_csv(file_path, sep=r'\s+', 
                           comment='#', names=column_names,
                           dtype={'WVHT': float}, na_values=99.0)
            
            # Create datetime from components
            df['datetime'] = pd.to_datetime(
                df[['YY', 'MM', 'DD', 'hh', 'mm']].rename(
                    columns={'YY': 'year', 'MM': 'month', 'DD': 'day', 
                            'hh': 'hour', 'mm': 'minute'}
                )
            )
            
            dfs.append(df)
            print(f"  Loaded: {os.path.basename(file_path)}")
            
        except Exception as e:
            print(f"  Error reading {file_path}: {e}")
    
    # Concatenate all dataframes for this buoy
    if dfs:
        meteo_data = pd.concat(dfs, ignore_index=True)
        meteo_data = meteo_data.sort_values('datetime').reset_index(drop=True)
        buoy_data[buoy_id] = {
            'time': meteo_data['datetime'].values,
            'Hs': meteo_data['WVHT'].values,
            'buoy_id': buoy_id.replace('NOAA_', '')
        }

# Create plot with both time series
plt.figure(figsize=(14, 6))
plt.clf()

# Plot both buoys
colors = ['#1f77b4', '#ff7f0e']  # blue and orange
for idx, (buoy_name, data) in enumerate(buoy_data.items()):
    # Convert numpy datetime to matplotlib-compatible format
    datetimes = pd.to_datetime(data['time'])
    plt.plot(datetimes, data['Hs'], label=f"NOAA Buoy {data['buoy_id']}", 
            linewidth=1.5, color=colors[idx], alpha=0.8)

plt.grid(True, alpha=0.3)
plt.xlabel('Year', fontweight='bold', fontsize=12)
plt.ylabel('H$_{s,0}$ [m]', fontweight='bold', fontsize=12)
plt.title('NOAA Buoys 41110 and 41159 - Significant Wave Height', 
         fontweight='bold', fontsize=14)
plt.legend(fontsize=11, loc='upper left')
plt.tight_layout()

# Ensure output directory exists
os.makedirs(plot_dir, exist_ok=True)

# Save figure
output_file = os.path.join(plot_dir, 'NOAA_Buoys_41110_41159_Wave_Height.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_file}")

plt.show()

print("\nPlotting complete!")
