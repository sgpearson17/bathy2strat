from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


def validate_directories(data_dirs, plot_dir):
    """Validate and normalize input directories."""
    data_dirs = [Path(d).resolve() for d in data_dirs]
    plot_dir = Path(plot_dir).resolve()

    for data_dir in data_dirs:
        if not data_dir.exists():
            raise FileNotFoundError(f"Data directory does not exist: {data_dir}")
        if not data_dir.is_dir():
            raise NotADirectoryError(f"Path is not a directory: {data_dir}")

    plot_dir.mkdir(parents=True, exist_ok=True)
    return data_dirs, plot_dir


def setup_plot_style():
    """Set up matplotlib with a readable bold italic font style."""
    available_fonts = set(matplotlib.font_manager.get_font_names())

    if 'Arial' in available_fonts or 'DejaVu Sans' in available_fonts:
        if 'Arial' in available_fonts:
            font_name = 'Arial'
        else:
            font_name = 'DejaVu Sans'
            print("Warning: Arial font not found. Using DejaVu Sans as fallback.")
    else:
        font_name = 'sans-serif'
        print("Warning: Arial and DejaVu fonts not found. Using default sans-serif.")

    plt.rcParams['font.family'] = font_name
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['font.style'] = 'italic'
    plt.rcParams['font.size'] = 10

    print(f"Figure style set to: {font_name}, Bold, Italic")


def load_and_plot_wave_data(data_dirs, plot_dir, show_plot=True, running_in_jupyter=False):
    """
    Load NOAA meteorological data from multiple directories and plot wave heights.
    """
    column_names = ['YY', 'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST', 'WVHT',
                    'DPD', 'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']

    buoy_data = {}

    for data_dir in data_dirs:
        data_dir = Path(data_dir)
        buoy_id = data_dir.name
        txt_files = sorted(data_dir.glob('*.txt'))
        txt_files = [f for f in txt_files if 'README' not in f.name.upper()]

        print(f"Processing {buoy_id}: Found {len(txt_files)} data files")

        dfs = []
        for file_path in txt_files:
            try:
                df = pd.read_csv(
                    file_path,
                    sep=r'\s+',
                    comment='#',
                    names=column_names,
                    dtype={'WVHT': float},
                    na_values=['99.0', '99', '99.00'],
                )
                df['datetime'] = pd.to_datetime(
                    df[['YY', 'MM', 'DD', 'hh', 'mm']].rename(
                        columns={'YY': 'year', 'MM': 'month', 'DD': 'day', 'hh': 'hour', 'mm': 'minute'}
                    )
                )
                dfs.append(df)
                print(f"  Loaded: {file_path.name}")
            except Exception as e:
                print(f"  Error reading {file_path.name}: {e}")

        if dfs:
            meteo_data = pd.concat(dfs, ignore_index=True)
            meteo_data = meteo_data.sort_values('datetime').reset_index(drop=True)
            buoy_data[buoy_id] = {
                'time': meteo_data['datetime'].values,
                'Hs': meteo_data['WVHT'].values,
                'buoy_id': buoy_id.replace('NOAA_', '')
            }

    fig = plt.figure(figsize=(14, 6))
    colors = ['#1f77b4', '#ff7f0e']
    for idx, (buoy_name, data) in enumerate(buoy_data.items()):
        datetimes = pd.to_datetime(data['time'])
        plt.plot(datetimes, data['Hs'], label=f"NOAA Buoy {data['buoy_id']}",
                 linewidth=1.5, color=colors[idx], alpha=0.8)

    plt.grid(True, alpha=0.3)
    plt.xlabel('Year', fontweight='bold', fontsize=12)
    plt.ylabel('H$_{s,0}$ [m]', fontweight='bold', fontsize=12)
    plt.title('NOAA Buoys 41110 and 41159 - Significant Wave Height', fontweight='bold', fontsize=14)
    plt.legend(fontsize=11, loc='upper left')
    plt.tight_layout()

    output_file = plot_dir / 'NOAA_Buoys_41110_41159_Wave_Height.png'
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {output_file}")

    if show_plot:
        plt.show()
    else:
        plt.close()

    return buoy_data, fig
