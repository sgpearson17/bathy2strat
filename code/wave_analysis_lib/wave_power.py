from __future__ import annotations

import numpy as np
import pandas as pd


def load_survey_dates(survey_file_path):
    """Load survey dates from text file (one date per line, ISO format)."""
    dates = []
    with open(survey_file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    date = pd.to_datetime(line)
                    dates.append(date)
                except Exception:
                    print(f"Warning: Could not parse date: {line}")
    return sorted(dates)


def reconstruct_buoy_dataframe(buoy_data_dict):
    """
    Reconstruct DataFrame from buoy_data dictionary with wave period estimation.
    Need Tp for wave power calculation per Splinter et al.
    """
    dfs = {}
    for buoy_name, data in buoy_data_dict.items():
        df = pd.DataFrame({"datetime": data["time"], "Hs": data["Hs"]})
        if "Tp" in data:
            df["Tp"] = data["Tp"]
        else:
            # Estimate Tp from Hs using a simple empirical relationship.
            df["Tp"] = 3.5 + 0.6 * df["Hs"]
        dfs[buoy_name] = df
    return dfs


def calculate_cumulative_wave_power(df, date_start, date_end, date_label=None):
    """
    Calculate cumulative wave power per Splinter et al. (2014).

    Equation 3: sum(P) = integral( (rho*g^2/64*pi) * Hs^2 * Tp * dt )
    """
    rho = 1025
    g = 9.81

    mask = (df["datetime"] >= date_start) & (df["datetime"] <= date_end)
    df_period = df[mask].copy()

    if len(df_period) == 0:
        return None, 0, {"mean_Hs": np.nan, "max_Hs": np.nan, "n_records": 0}

    df_period = df_period.dropna(subset=["Hs", "Tp"])
    df_period["time_diff"] = df_period["datetime"].diff().dt.total_seconds() / 3600
    df_period["power"] = (rho * g**2) / (64 * np.pi) * df_period["Hs"] ** 2 * df_period["Tp"] * df_period["time_diff"]

    cum_power = df_period["power"].sum() / 1e6

    hs_stats = {
        "mean_Hs": df_period["Hs"].mean(),
        "max_Hs": df_period["Hs"].max(),
        "n_records": len(df_period),
        "n_days": (date_end - date_start).days,
    }

    return cum_power, len(df_period), hs_stats


def calculate_wave_power_between_surveys(buoy_dfs, survey_dates, buoy_names=None):
    """Calculate cumulative wave power between consecutive survey dates."""
    results = []

    for i in range(len(survey_dates) - 1):
        start_date = survey_dates[i]
        end_date = survey_dates[i + 1]

        for buoy_name, df in buoy_dfs.items():
            cum_power, _n_records, hs_stats = calculate_cumulative_wave_power(
                df, start_date, end_date, f"{buoy_name}_{i}"
            )

            if cum_power is not None:
                results.append(
                    {
                        "buoy": buoy_name,
                        "survey_idx": i,
                        "start_date": start_date,
                        "end_date": end_date,
                        "period_days": hs_stats["n_days"],
                        "cum_wave_power_MWh_m": cum_power,
                        "mean_Hs": hs_stats["mean_Hs"],
                        "max_Hs": hs_stats["max_Hs"],
                        "n_records": hs_stats["n_records"],
                    }
                )

    return pd.DataFrame(results)


def save_wave_power_results(wave_power_df, output_path):
    """Save cumulative wave power results in tab-delimited text format."""
    if wave_power_df.empty:
        print("No wave power results to save.")
        return

    output_df = wave_power_df.copy().sort_values(["survey_idx", "buoy"])
    output_df["start_date"] = pd.to_datetime(output_df["start_date"]).dt.strftime("%Y-%m-%d")
    output_df["end_date"] = pd.to_datetime(output_df["end_date"]).dt.strftime("%Y-%m-%d")

    output_df.to_csv(output_path, sep="\t", index=False, float_format="%.6f")
    print(f"Wave power statistics saved to: {output_path}")
