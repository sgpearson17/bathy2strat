from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def compare_buoys_overlapping(buoy_data):
    """
    Find overlapping time period and align data for both buoys.
    Returns aligned DataFrames for analysis.
    """
    if len(buoy_data) < 2:
        print("Need data from at least 2 buoys")
        return None, None, None

    buoys = list(buoy_data.values())
    buoy_names = list(buoy_data.keys())

    df_buoy1 = pd.DataFrame({"datetime": buoys[0]["time"], "Hs": buoys[0]["Hs"]})
    df_buoy2 = pd.DataFrame({"datetime": buoys[1]["time"], "Hs": buoys[1]["Hs"]})

    overlap_start = max(df_buoy1["datetime"].min(), df_buoy2["datetime"].min())
    overlap_end = min(df_buoy1["datetime"].max(), df_buoy2["datetime"].max())

    df_buoy1_overlap = df_buoy1[
        (df_buoy1["datetime"] >= overlap_start) & (df_buoy1["datetime"] <= overlap_end)
    ].reset_index(drop=True)
    df_buoy2_overlap = df_buoy2[
        (df_buoy2["datetime"] >= overlap_start) & (df_buoy2["datetime"] <= overlap_end)
    ].reset_index(drop=True)

    df_merged = pd.merge(
        df_buoy1_overlap.rename(columns={"Hs": "Hs_buoy1"}),
        df_buoy2_overlap.rename(columns={"Hs": "Hs_buoy2"}),
        on="datetime",
        how="inner",
    )

    print(f"Overlapping period: {overlap_start.date()} to {overlap_end.date()}")
    print(f"  {buoy_names[0]}: {len(df_buoy1_overlap)} records (before alignment)")
    print(f"  {buoy_names[1]}: {len(df_buoy2_overlap)} records (before alignment)")
    print(f"  Aligned data: {len(df_merged)} common time points")

    df_buoy1_overlap = df_merged[["datetime", "Hs_buoy1"]].rename(columns={"Hs_buoy1": "Hs"})
    df_buoy2_overlap = df_merged[["datetime", "Hs_buoy2"]].rename(columns={"Hs_buoy2": "Hs"})

    return df_buoy1_overlap, df_buoy2_overlap, buoy_names


def calculate_buoy_statistics(df1, df2, buoy_names):
    """Calculate comprehensive statistical comparison between buoys."""
    valid_idx = ~(df1["Hs"].isna() | df2["Hs"].isna())
    hs1_clean = df1.loc[valid_idx, "Hs"].values
    hs2_clean = df2.loc[valid_idx, "Hs"].values

    stats_dict = {}

    corr_pearson, pval_pearson = stats.pearsonr(hs1_clean, hs2_clean)
    stats_dict["pearson_r"] = corr_pearson
    stats_dict["pearson_pval"] = pval_pearson

    corr_spearman, pval_spearman = stats.spearmanr(hs1_clean, hs2_clean)
    stats_dict["spearman_r"] = corr_spearman
    stats_dict["spearman_pval"] = pval_spearman

    slope, intercept, r_value, p_value, _std_err = stats.linregress(hs1_clean, hs2_clean)
    stats_dict["regression_slope"] = slope
    stats_dict["regression_intercept"] = intercept
    stats_dict["regression_r2"] = r_value**2
    stats_dict["regression_pval"] = p_value

    stats_dict["buoy1_mean"] = np.mean(hs1_clean)
    stats_dict["buoy1_std"] = np.std(hs1_clean)
    stats_dict["buoy1_median"] = np.median(hs1_clean)
    stats_dict["buoy1_max"] = np.max(hs1_clean)

    stats_dict["buoy2_mean"] = np.mean(hs2_clean)
    stats_dict["buoy2_std"] = np.std(hs2_clean)
    stats_dict["buoy2_median"] = np.median(hs2_clean)
    stats_dict["buoy2_max"] = np.max(hs2_clean)

    ks_stat, ks_pval = stats.ks_2samp(hs1_clean, hs2_clean)
    stats_dict["ks_statistic"] = ks_stat
    stats_dict["ks_pval"] = ks_pval

    storm_threshold = 2.0
    storm_count_1 = np.sum(hs1_clean > storm_threshold)
    storm_count_2 = np.sum(hs2_clean > storm_threshold)
    stats_dict["storms_above_2m_buoy1"] = storm_count_1
    stats_dict["storms_above_2m_buoy2"] = storm_count_2
    stats_dict["storms_above_2m_pct_buoy1"] = 100 * storm_count_1 / len(hs1_clean)
    stats_dict["storms_above_2m_pct_buoy2"] = 100 * storm_count_2 / len(hs2_clean)

    return stats_dict


def print_comparison_statistics(stats_dict, buoy_names):
    """Print formatted statistical comparison."""
    print("\n" + "=" * 70)
    print(f"BUOY COMPARISON STATISTICS: {buoy_names[0]} vs {buoy_names[1]}")
    print("=" * 70)

    print("\n--- CORRELATION ---")
    print(f"Pearson r: {stats_dict['pearson_r']:.4f} (p-value: {stats_dict['pearson_pval']:.2e})")
    print(f"Spearman rho: {stats_dict['spearman_r']:.4f} (p-value: {stats_dict['spearman_pval']:.2e})")

    print("\n--- LINEAR REGRESSION: Hs(41159) = a*Hs(41110) + b ---")
    print(f"Slope (a): {stats_dict['regression_slope']:.4f}")
    print(f"Intercept (b): {stats_dict['regression_intercept']:.4f}")
    print(f"R^2: {stats_dict['regression_r2']:.4f} (p-value: {stats_dict['regression_pval']:.2e})")

    print("\n--- BASIC STATISTICS ---")
    print(f"{buoy_names[0]:20s} | {buoy_names[1]:20s}")
    print(f"Mean Hs: {stats_dict['buoy1_mean']:6.2f} m    | {stats_dict['buoy2_mean']:6.2f} m")
    print(f"Std Dev: {stats_dict['buoy1_std']:6.2f} m    | {stats_dict['buoy2_std']:6.2f} m")
    print(f"Median:  {stats_dict['buoy1_median']:6.2f} m    | {stats_dict['buoy2_median']:6.2f} m")
    print(f"Max:     {stats_dict['buoy1_max']:6.2f} m    | {stats_dict['buoy2_max']:6.2f} m")

    print("\n--- DISTRIBUTION TEST (Kolmogorov-Smirnov) ---")
    print(f"KS statistic: {stats_dict['ks_statistic']:.4f} (p-value: {stats_dict['ks_pval']:.2e})")
    if stats_dict["ks_pval"] > 0.05:
        print("  -> Distributions are NOT significantly different (fail to reject H0)")
    else:
        print("  -> Distributions ARE significantly different (reject H0)")

    print("\n--- STORM EXCEEDANCE (Hs > 2m, per Splinter et al.) ---")
    print(
        f"{buoy_names[0]}: {stats_dict['storms_above_2m_buoy1']:.0f} events "
        f"({stats_dict['storms_above_2m_pct_buoy1']:.1f}%)"
    )
    print(
        f"{buoy_names[1]}: {stats_dict['storms_above_2m_buoy2']:.0f} events "
        f"({stats_dict['storms_above_2m_pct_buoy2']:.1f}%)"
    )
    print("=" * 70 + "\n")
