from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_buoy_comparison(df1, df2, buoy_names, stats_dict, plot_dir, running_in_jupyter=True):
    """Create comprehensive comparison plots."""
    valid_idx = ~(df1["Hs"].isna() | df2["Hs"].isna())
    hs1 = df1.loc[valid_idx, "Hs"].values
    hs2 = df2.loc[valid_idx, "Hs"].values
    time1 = df1.loc[valid_idx, "datetime"].values

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle(f"Buoy Comparison: {buoy_names[0]} vs {buoy_names[1]}", fontsize=14, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(time1, hs1, "o-", label=buoy_names[0], linewidth=1.5, markersize=3, alpha=0.7)
    ax.plot(time1, hs2, "s-", label=buoy_names[1], linewidth=1.5, markersize=3, alpha=0.7)
    ax.axhline(2, color="red", linestyle="--", alpha=0.5, label="Storm threshold (2 m)")
    ax.set_xlabel("Date", fontweight="bold")
    ax.set_ylabel("Hs [m]", fontweight="bold")
    ax.set_title("Overlapping Time Series")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.scatter(hs1, hs2, alpha=0.5, s=20)
    x_line = np.array([hs1.min(), hs1.max()])
    y_line = stats_dict["regression_slope"] * x_line + stats_dict["regression_intercept"]
    ax.plot(x_line, y_line, "r--", linewidth=2, label=f"r^2={stats_dict['regression_r2']:.3f}")
    lim = max(hs1.max(), hs2.max())
    ax.plot([0, lim], [0, lim], "k--", alpha=0.3, linewidth=1, label="1:1 line")
    ax.set_xlabel(f"{buoy_names[0]} Hs [m]", fontweight="bold")
    ax.set_ylabel(f"{buoy_names[1]} Hs [m]", fontweight="bold")
    ax.set_title(f"Scatter Plot (r={stats_dict['pearson_r']:.3f}, p<0.001)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")

    ax = axes[1, 0]
    sorted_hs1 = np.sort(hs1)
    sorted_hs2 = np.sort(hs2)
    cdf1 = np.arange(1, len(sorted_hs1) + 1) / len(sorted_hs1)
    cdf2 = np.arange(1, len(sorted_hs2) + 1) / len(sorted_hs2)
    ax.plot(sorted_hs1, cdf1, linewidth=2, label=buoy_names[0])
    ax.plot(sorted_hs2, cdf2, linewidth=2, label=buoy_names[1])
    ax.axvline(2, color="red", linestyle="--", alpha=0.5, label="Storm threshold (2 m)")
    ax.set_xlabel("Hs [m]", fontweight="bold")
    ax.set_ylabel("Cumulative Probability", fontweight="bold")
    ax.set_title(f"CDF Comparison (KS stat={stats_dict['ks_statistic']:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    bins = np.linspace(0, max(hs1.max(), hs2.max()), 30)
    ax.hist(hs1, bins=bins, alpha=0.6, label=buoy_names[0], edgecolor="black")
    ax.hist(hs2, bins=bins, alpha=0.6, label=buoy_names[1], edgecolor="black")
    ax.axvline(2, color="red", linestyle="--", alpha=0.5, linewidth=2, label="Storm threshold")
    ax.set_xlabel("Hs [m]", fontweight="bold")
    ax.set_ylabel("Frequency", fontweight="bold")
    ax.set_title("Distribution Comparison")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    output_file = plot_dir / "Buoy_Comparison_Statistics.png"
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {output_file}")

    if running_in_jupyter:
        plt.show()
    else:
        plt.close()


def plot_wave_power_by_survey(wave_power_df, plot_dir, running_in_jupyter=True):
    """Plot cumulative wave power as side-by-side bars with max Hs overlays."""
    if wave_power_df.empty:
        print("No wave power data to plot")
        return

    buoys = list(wave_power_df["buoy"].unique())
    period_df = (
        wave_power_df[["survey_idx", "start_date", "end_date"]]
        .drop_duplicates()
        .sort_values("survey_idx")
        .reset_index(drop=True)
    )
    master_idx = period_df["survey_idx"].to_numpy()
    x = np.arange(len(master_idx))
    x_labels = [
        f"{row['start_date'].strftime('%Y-%m')}\nto\n{row['end_date'].strftime('%Y-%m')}"
        for _, row in period_df.iterrows()
    ]

    fig, ax = plt.subplots(figsize=(max(16, len(master_idx) * 0.55), 7))
    fig.suptitle(
        "Cumulative Wave Power Between Survey Dates\n(per Splinter et al. 2014)",
        fontsize=14,
        fontweight="bold",
    )

    ax2 = ax.twinx()
    width = 0.8 / max(len(buoys), 1)
    colors = ["steelblue", "darkorange", "seagreen", "slateblue"]

    for i, buoy_name in enumerate(buoys):
        color = colors[i % len(colors)]
        buoy_df = wave_power_df[wave_power_df["buoy"] == buoy_name].set_index("survey_idx").reindex(master_idx)

        power_series = buoy_df["cum_wave_power_MWh_m"]
        maxhs_series = buoy_df["max_Hs"]

        offsets = x - 0.4 + (i + 0.5) * width
        valid_mask = power_series.notna().to_numpy()
        bars = ax.bar(
            offsets[valid_mask],
            power_series.to_numpy()[valid_mask],
            width=width,
            color=color,
            edgecolor="black",
            alpha=0.8,
            label=f"{buoy_name} wave power",
        )

        for bar, val in zip(bars, power_series.to_numpy()[valid_mask]):
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height(),
                f"{val:.1f}",
                ha="center",
                va="bottom",
                fontsize=7,
                rotation=90,
            )

        ax2.plot(
            x,
            maxhs_series.to_numpy(),
            linestyle="--",
            marker="o",
            linewidth=1.8,
            markersize=4,
            color=color,
            alpha=0.9,
            label=f"{buoy_name} max Hs",
        )

    ax.set_ylabel("Wave Power [MWh/m]", fontweight="bold", fontsize=11)
    ax.set_xlabel("Survey Interval", fontweight="bold")
    ax.set_title("Side-by-Side Period-Integrated Wave Power by Buoy", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    ax2.set_ylabel("Max Hs [m]", fontweight="bold", fontsize=11)

    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, loc="upper left", fontsize=8)

    plt.tight_layout()
    output_file = plot_dir / "Cumulative_Wave_Power_Survey_Periods.png"
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {output_file}")

    if running_in_jupyter:
        plt.show()
    else:
        plt.close()
