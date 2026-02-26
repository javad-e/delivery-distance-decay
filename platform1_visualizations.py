import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from scipy.interpolate import UnivariateSpline
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d

plt.style.use("seaborn-v0_8-darkgrid")
sns.set_palette("husl")
plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.size"] = 11

COLORS = {"no_incentive": "#2E86AB", "with_incentive": "#A23B72"}
LABELS = {"no_incentive": "No Incentive", "with_incentive": "With Incentive"}

EXPORT_DIR = "./grouped_exports/"


# ---------------------------------------------------------------------------
# Density plot (population-level)
# ---------------------------------------------------------------------------

def plot_density(df: pd.DataFrame, xlim=None):
    fig, ax = plt.subplots(figsize=(11, 6))

    for cat in ["no_incentive", "with_incentive"]:
        data = df[df["incentive_cat"] == cat]["distance_km"].values
        if len(data) == 0:
            continue

        kde = gaussian_kde(data, bw_method=0.2)
        x = np.linspace(data.min(), data.max(), 1000)
        spline = UnivariateSpline(x, kde(x), s=0.001, k=3)
        xs = np.linspace(data.min(), data.max(), 2000)

        ax.plot(xs, spline(xs), linewidth=2.5, label=LABELS[cat], color=COLORS[cat], alpha=0.8)
        ax.axvline(data.mean(), color=COLORS[cat], linestyle="--", linewidth=1.5, alpha=0.6)

    ax.set_xlabel("Delivery Distance (km)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Density", fontsize=13, fontweight="bold")
    ax.set_title("Distribution of Delivery Distances by Incentive Status", fontsize=15, fontweight="bold", pad=20)
    if xlim:
        ax.set_xlim(xlim)
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# Within-customer density (demeaned)
# ---------------------------------------------------------------------------

def plot_within_customer_density(df: pd.DataFrame, xlim=None, curve_output_path: str = EXPORT_DIR + "within_customer_curves.csv"):
    df_plot = df.copy()
    df_plot["distance_demeaned"] = df_plot.groupby("account_id")["distance_km"].transform(lambda x: x - x.mean())

    fig, ax = plt.subplots(figsize=(11, 6))
    curve_list = []

    for cat in ["no_incentive", "with_incentive"]:
        data = df_plot[df_plot["incentive_cat"] == cat]["distance_demeaned"].values
        if len(data) == 0:
            continue

        kde = gaussian_kde(data, bw_method=0.2)
        x = np.linspace(data.min(), data.max(), 1000)
        spline = UnivariateSpline(x, kde(x), s=0.001, k=3)
        xs = np.linspace(data.min(), data.max(), 2000)
        ys = spline(xs)

        curve_list.append(pd.DataFrame({"incentive_cat": cat, "x": xs, "y": ys}))
        ax.plot(xs, ys, linewidth=2.5, label=LABELS[cat], color=COLORS[cat], alpha=0.8)
        ax.axvline(data.mean(), color=COLORS[cat], linestyle="--", linewidth=1.5, alpha=0.6)

    if curve_list:
        pd.concat(curve_list, ignore_index=True).to_csv(curve_output_path, index=False)

    ax.axvline(0, color="black", linestyle=":", linewidth=1.5, alpha=0.5)
    ax.set_xlabel("Distance Deviation from Customer Mean (km)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Density", fontsize=13, fontweight="bold")
    ax.set_title("Distribution of Within-Customer Distance Deviations", fontsize=15, fontweight="bold", pad=20)
    if xlim:
        ax.set_xlim(xlim)
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# CDF plot
# ---------------------------------------------------------------------------

def plot_cdf(df: pd.DataFrame, xlim=None):
    fig, ax = plt.subplots(figsize=(11, 6))

    for cat in ["no_incentive", "with_incentive"]:
        data = np.sort(df[df["incentive_cat"] == cat]["distance_km"].values)
        if len(data) == 0:
            continue

        n_points = min(10000, len(data))
        indices = np.linspace(0, len(data) - 1, n_points, dtype=int)
        x_plot = data[indices]
        y_smooth = gaussian_filter1d((indices + 1) / len(data), sigma=2)

        ax.plot(x_plot, y_smooth, linewidth=2.5, label=LABELS[cat], color=COLORS[cat], alpha=0.9)

    for p in [0.25, 0.5, 0.75]:
        ax.axhline(y=p, color="gray", linestyle=":", linewidth=1, alpha=0.5)

    ax.set_xlabel("Delivery Distance (km)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Cumulative Probability", fontsize=13, fontweight="bold")
    ax.set_title("Cumulative Distribution of Delivery Distances", fontsize=15, fontweight="bold", pad=20)
    if xlim:
        ax.set_xlim(xlim)
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# Side-by-side bar chart
# ---------------------------------------------------------------------------

def plot_grouped_bars(
    df_means: pd.DataFrame,
    title: str = "",
    output_path: str = None,
    highlight_groups=None,
    comparison_text: str = None,
):
    fig, ax = plt.subplots(figsize=(12, 6))
    groups = df_means["group"].values
    x = np.arange(len(groups))
    width = 0.35

    bars1 = ax.bar(x - width / 2, df_means["mean_no_incentive"], width, label="No Incentive", color=COLORS["no_incentive"], alpha=0.8, edgecolor="black", linewidth=1.5)
    bars2 = ax.bar(x + width / 2, df_means["mean_with_incentive"], width, label="With Incentive", color=COLORS["with_incentive"], alpha=0.8, edgecolor="black", linewidth=1.5)

    for bars in [bars1, bars2]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height(), f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    if highlight_groups:
        for i, g in enumerate(groups):
            if g in highlight_groups:
                rect = plt.Rectangle((i - 0.45, ax.get_ylim()[0]), 0.9, ax.get_ylim()[1] - ax.get_ylim()[0], facecolor="yellow", alpha=0.2, zorder=0)
                ax.add_patch(rect)

    if comparison_text:
        ax.text(0.98, 0.98, comparison_text, transform=ax.transAxes, fontsize=10, va="top", ha="right", bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9, edgecolor="black"))

    ax.set_xlabel("Customer Group", fontsize=13, fontweight="bold")
    ax.set_ylabel("Mean Distance (km)", fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=15, fontweight="bold", pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.legend(loc="upper left", framealpha=0.95)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()

    if output_path:
        df_means.to_csv(output_path, index=False)

    return fig, ax


# ---------------------------------------------------------------------------
# Split violin plot
# ---------------------------------------------------------------------------

def _prepare_violin_data(df: pd.DataFrame, group_col: str, group_values: list) -> pd.DataFrame:
    rows = []
    for val in group_values:
        sub = df[df[group_col] == val]
        for cat in ["no_incentive", "with_incentive"]:
            for dist in sub[sub["incentive_cat"] == cat]["distance_km"].values:
                rows.append({"group": val, "distance_km": dist, "incentive": cat})
    return pd.DataFrame(rows)


def plot_split_violin(
    df: pd.DataFrame,
    group_col: str,
    group_values: list,
    title: str = "",
    output_path: str = None,
    highlight_groups=None,
    comparison_text: str = None,
    ylim=None,
):
    df_plot = _prepare_violin_data(df, group_col, group_values)

    if output_path:
        df_plot.to_csv(output_path, index=False)

    fig, ax = plt.subplots(figsize=(12, 7))

    sns.violinplot(
        data=df_plot, x="group", y="distance_km", hue="incentive",
        split=True, inner="quartile",
        palette={"no_incentive": COLORS["no_incentive"], "with_incentive": COLORS["with_incentive"]},
        ax=ax, alpha=0.8, cut=0,
    )

    upper = ylim[1] if ylim else df_plot["distance_km"].quantile(0.95) * 1.1
    ax.set_ylim(ylim or (0, upper))

    if highlight_groups:
        for i, g in enumerate(group_values):
            if g in highlight_groups:
                rect = Rectangle((i - 0.4, ax.get_ylim()[0]), 0.8, ax.get_ylim()[1] - ax.get_ylim()[0], facecolor="yellow", alpha=0.15, zorder=0)
                ax.add_patch(rect)

    if comparison_text:
        ax.text(0.98, 0.98, comparison_text, transform=ax.transAxes, fontsize=10, va="top", ha="right", bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.9, edgecolor="black"))

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ["No Incentive", "With Incentive"], loc="upper left", framealpha=0.95)
    ax.set_xlabel("Customer Group", fontsize=13, fontweight="bold")
    ax.set_ylabel("Delivery Distance (km)", fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=15, fontweight="bold", pad=20)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    return fig, ax


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from analyze_incentives import (
        prepare_data, run_heterogeneity_by_income, run_heterogeneity_by_distance,
        test_coefficient_difference, format_comparison_text,
    )
    from preprocessing import load_income_geodata, load_orders, load_incentives, merge_orders_incentives, engineer_features
    import os

    os.makedirs(EXPORT_DIR, exist_ok=True)

    gdf_webmerc = load_income_geodata()
    orders = load_orders()
    incentives = load_incentives()
    dubai_df = engineer_features(merge_orders_incentives(orders, incentives))
    df_food = dubai_df[dubai_df["vertical_class"] == "food"].drop_duplicates(subset="order_id").reset_index(drop=True)
    df_analysis = prepare_data(df_food, gdf_webmerc)

    xlim = (0, df_analysis["distance_km"].quantile(0.95))
    plot_density(df_analysis, xlim)
    plot_cdf(df_analysis, xlim)

    xlim_within = (-df_analysis["distance_km"].std() * 2, df_analysis["distance_km"].std() * 2)
    plot_within_customer_density(df_analysis, xlim=xlim_within)

    df_income_results, df_income_means = run_heterogeneity_by_income(df_analysis)
    df_distance_results, df_distance_means = run_heterogeneity_by_distance(df_analysis)

    income_test = test_coefficient_difference(df_income_results, "Income 1")
    distance_test = test_coefficient_difference(df_distance_results, "Q4")

    plot_grouped_bars(df_income_means, title="Mean Delivery Distance by Income Category", output_path=EXPORT_DIR + "income_bar_data.csv", highlight_groups=["Income 1"], comparison_text=format_comparison_text(income_test))
    plot_grouped_bars(df_distance_means, title="Mean Delivery Distance by Distance Quartile", output_path=EXPORT_DIR + "distance_bar_data.csv", highlight_groups=["Q4"], comparison_text=format_comparison_text(distance_test))

    income_cats = sorted(df_analysis["income_category"].dropna().unique())
    plot_split_violin(df_analysis, "income_category", income_cats, title="Delivery Distances by Income Category", output_path=EXPORT_DIR + "income_violin_data.csv", highlight_groups=[1], comparison_text=format_comparison_text(income_test))
    plot_split_violin(df_analysis, "customer_dist_quartile", ["Q1", "Q2", "Q3", "Q4"], title="Delivery Distances by Distance Quartile", output_path=EXPORT_DIR + "distance_violin_data.csv", highlight_groups=["Q4"], comparison_text=format_comparison_text(distance_test))

    plt.show()
