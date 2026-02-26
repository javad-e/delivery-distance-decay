import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.stats import norm
from linearmodels.panel import PanelOLS

warnings.filterwarnings("ignore")

OUTPUT_PATH = "~/imperial/Reach of Last Mile/"


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def haversine_distance(lat1, lon1, lat2, lon2) -> np.ndarray:
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 6371 * 2 * np.arcsin(np.sqrt(a))


def assign_income_category(df: pd.DataFrame, gdf_income: gpd.GeoDataFrame) -> pd.DataFrame:
    gdf_customers = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.customer_lng, df.customer_lat),
        crs="EPSG:4326",
    ).to_crs("EPSG:3857")

    merged = gpd.sjoin_nearest(
        gdf_customers,
        gdf_income[["income_category", "geometry"]],
        how="left",
        distance_col="dist_to_polygon",
    )
    merged["income_category"] = merged["income_category"].fillna(
        merged.groupby("account_id")["income_category"].transform(
            lambda x: x.mode()[0] if x.notna().any() else np.nan
        )
    )
    return merged.drop(columns=["geometry", "index_right", "dist_to_polygon"]).reset_index(drop=True)


def prepare_data(df: pd.DataFrame, gdf_income: gpd.GeoDataFrame, min_orders: int = 3) -> pd.DataFrame:
    print("Preparing data...")
    df = df.copy()

    df["distance_km"] = haversine_distance(
        df["vendor_lat"], df["vendor_lng"],
        df["customer_lat"], df["customer_lng"],
    )
    df = df[df["distance_km"] > 0].copy()

    df["incentive_cat"] = df["delivery_incentive_lc"].apply(
        lambda x: "with_incentive" if x > 0 else "no_incentive"
    )

    df = assign_income_category(df, gdf_income)

    df["fee_per_km"] = df["delivery_fee_amount_lc"] / df["distance_km"]
    df["avg_fee_per_km"] = df.groupby(["vendor_id", "time_period"])["fee_per_km"].transform("mean")

    order_counts = df.groupby("account_id").size()
    df = df[df["account_id"].isin(order_counts[order_counts >= min_orders].index)].copy()

    has_variation = df.groupby("account_id")["incentive_cat"].apply(lambda x: x.nunique()) > 1
    df = df[df["account_id"].isin(has_variation[has_variation].index)].copy()

    customer_avg_dist = df.groupby("account_id")["distance_km"].mean()
    df["customer_dist_quartile"] = df["account_id"].map(
        pd.qcut(customer_avg_dist, q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    )

    df["log_distance"] = np.log(df["distance_km"])
    df["log_basket"] = np.log(df["basket_amount_lc"] + 1)
    df["is_food"] = (df["vertical_class"] == "food").astype(int)
    df["is_weekend"] = (df["day_type"] == "weekend").astype(int)
    df["has_incentive"] = (df["incentive_cat"] == "with_incentive").astype(int)

    print(f"Final sample: {len(df):,} orders, {df['account_id'].nunique():,} customers")
    return df


# ---------------------------------------------------------------------------
# Estimation
# ---------------------------------------------------------------------------

def _export_results(results, model_name: str, output_path: str = OUTPUT_PATH):
    coef_data = []
    for var in results.params.index:
        pval = results.pvalues[var]
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        coef_data.append({
            "variable": var,
            "coefficient": results.params[var],
            "std_error": results.std_errors[var],
            "t_statistic": results.tstats[var],
            "p_value": pval,
            "significance": sig,
        })

    stats_data = {
        "model_name": [model_name],
        "N_observations": [results.nobs],
        "N_entities": [results.entity_info["total"]],
        "R_squared": [results.rsquared],
        "R_squared_within": [results.rsquared_within],
        "R_squared_between": [results.rsquared_between],
        "R_squared_overall": [results.rsquared_overall],
        "F_statistic": [results.f_statistic.stat if hasattr(results.f_statistic, "stat") else np.nan],
        "F_pvalue": [results.f_statistic.pval if hasattr(results.f_statistic, "pval") else np.nan],
    }

    pd.DataFrame(coef_data).to_csv(output_path + "regression_coefficients.csv", index=False)
    pd.DataFrame(stats_data).to_csv(output_path + "regression_stats.csv", index=False)


def estimate_model(df: pd.DataFrame, controls: bool = True, subset_filter=None, subset_name: str = "Full Sample"):
    df_model = df if subset_filter is None else df[subset_filter].copy()

    formula = "log_distance ~ has_incentive"
    if controls:
        formula += " + log_basket + avg_fee_per_km + is_food + is_weekend + C(time_period) + C(quarter)"
    formula += " + EntityEffects"

    df_panel = df_model.set_index(["account_id", df_model.groupby("account_id").cumcount()])
    results = PanelOLS.from_formula(formula, data=df_panel).fit(cov_type="clustered", cluster_entity=True)

    coef = results.params.get("has_incentive", np.nan)
    se = results.std_errors.get("has_incentive", np.nan)
    pval = results.pvalues.get("has_incentive", np.nan)

    summary = {
        "subset": subset_name,
        "n_obs": results.nobs,
        "n_customers": results.entity_info["total"],
        "coef": coef,
        "se": se,
        "pval": pval,
        "pct_change": (np.exp(coef) - 1) * 100 if not np.isnan(coef) else np.nan,
        "r2": results.rsquared,
        "r2_within": results.rsquared_within,
    }

    _export_results(results, subset_name)
    return results, summary


def run_baseline_analysis(df: pd.DataFrame):
    print("\nBaseline estimation...")
    results, summary = estimate_model(df, controls=True)
    print(f"N = {summary['n_obs']:,.0f} orders, {summary['n_customers']:,.0f} customers")
    print(f"Incentive effect: {summary['coef']:.4f} (SE: {summary['se']:.4f}), p={summary['pval']:.4f}")
    print(f"  → {summary['pct_change']:+.2f}% change in distance")
    print(f"R² (within): {summary['r2_within']:.4f}")
    return results, summary


# ---------------------------------------------------------------------------
# Heterogeneity analysis
# ---------------------------------------------------------------------------

def _group_means(df: pd.DataFrame, group_col: str, group_values: list) -> pd.DataFrame:
    rows = []
    for val in group_values:
        sub = df[df[group_col] == val]
        rows.append({
            "group": val,
            "mean_no_incentive": sub[sub["incentive_cat"] == "no_incentive"]["distance_km"].mean(),
            "mean_with_incentive": sub[sub["incentive_cat"] == "with_incentive"]["distance_km"].mean(),
            "n_obs": len(sub),
        })
    return pd.DataFrame(rows)


def run_heterogeneity_by_income(df: pd.DataFrame):
    print("\nHeterogeneity by income category...")
    income_cats = sorted(df["income_category"].dropna().unique())
    results_list = []
    for cat in income_cats:
        _, summary = estimate_model(df, controls=True, subset_filter=df["income_category"] == cat, subset_name=f"Income {cat}")
        results_list.append(summary)
        print(f"  Income {cat}: N={summary['n_obs']:,.0f}, effect={summary['pct_change']:+.2f}%")

    df_results = pd.DataFrame(results_list)
    df_means = _group_means(df, "income_category", income_cats)
    df_means["group"] = "Income " + df_means["group"].astype(str)
    return df_results, df_means


def run_heterogeneity_by_distance(df: pd.DataFrame):
    print("\nHeterogeneity by distance quartile...")
    quartiles = ["Q1", "Q2", "Q3", "Q4"]
    results_list = []
    for q in quartiles:
        _, summary = estimate_model(df, controls=True, subset_filter=df["customer_dist_quartile"] == q, subset_name=q)
        results_list.append(summary)
        print(f"  {q}: N={summary['n_obs']:,.0f}, effect={summary['pct_change']:+.2f}%")

    df_results = pd.DataFrame(results_list)
    df_means = _group_means(df, "customer_dist_quartile", quartiles)
    return df_results, df_means


def test_coefficient_difference(df_results: pd.DataFrame, target_group: str, group_col: str = "subset") -> dict:
    target = df_results[df_results[group_col] == target_group].iloc[0]
    others = df_results[df_results[group_col] != target_group]

    weights = 1 / (others["se"] ** 2)
    coef_others = np.average(others["coef"], weights=weights)
    se_others = np.sqrt(1 / weights.sum())

    diff = target["coef"] - coef_others
    t_stat = diff / np.sqrt(target["se"] ** 2 + se_others ** 2)
    p_val = 2 * (1 - norm.cdf(abs(t_stat)))

    return {
        "target_group": target_group,
        "target_pct": (np.exp(target["coef"]) - 1) * 100,
        "others_pct": (np.exp(coef_others) - 1) * 100,
        "difference": ((np.exp(target["coef"]) - 1) - (np.exp(coef_others) - 1)) * 100,
        "t_stat": t_stat,
        "p_value": p_val,
        "significant": p_val < 0.05,
    }


def format_comparison_text(test_result: dict) -> str:
    sig = "***" if test_result["p_value"] < 0.01 else "**" if test_result["p_value"] < 0.05 else "*" if test_result["p_value"] < 0.10 else ""
    g = test_result["target_group"]
    return (
        f"{g} vs Others:\n"
        f"  {g}: {test_result['target_pct']:+.2f}%\n"
        f"  Others: {test_result['others_pct']:+.2f}%\n"
        f"  Difference: {test_result['difference']:+.2f}%\n"
        f"  p-value: {test_result['p_value']:.3f} {sig}"
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from preprocessing import load_income_geodata, load_orders, load_incentives, merge_orders_incentives, engineer_features

    gdf_webmerc = load_income_geodata()
    orders = load_orders()
    incentives = load_incentives()
    dubai_df = engineer_features(merge_orders_incentives(orders, incentives))

    df_food = dubai_df[dubai_df["vertical_class"] == "food"].drop_duplicates(subset="order_id").reset_index(drop=True)
    df_analysis = prepare_data(df_food, gdf_webmerc)

    results_baseline, summary_baseline = run_baseline_analysis(df_analysis)

    df_income_results, df_income_means = run_heterogeneity_by_income(df_analysis)
    df_distance_results, df_distance_means = run_heterogeneity_by_distance(df_analysis)

    income_test = test_coefficient_difference(df_income_results, "Income 1")
    distance_test = test_coefficient_difference(df_distance_results, "Q4")

    print("\nIncome 1 vs others:", format_comparison_text(income_test))
    print("\nQ4 vs others:", format_comparison_text(distance_test))
