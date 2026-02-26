import os
import pickle
import numpy as np
import pandas as pd
import geopandas as gpd
import h3
import osmnx as ox
import networkx as nx
import requests
from math import radians, cos, sin, asin, sqrt
from shapely.geometry import Point, Polygon, box
from tqdm import tqdm
from time import sleep

RESOLUTION = 7
OSRM_URL = "http://router.project-osrm.org"
DISTANCE_CACHE_FILE = "distance_cache.pkl"

DUBAI_BBOX = dict(min_lon=55.0, min_lat=24.9, max_lon=55.45, max_lat=25.4, buffer=0.1)


# ---------------------------------------------------------------------------
# Distance utilities
# ---------------------------------------------------------------------------

def haversine(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    return 6371 * 2 * asin(sqrt(a))


def haversine_vectorized(lon1, lat1, lon2, lat2) -> np.ndarray:
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 6371 * 2 * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------------------
# Basic metrics
# ---------------------------------------------------------------------------

def compute_basic_metrics(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["picked_up_at"] = pd.to_datetime(df["picked_up_at"], utc=True, errors="coerce")
    df["delivered_at"] = pd.to_datetime(df["delivered_at"], utc=True, errors="coerce")
    df["travel_time"] = (df["delivered_at"] - df["picked_up_at"]).dt.total_seconds() / 60.0
    df["euclidean_km"] = haversine_vectorized(
        df["vendor_lng"], df["vendor_lat"],
        df["customer_lng"], df["customer_lat"],
    )
    df["min_to_km_ratio"] = df["travel_time"] / df["euclidean_km"]
    df["perc_discounted"] = df["delivery_incentive_lc"].fillna(0) / df["delivery_fee_amount_lc"]
    df["fee_per_minute"] = df["delivery_fee_amount_lc"] / df["travel_time"]
    df["fee_per_km"] = df["delivery_fee_amount_lc"] / df["euclidean_km"]
    return df


# ---------------------------------------------------------------------------
# Bounding box + H3 grid
# ---------------------------------------------------------------------------

def create_bbox(min_lon=55.0, min_lat=24.9, max_lon=55.45, max_lat=25.4, buffer=0.1):
    return box(min_lon - buffer, min_lat - buffer, max_lon + buffer, max_lat + buffer)


def generate_h3_cells(bbox, h3_res: int) -> gpd.GeoDataFrame:
    outer_ring = [[y, x] for x, y in bbox.exterior.coords[:]]
    latlng_poly = h3.LatLngPoly(outer_ring)
    all_cells = set(h3.polygon_to_cells(latlng_poly, res=h3_res))

    filtered, polygons = [], []
    for cell in all_cells:
        poly_coords = [(lng, lat) for lat, lng in h3.cell_to_boundary(cell)]
        poly = Polygon(poly_coords)
        if bbox.contains(poly.centroid):
            filtered.append(cell)
            polygons.append(poly)

    return gpd.GeoDataFrame({"h3_index": filtered, "geometry": polygons}, crs="EPSG:4326")


# ---------------------------------------------------------------------------
# Trip filtering + H3 assignment
# ---------------------------------------------------------------------------

def filter_trips_to_bbox(df: pd.DataFrame, bbox, input_crs="EPSG:4326") -> pd.DataFrame:
    vendor_pts = gpd.GeoSeries(
        [Point(xy) for xy in zip(df["vendor_lng"], df["vendor_lat"])], crs=input_crs
    )
    customer_pts = gpd.GeoSeries(
        [Point(xy) for xy in zip(df["customer_lng"], df["customer_lat"])], crs=input_crs
    )
    vendor_inside = vendor_pts.within(bbox).to_numpy()
    customer_inside = customer_pts.within(bbox).to_numpy()
    return df[vendor_inside & customer_inside].reset_index(drop=True)


def filter_income_polygons(gdf_webmerc: gpd.GeoDataFrame, bbox) -> gpd.GeoDataFrame:
    gdf = gdf_webmerc.to_crs("EPSG:4326")
    return gdf[gdf.intersects(bbox)].copy()


def assign_h3_indices(df: pd.DataFrame, h3_res: int) -> pd.DataFrame:
    df["origin"] = df.apply(
        lambda r: h3.latlng_to_cell(r["vendor_lat"], r["vendor_lng"], h3_res), axis=1
    )
    df["destination"] = df.apply(
        lambda r: h3.latlng_to_cell(r["customer_lat"], r["customer_lng"], h3_res), axis=1
    )
    return df


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def _entropy(series: pd.Series) -> float:
    counts = series.value_counts()
    p = counts / len(series)
    return -np.sum(p * np.log(p + 1e-10))


def compute_vendor_entropy(df: pd.DataFrame, groupby_col: str) -> pd.DataFrame:
    return (
        df.groupby(groupby_col)["vendor_id"]
        .apply(_entropy)
        .reset_index()
        .rename(columns={"vendor_id": f"vendor_entropy_{groupby_col}"})
    )


def compute_aggregated_counts(df: pd.DataFrame) -> dict:
    counts = {
        "num_vendors_origin": df.groupby("origin")["vendor_id"].nunique().reset_index().rename(columns={"vendor_id": "num_vendors_origin"}),
        "num_vendors_dest": df.groupby("destination")["vendor_id"].nunique().reset_index().rename(columns={"vendor_id": "num_vendors_dest"}),
        "num_customers_origin": df.groupby("origin")["account_id"].nunique().reset_index().rename(columns={"account_id": "num_customers_origin"}),
        "num_customers_dest": df.groupby("destination")["account_id"].nunique().reset_index().rename(columns={"account_id": "num_customers_dest"}),
        "vendor_entropy_origin": compute_vendor_entropy(df, "origin"),
        "vendor_entropy_dest": compute_vendor_entropy(df, "destination"),
    }
    return counts


def aggregate_flows(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["origin", "destination"])
        .agg(
            flow=("origin", "size"),
            euclidean_km=("euclidean_km", "mean"),
            travel_time=("travel_time", "mean"),
            avg_delivery_fee=("delivery_fee_amount_lc", "mean"),
            incentive=("delivery_incentive_lc", "mean"),
            fee_per_minute=("fee_per_minute", "mean"),
            fee_per_km=("fee_per_km", "mean"),
            pct_incentivized=("delivery_incentive_lc", lambda x: (x.fillna(0) > 0).sum() / len(x) * 100),
        )
        .reset_index()
    )


def merge_all_counts(flow_df: pd.DataFrame, counts: dict) -> pd.DataFrame:
    result = flow_df.copy()
    for df in counts.values():
        merge_col = "origin" if "origin" in df.columns else "destination"
        result = result.merge(df, on=merge_col, how="left")
    return result


# ---------------------------------------------------------------------------
# Income + population assignment
# ---------------------------------------------------------------------------

def assign_income_categories(flow_df: pd.DataFrame, gdf_income: gpd.GeoDataFrame) -> pd.DataFrame:
    for side, col in [("origin", "income_cat_origin"), ("destination", "income_cat_dest")]:
        centroids = [
            Point(h3.cell_to_latlng(cell)[1], h3.cell_to_latlng(cell)[0])
            for cell in flow_df[side]
        ]
        pts_gdf = gpd.GeoDataFrame(flow_df[[side]], geometry=centroids, crs="EPSG:4326")
        joined = gpd.sjoin(pts_gdf, gdf_income[["income_category", "geometry"]], how="left", predicate="within")
        flow_df[col] = joined["income_category"].values

    return flow_df


def assign_population(flow_df: pd.DataFrame, city_df: pd.DataFrame) -> pd.DataFrame:
    pop = (
        city_df.groupby("destination")["account_id"]
        .nunique()
        .reset_index()
        .rename(columns={"account_id": "population_dest"})
    )
    flow_df = flow_df.merge(pop, on="destination", how="left")
    flow_df["population_dest"] = flow_df["population_dest"].fillna(0).astype(int)
    return flow_df


# ---------------------------------------------------------------------------
# H3 centroid distances
# ---------------------------------------------------------------------------

def compute_h3_centroid_distance(flow_df: pd.DataFrame) -> pd.DataFrame:
    print("Computing H3 centroid euclidean distances...")
    distances = []
    for _, row in flow_df.iterrows():
        origin_lat, origin_lon = h3.cell_to_latlng(row["origin"])
        dest_lat, dest_lon = h3.cell_to_latlng(row["destination"])
        distances.append(haversine(origin_lon, origin_lat, dest_lon, dest_lat))
    flow_df["h3_euclidean_km"] = distances
    return flow_df


# ---------------------------------------------------------------------------
# OSRM driving distances
# ---------------------------------------------------------------------------

def _load_cache(cache_file: str) -> dict:
    if os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            cache = pickle.load(f)
        print(f"Loaded {len(cache)} cached distances from {cache_file}")
        return cache
    return {}


def _save_cache(cache: dict, cache_file: str) -> None:
    with open(cache_file, "wb") as f:
        pickle.dump(cache, f)


def _get_osrm_distance(origin_coords, dest_coords, osrm_url: str = OSRM_URL):
    try:
        url = f"{osrm_url}/route/v1/driving/{origin_coords[1]},{origin_coords[0]};{dest_coords[1]},{dest_coords[0]}"
        r = requests.get(url, params={"overview": "false"}, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if data["code"] == "Ok" and data["routes"]:
                return data["routes"][0]["distance"] / 1000.0
    except Exception:
        pass
    return None


def compute_driving_distances_osrm(
    flow_df: pd.DataFrame,
    batch_size: int = 100,
    delay: float = 1.0,
    osrm_url: str = OSRM_URL,
    cache_file: str = DISTANCE_CACHE_FILE,
) -> pd.DataFrame:
    print("Computing driving distances via OSRM...")
    cache = _load_cache(cache_file)
    unique_pairs = flow_df[["origin", "destination"]].drop_duplicates()
    to_query = [row for _, row in unique_pairs.iterrows() if (row["origin"], row["destination"]) not in cache]

    for i, row in enumerate(tqdm(to_query, desc="OSRM queries")):
        origin_coords = h3.cell_to_latlng(row["origin"])
        dest_coords = h3.cell_to_latlng(row["destination"])
        cache[(row["origin"], row["destination"])] = _get_osrm_distance(origin_coords, dest_coords, osrm_url)
        if (i + 1) % batch_size == 0:
            sleep(delay)
        if (i + 1) % (batch_size * 5) == 0:
            _save_cache(cache, cache_file)

    _save_cache(cache, cache_file)

    flow_df["driving_km"] = flow_df.apply(
        lambda r: cache.get((r["origin"], r["destination"]), np.nan), axis=1
    )
    flow_df["detour_factor"] = flow_df["driving_km"] / flow_df["h3_euclidean_km"]
    return flow_df


# ---------------------------------------------------------------------------
# OSMnx driving distances (alternative)
# ---------------------------------------------------------------------------

def load_osm_network(bbox, network_type="drive", cache_path=None):
    if cache_path and os.path.exists(cache_path):
        G = ox.load_graphml(cache_path)
        return G

    north, south, east, west = bbox.bounds[3], bbox.bounds[1], bbox.bounds[2], bbox.bounds[0]
    G = ox.graph_from_bbox(bbox=(north, south, east, west), network_type=network_type, simplify=True)
    G = ox.add_edge_speeds(G)
    G = ox.add_edge_travel_times(G)

    if cache_path:
        ox.save_graphml(G, cache_path)
    return G


def compute_driving_distances_osmnx(flow_df: pd.DataFrame, G) -> pd.DataFrame:
    print("Computing driving distances via OSMnx...")
    unique_pairs = flow_df[["origin", "destination"]].drop_duplicates()
    lookup = {}

    for _, row in tqdm(unique_pairs.iterrows(), total=len(unique_pairs)):
        o_lat, o_lon = h3.cell_to_latlng(row["origin"])
        d_lat, d_lon = h3.cell_to_latlng(row["destination"])
        try:
            o_node = ox.distance.nearest_nodes(G, o_lon, o_lat)
            d_node = ox.distance.nearest_nodes(G, d_lon, d_lat)
            route = nx.shortest_path(G, o_node, d_node, weight="length")
            dist = sum(G[route[i]][route[i + 1]][0]["length"] for i in range(len(route) - 1)) / 1000.0
            lookup[(row["origin"], row["destination"])] = dist
        except Exception:
            lookup[(row["origin"], row["destination"])] = None

    flow_df["driving_km"] = flow_df.apply(lambda r: lookup.get((r["origin"], r["destination"]), np.nan), axis=1)
    flow_df["detour_factor"] = flow_df["driving_km"] / flow_df["h3_euclidean_km"]
    return flow_df


# ---------------------------------------------------------------------------
# Gap filling
# ---------------------------------------------------------------------------

def fill_income_from_neighbors(flow_df: pd.DataFrame, col: str) -> pd.DataFrame:
    income_col = "income_cat_origin" if col == "origin" else "income_cat_dest"

    for idx, row in flow_df[flow_df[income_col].isna()].iterrows():
        cell = row[col]
        neighbors = list(h3.grid_ring(cell, 1))
        neighbor_incomes = flow_df.loc[flow_df[col].isin(neighbors), income_col].dropna()

        if len(neighbor_incomes) == 0:
            flow_df.at[idx, income_col] = 3
        else:
            flow_df.at[idx, income_col] = float(neighbor_incomes.median())

    return flow_df


def fill_travel_time_with_nearest(flow_df: pd.DataFrame) -> pd.DataFrame:
    known = flow_df[flow_df["travel_time"].notna()][["driving_km", "travel_time"]].copy()
    missing_idx = flow_df[flow_df["travel_time"].isna()].index

    for idx in missing_idx:
        nearest = known.loc[np.abs(known["driving_km"] - flow_df.at[idx, "driving_km"]).idxmin(), "travel_time"]
        flow_df.at[idx, "travel_time"] = nearest

    return flow_df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def build_dubai_flows(
    raw_df: pd.DataFrame,
    gdf_webmerc: gpd.GeoDataFrame,
    h3_res: int = RESOLUTION,
    use_osrm: bool = True,
    osrm_url: str = OSRM_URL,
    osm_network=None,
) -> pd.DataFrame:
    print("Step 1: Computing basic metrics...")
    df = compute_basic_metrics(raw_df)

    print("Step 2: Creating bounding box...")
    bbox = create_bbox(**DUBAI_BBOX)

    print("Step 3: Filtering income polygons...")
    gdf_filtered = filter_income_polygons(gdf_webmerc, bbox)

    print("Step 4: Filtering trips to bbox...")
    city_df = filter_trips_to_bbox(df, bbox)
    print(f"  {len(city_df)} trips within bbox")

    print("Step 5: Assigning H3 indices...")
    city_df = assign_h3_indices(city_df, h3_res)

    print("Step 6: Aggregating counts and flows...")
    counts = compute_aggregated_counts(city_df)
    flow_df = aggregate_flows(city_df)
    flow_df = merge_all_counts(flow_df, counts)

    print("Step 7: Assigning income categories and population...")
    flow_df = assign_income_categories(flow_df, gdf_filtered)
    flow_df = assign_population(flow_df, city_df)

    print("Step 8: Computing distances...")
    flow_df = compute_h3_centroid_distance(flow_df)

    if use_osrm:
        flow_df = compute_driving_distances_osrm(flow_df, osrm_url=osrm_url)
    else:
        G = osm_network or load_osm_network(bbox)
        flow_df = compute_driving_distances_osmnx(flow_df, G)

    print(f"\nDone. {len(flow_df)} OD pairs, columns: {list(flow_df.columns)}")
    return flow_df
