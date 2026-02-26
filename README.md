# Distance Decay in Online Deliveries and the Sustainability Impact of Pricing Incentives

**Javad Eshtiyagh¹·³ · Anupriya¹ · Daniel Hörcher¹·² · Paolo Santi³ · Daniel J. Graham¹***

¹ Transport Strategy Centre, Imperial College London, London, UK  
² Institute for Advanced Studies, Corvinus University of Budapest, Budapest, Hungary  
³ Senseable City Lab, Massachusetts Institute of Technology, Cambridge MA, USA  
\* Corresponding author: d.j.graham@imperial.ac.uk

**Preprint:** [Link to be added]

---

## Abstract

The rapid expansion of last-mile delivery has reshaped urban logistics, yet the spatial patterns governing how delivery flows vary with distance and travel time remain largely unexplored, hindering efficient network design and evidence-based policy interventions. Using a novel large-scale dataset comprising tens of millions of records from multiple delivery platforms across different continents, we provide the first comprehensive empirical estimation of flow-distance relationships in last-mile delivery. Employing semi-parametric spline models and Gaussian process regression, we demonstrate that express delivery orders, including food and grocery, exhibit a universal exponential decay with distance and travel time across all cities. In contrast, centralized scheduled parcel deliveries fulfilled by a small number of peripheral warehouses display near-uniform distance distributions, largely independent of distance. Moreover, we show that pricing mechanisms can meaningfully alter these spatial patterns. In response to delivery fee incentives, customers increase ordering distances by an average of 14%, contradicting the common assumption of treating order travel distance as exogenous. These findings indicate that spatial and environmental externalities are behaviorally responsive to pricing strategies, underscoring the need for interventions that preserve price signals reflecting the true social cost of distance.

---

## Repository Structure

```
📦 repository root
├── README.md
├── code/
│   ├── 0_data_preparation/
│   │   └── flow_table_builder.py       # Builds flow tables from raw data
│   ├── 1_curve_estimation/
│   │   ├── platform1_curve_estimation.R
│   │   ├── platform1_exponential.R
│   │   ├── platform1_powerlaw.R
│   │   ├── platform1_gaussian_process.R
│   │   ├── platform2_curve_estimation.R
│   │   ├── platform3_curve_estimation.R
│   │   └── platform4_curve_estimation.R
│   ├── 2_heterogeneity_and_incentives/
│   │   ├── heterogeneity_analysis.R    # Subgroup and cross-city heterogeneity
│   │   └── platform1_incentives.py     # Pricing incentive analysis
│   ├── 3_robustness/
│   │   ├── alternative_cost_measures.R
│   │   └── alternative_spline_settings.R
│   └── 4_visualization/
│       ├── multicity_grid_curves.R
│       ├── multicity_grid_exponential.R
│       ├── stacked_bar_chart-distance.R
│       ├── stacked_bar_chart-income.R
│       └── platform1_visualizations.py
└── figures/                            # All output figures
```

---

## Data Availability

Socioeconomic data and data from **Platform #3** and **Platform #4** are publicly available:

- **Platform #3** data: Eshtiyagh et al. (2023), available at [https://arxiv.org/abs/2306.10675](https://arxiv.org/abs/2306.10675)
- **Platform #4** data: Lei et al. (2022), available at [https://pubsonline.informs.org/doi/10.1287/trsc.2022.1173](https://pubsonline.informs.org/doi/10.1287/trsc.2022.1173)

To protect privacy and commercial confidentiality, the proprietary datasets from **Platform #1** and **Platform #2** cannot be shared in full. However, de-identified and aggregated samples are available from the corresponding author upon reasonable request (d.j.graham@imperial.ac.uk).

---

## Reproducing the Analysis

Run the scripts in the following order:

1. **`0_data_preparation/`** — Build flow tables from raw input data
2. **`1_curve_estimation/`** — Estimate distance-decay curves per platform; fit exponential, power-law, and Gaussian process models
3. **`2_heterogeneity_and_incentives/`** — Run heterogeneity analysis and pricing incentive regressions
4. **`3_robustness/`** — Reproduce robustness checks and sensitivity analyses
5. **`4_visualization/`** — Generate all figures (saved to `figures/`)

---

## Dependencies

**R packages**

```r
install.packages(c("lme4", "splines", "ggplot2", "dplyr", "tidyr"))
```

**Python packages**

```bash
pip install contextily geopandas geopy h3 linearmodels matplotlib networkx \
            numpy osmnx pandas pyproj requests scipy seaborn shapely \
            statsmodels tqdm
```

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Citation

If you use this code or data in your work, please cite:

```bibtex
@article{eshtiyagh2025distancedecay,
  title   = {Distance Decay in Online Deliveries and the Sustainability Impact of Pricing Incentives},
  author  = {Eshtiyagh, Javad and Anupriya and H{\"o}rcher, Daniel and Santi, Paolo and Graham, Daniel J.},
  journal = {[Under review]},
  year    = {2025},
  url     = {[Preprint URL to be added]}
}
```
