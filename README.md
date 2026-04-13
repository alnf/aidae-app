# exprs-app-r

Universal application for expression analysis.

## Running the app

From the repository root, start the Shiny app (for example):

```r
shiny::runApp("app.R")
```

Configuration lives in `config.yaml` (app title and study list). Per-study settings and data paths are under `data/<study_id>/config.yaml`. Count matrices and metadata are expected in the locations referenced there.

## R package dependencies

The dashboard is not shipped as a formal R package (there is no `DESCRIPTION`), so dependencies are listed here.

### Required for the main app

These are loaded from `app.R` and are needed for heatmaps, tables, layout, and authentication:

| Package | Role |
|--------|------|
| **shiny** | App framework |
| **bs4Dash** | UI layout |
| **shinymanager** | Login / auth |
| **InteractiveComplexHeatmap**, **ComplexHeatmap**, **circlize** | Interactive heatmaps |
| **yaml** | Reading `config.yaml` and study configs |
| **DT** | Result tables |
| **GetoptLong** | Used in the app (brush / utilities) |

Install missing packages from CRAN, for example:

```r
install.packages(c(
  "shiny", "bs4Dash", "shinymanager",
  "InteractiveComplexHeatmap", "ComplexHeatmap", "circlize",
  "yaml", "DT", "GetoptLong"
))
```

**ComplexHeatmap** and **InteractiveComplexHeatmap** may need [Bioconductor](https://bioconductor.org/) if you install them with `BiocManager::install()` rather than CRAN mirrors that carry them.

### Gene tab (optional boxplots)

The **Gene** tab draws per-study gene boxplots with `ggplot2` / **ggpubr** / **ggnewscale**. If these are missing, the tab still loads, but plots show a short message instead of the figure.

```r
install.packages(c("ggplot2", "ggpubr", "ggnewscale"))
```

`htmltools` (used for escaping in gene info HTML) is installed as a dependency of **shiny** in typical setups.

### ORA tab (overrepresentation)

The **ORA** tab uses **clusterProfiler** (`enricher`) on a chosen pathway database file under `databases/pathways/` (all pathways in that file). For each study, enrichment runs **for every DEG list** in that study’s `config.yaml`, using **thresholds from the config** (study-level `thresholds` and optional per-list `deg_lists[].thresholds`), merged with the same defaults as the heatmap (`deg_list_threshold_defaults()`). A **single ggplot** uses **facets by study** so you can compare studies side by side: **comparisons (DEG list labels) on the x-axis**, **pathways on the y-axis**, **Gene ratio** as a **blue–red** colour scale, and **overlap count** as point size (p-values are not mapped to aesthetics; they are only used internally to pick the top pathways to display). Install from Bioconductor, for example:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

**ggplot2** is required for plots (often already installed with the Gene tab).

### Gene tab: precomputed DE long file (`gdegs_file`)

Optional per-study key in `data/<study_id>/config.yaml`: **`gdegs_file`** — path **relative to that study directory** (e.g. `degs/gdegs_long.tsv`). The app loads this file for p-value brackets on Gene-tab boxplots; if the key is missing or the file is absent, plots still work without brackets.

Build or refresh the file **outside** the running app from the repository root, for example:

```r
source("scripts/gdf_utils.R")
df <- build_gene_deg_long("your_study_id")
write_gene_deg_long(df, "your_study_id")
```

`write_gene_deg_long()` writes to `data/<study_id>/<gdegs_file>` from config. Optional **`gene_tab_facet`** in the same config controls metadata column used for faceting (see `scripts/gdf_utils.R`).

### Password hashing utility

`scripts/hash_password.R` can use the **scrypt** package if present; it falls back otherwise. Only relevant if you manage `shinymanager` credentials with that script.
