# exprs-app-r

Universal application for expression analysis.

## Running the app

From the repository root, start the Shiny app (for example):

```r
shiny::runApp("app.R")
```

Configuration lives in `config.yaml` (app title and study list). Per-study settings and data paths are under `data/<study_id>/config.yaml`. Count matrices and metadata are expected in the locations referenced there.

Per-study count matrix config supports both:

- top-level `counts_file` (legacy and still supported), and
- optional per-comparison override `deg_lists[].counts_file`.

When `deg_lists[].counts_file` is set, DEGs and ORA use that matrix for the selected DEG list; otherwise they fall back to top-level `counts_file`.

For the Gene tab:

- if top-level `counts_file` exists and is valid, the app uses that single matrix for the study;
- otherwise the app groups DEG lists by resolved count matrix and renders one plot panel per unique matrix (panel title is DEG `label` list joined with commas).

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

### Optional raster backends for heatmaps

The heatmap code uses raster rendering for performance and selects a backend at runtime:

- If **ragg** is installed, it uses `agg_png`.
- Otherwise it falls back to base `png`.

No ImageMagick / **magick** package is required for this path.

### Gene tab (optional boxplots)

The **Gene** tab draws per-study gene boxplots with `ggplot2` / **ggpubr** / **ggnewscale**. If these are missing, the tab still loads, but plots show a short message instead of the figure.

```r
install.packages(c("ggplot2", "ggpubr", "ggnewscale"))
```

`htmltools` (used for escaping in gene info HTML) is installed as a dependency of **shiny** in typical setups.

### ORA tab (overrepresentation)

The **ORA** tab uses **clusterProfiler** (`enricher`) on a chosen pathway database file under `databases/pathways/` (all pathways in that file). For each study, enrichment runs **for every DEG list** in that study’s `config.yaml`, using **thresholds from the config** (study-level `thresholds` and optional per-list `deg_lists[].thresholds`), merged with the same defaults as the heatmap (`deg_list_threshold_defaults()`). A faceted ORA dot plot shows **comparisons (DEG list labels) on the x-axis**, **pathways on the y-axis**, **Gene ratio** as a blue-red colour scale, and **overlap count** as point size (p-values are not mapped to aesthetics; they are used internally to rank pathways).

The ORA dot plot is interactive with **ggiraph**:
- click a dot and interpret it as **Pathway (row intent)** to build a heatplot of pathway genes x all study comparisons (cell = `log2FC`);
- click a dot and interpret it as **Comparison (column intent)** to build a classical heatplot of pathways x genes for the selected comparison (cell = `log2FC`).

Sidebar filters for ORA include minimum overlap count, minimum pathway size, minimum gene ratio, and **maximum adjusted p-value (FDR)**. Default FDR cutoff is `1` (no FDR filtering), so behavior stays as before unless you choose a stricter threshold. Top-`N` pathway display still follows existing ranking by enrichment significance (`p_adj` / p-value order).

Install from Bioconductor, for example:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

**ggplot2** and **ggiraph** are required for ORA plotting/interactivity:

```r
install.packages(c("ggplot2", "ggiraph"))
```

**Custom ontology (ORA sidebar):** upload an `.xlsx` in **long** format: **column 1 = gene symbol**, **column 2 = category** (pathway). This differs from Enrichr-style `.txt` files in the dropdown (one line per pathway, genes across columns). Here each row assigns one gene to one category; the app derives **one gene set per category** (all symbols in rows with that category), then builds the same `TERM2GENE` table `enricher` expects. **Gene symbols may be lower- or mixed-case in the file; the app uppercases them** so matching to DEG tables is case-insensitive. Category text is kept as in the sheet (trimmed). Named columns `symbol` / `category` are preferred (`gene` / `pathway` allowed; otherwise the first two columns are used). The app runs the same `clusterProfiler::enricher` pipeline as for built-in pathway files (no RDS cache for this mode). Install **readxl**:

```r
install.packages("readxl")
```

Optional key **`pathways_list`** can be set in root `config.yaml` and/or per-study `data/<study_id>/config.yaml`.

- In root `config.yaml`, `pathways_list` can be either:
  - a vector of pathway filenames (e.g. `KEGG_2019_Mouse.txt`), or
  - a path to a list file (e.g. `databases/pathways_list.yaml`).
- In per-study config, `pathways_list` is a vector of pathway filenames.

If any configured list is present, the ORA pathway dropdown uses the union of those names and does not scan all files under `databases/pathways/`. If none are configured, the app falls back to scanning `databases/pathways/*.txt` as before.

To generate a ready-to-use list file from the current pathway folder, run:

```r
Rscript scripts/generate_pathways_list.R --out databases/pathways_list.yaml
```

This writes `databases/pathways_list.yaml` with a `pathways_list:` key. You can reference it directly from root config:

```yaml
pathways_list: databases/pathways_list.yaml
```

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

### Custom modules

Module: Mechanism-driven mitochondrial gene sets

Description:
Genes involved in the formation and maintenance of mitochondrial cristae structure,
including MICOS complex components and ATP synthase–mediated membrane curvature.

Sources:
- Pfanner et al., Nat Rev Mol Cell Biol (2014)
- Rampelt et al., J Cell Biol (2017)
- MitoCarta3.0 (for gene localization)