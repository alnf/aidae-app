# exprs-app-r

Universal application for expression analysis.

## Running the app

From the repository root, start the Shiny app (for example):

```r
shiny::runApp("app.R")
```

Configuration lives in `config.yaml` (app title and study list). Per-study settings and data paths are under `data/<study_id>/config.yaml`. Count matrices and metadata are expected in the locations referenced there. On the **UpSet** tab, studies are coloured with the ColorBrewer **Set3** qualitative palette in the same order as studies appear in the main `config.yaml` list (row stripes and Study swatch legend).

Per-study count matrix config supports both:

- top-level `counts_file` (legacy and still supported), and
- optional per-comparison override `deg_lists[].counts_file`.

When `deg_lists[].counts_file` is set, DEGs and ORA use that matrix for the selected DEG list; otherwise they fall back to top-level `counts_file`.

For the Gene tab:

- if top-level `counts_file` exists and is valid, the app uses that single matrix for the study;
- otherwise the app groups DEG lists by resolved count matrix and renders one plot panel per unique matrix (panel title is DEG `label` list joined with commas).

## Posit Connect deploy

Per-app definitions live under [`deploy/apps/`](deploy/apps/) as YAML. **Real files** such as `all.yaml`, `heart.yaml`, and `placenta.yaml` are **gitignored** because they include Posit Connect identifiers (`deploy.account`, `deploy.title`, `deploy.app_id`). In git you only have **templates**: `all.example.yaml`, `heart.example.yaml`, `placenta.example.yaml` (placeholders, safe to commit). Locally, copy a template to the real name and fill in `deploy:`:

```bash
cp deploy/apps/all.example.yaml deploy/apps/all.yaml
# edit deploy/apps/all.yaml — set account, title, and app_id from the Connect UI (after first publish)
```

Each file combines **runtime** keys (same idea as root `config.yaml`: `title`, `studies`, `pathways_list`, …), optional `extends: config.yaml`, a **`deploy:`** block (`account`, `title`, and `app_id` for `rsconnect deploy`; `app_id` is optional when using **`--new`** below), and a **`bundle:`** block (omit `databases/`, `orig/`, swap redundant count TSV for sibling `.rds`/`.eds` when present, ORA validation mode).

- **Build manifest and deploy:** `./deploy/deploy.sh heart` or `./deploy/deploy.sh deploy/apps/heart.yaml`. Add `--dry-run` to print bundle paths and approximate size without writing `manifest.json`. Add **`--verbose`** or **`-v`** for per-study bundle messages and `rsconnect::writeManifest(verbose = TRUE)` (dependency capture is clearer but still slow). Add **`--debug`** for a full bundle path list and `set -x` during the upload step. For a **first-time** publish to Connect (no content GUID yet), run with **`--new`** (for example `./deploy/deploy.sh newstudy --new --debug`); rsconnect creates new content and prints its URL—copy the **`app_id`** into the YAML for later updates without `--new`.
- **“All studies” wrapper:** [`rsconnect.sh`](rsconnect.sh) runs `./deploy/deploy.sh all` (expects a local **`deploy/apps/all.yaml`**, typically created from [`deploy/apps/all.example.yaml`](deploy/apps/all.example.yaml)).
- **Runtime config on the server:** set environment variable **`EXPRS_MAIN_CONFIG`** to the path of the same YAML file **inside the deployed bundle** (for example `deploy/apps/heart.yaml`). The app falls back to `config.yaml` when unset.
- **Thin bundles (`bundle.omit_databases: true`):** do not upload `databases/pathways/`. Use an **inline `pathways_list`** in the deploy YAML (or a repo-relative list file such as `databases/pathways_list.yaml`, which is still bundled) and **precompute** per-study ORA (`ora_file`, default `ora/enrichment.rds`) for every pathway in that list. The manifest step validates RDS coverage (`bundle.ora_validate`: `strict`, `warn`, or `skip`). **Custom ontology** (.xlsx) still runs live `enricher` and does not require pathway files on disk. If the deploy YAML uses **`extends:`**, the base file (e.g. `config.yaml`) is included in the bundle so the app can merge the same keys at runtime.
- **Secrets:** store Connect API keys outside the repo (environment variables or `rsconnect` account configuration). If an API key was ever committed, rotate it on the server.

You need the **rsconnect** R package for `deploy/write_manifest.R`. The first run may take several minutes while dependencies are captured. From the repo root, `EXPRS_DEPLOY_REPO` defaults to the current working directory.

## R package dependencies

The dashboard is not shipped as a formal R package (there is no `DESCRIPTION`), so dependencies are listed here.

### Required for the main app

These are loaded from `app.R` and are needed for heatmaps, tables, layout, and authentication:

| Package | Role |
|--------|------|
| **shiny** | App framework |
| **bs4Dash** | UI layout |
| **shinymanager** | Login / auth |
| **InteractiveComplexHeatmap**, **ComplexHeatmap**, **circlize** | Interactive heatmaps and **UpSet** tab (`UpSet()`, `make_comb_mat()`, row annotations) |
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

**Parallel ORA (optional):** work is split with `parallel::mclapply` (forking on Unix/macOS; Windows stays effectively sequential). Set **`EXPRS_ORA_WORKERS`** to an integer ≥ `2`, or use **`scripts/precompute_ora.R --cores N`**. The effective worker count is capped by `parallel::detectCores()` (logical CPUs) and by the **number of parallel tasks** at the active layer (set workers lower if you want to leave cores free).

**Layers (only one layer runs in parallel at a time):**

- **Shiny ORA tab / custom ontology:** one pathway database (or one custom ontology) per run. If a study has several DEG lists, jobs are **parallel across studies × DEG lists** for that ontology. With **one DEG list per study**, there is only one job per study, so you still only get parallelism when **multiple studies** each contribute a task.
- **`scripts/precompute_ora.R`:** (1) If this run includes **several pathway / ontology files**, precompute uses **ontology-parallel** (`mclapply` over files); **DEG lists run sequentially** inside each ontology (avoids nested parallel). (2) If there is **exactly one** pathway file for the study, precompute can use **DEG-list parallel** when there are **≥2 DEG lists** (same idea as the Shiny ORA tab for one database). (3) On a **full multi-study batch** with **exactly one pathway file** in the run and **at most one DEG list per study**, it can run **studies in parallel**. Nested `mclapply` is deliberately avoided.
- **Incremental RDS (`--incremental` or `EXPRS_ORA_INCREMENTAL=1`):** the script updates the per-study `enrichment.rds` so a crash leaves partial work. **Single ontology:** merge-save **after each DEG list** (sequential ontologies). **Several ontologies + parallel:** for each DEG list, `enricher` runs **in parallel across all ontologies**, then the parent performs **one** read–merge–write for that comparison across all ontologies (no file locking). Re-run the same command to **resume** (`[resume-skip]`). For a full refresh, remove the RDS or use **`--no-incremental`**.
- **Heavy / light ontology buckets (`databases/pathways_list.yaml`):** optional keys `pathways_heavy` and `pathways_list` define the master list; **`--light`** runs `setdiff(pathways_list, pathways_heavy)` (or explicit `pathways_light` if set), **`--heavy`** runs `pathways_heavy` only. Logs include **`[ontology-done]`** / **`[deg-round]`** timings for profiling.

Per-step Shiny progress updates run only in single-core mode; with parallel ORA the progress bar advances once when the parallel phase finishes. **`--cores 1`** forces sequential precompute and overrides `EXPRS_ORA_WORKERS` for that run. R/BLAS may still use extra threads per process unless you set **`OMP_NUM_THREADS=1`** (and similar), so `htop` can show more busy CPUs than the ORA worker count.

Install from Bioconductor, for example:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
```

**ggplot2** and **ggiraph** are required for ORA plotting/interactivity:

```r
install.packages(c("ggplot2", "ggiraph"))
```

### UpSet tab (DEG list intersections)

The **UpSet** tab shows intersections of significant genes across **every DEG list** from **all studies** listed in the main `config.yaml`, using the same filtering as the heatmap (`load_study_data()` + `filter_heatmap_row_index()` + `deg_list_threshold_defaults()`). **Sidebar thresholds apply only to the UpSet tab** (they are not the DEGs-tab sliders). Use **Refresh thresholds for selected list** to store overrides for one list, **Apply current thresholds to all lists** to copy the sidebar values onto every configured DEG list for this session, or **Reset all lists to default thresholds** to drop every override and return to YAML defaults. **Row (set) order** above the plot can follow **study then list label** only, or **largest intersection overlap bar** each set participates in (descending; list label breaks ties). **Max lists per intersection** defaults to **2** and caps intersection degree (pairwise only when set to 2). The figure is drawn with **ComplexUpset** (inclusive intersection / “intersect” semantics, same idea as the previous ComplexHeatmap `make_comb_mat(..., mode = "intersect")` setup): degree **≥ 2** only; per-set totals appear in the set-size strip. **Click an active dot** on a DEG-list row to jump to the **DEGs** tab with that study and DEG list selected and the heatmap restricted to **genes in that intersection** (then adjust thresholds or click **Generate heatmap** as usual). A **gene multiplicity** chart (ggplot2) and **Copy genes** (intersection dropdown) use the same intersection keys as the plot.

Install **ComplexUpset** from CRAN if it is not already present; **ggiraph** is required for the interactive matrix (same as the ORA tab):

```r
install.packages(c("ComplexUpset", "ggplot2", "ggiraph"))
```

The app still uses **ComplexHeatmap** elsewhere; the UpSet tab no longer depends on it for the main intersection figure.

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