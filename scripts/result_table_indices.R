# Shared logic: row indices into a DEG table matching the DEGs tab result table.
# Depends on heatmap_utils.R (filter_heatmap_row_index, default_heatmap_thresholds).

#' Row indices for genes matching the DEGs result table rules (thresholds, brush, custom list).
#'
#' @param active_study_id Current sidebar `input$study` value.
#' @param this_study_id Study whose `res`/`mm` are being evaluated.
#' @param selected_rows Brush selection (`rv$selected_rows`); only used when
#'   `identical(active_study_id, this_study_id)`.
#' @param base_mean,svalue Numeric cutoffs; use defaults from
#'   `default_heatmap_thresholds()` when columns are absent or inputs are NULL.
#' @return Integer vector of row indices into `res`, or `integer(0)`, or `NULL` if `res`/`mm` invalid.
result_table_row_indices_for_study <- function(
    res,
    mm,
    active_study_id,
    this_study_id,
    selected_rows,
    custom_genes,
    custom_genes_study,
    lock_gene_list,
    fdr,
    log2fc,
    base_mean,
    svalue) {
  if (is.null(res) || is.null(mm)) return(NULL)
  if (!is.null(selected_rows) && length(selected_rows) > 0L &&
        identical(active_study_id, this_study_id)) {
    return(selected_rows)
  }
  sval <- if ("svalue" %in% colnames(res)) {
    as.numeric(svalue)
  } else {
    default_heatmap_thresholds()$svalue
  }
  bmean <- if ("baseMean" %in% colnames(res)) {
    as.numeric(base_mean)
  } else {
    0
  }
  use_custom <- !is.null(custom_genes) &&
    length(custom_genes) > 0 &&
    (isTRUE(lock_gene_list) ||
      (!is.null(custom_genes_study) && identical(active_study_id, custom_genes_study)))
  if (use_custom) {
    genes_vec <- unique(trimws(custom_genes))
    genes_vec <- genes_vec[nzchar(genes_vec)]
    if (length(genes_vec) == 0L) return(integer(0))
    sel <- tolower(res$symbol) %in% tolower(genes_vec)
    if (!any(sel)) return(integer(0))
    res_sub <- res[sel, , drop = FALSE]
    mm_sub <- mm[res_sub$ens_gene, , drop = FALSE]
    idx_sub <- filter_heatmap_row_index(res_sub, mm_sub, fdr, bmean, log2fc, sval)
    if (is.null(idx_sub) || length(idx_sub) == 0L) return(integer(0))
    return(which(sel)[idx_sub])
  }
  idx <- filter_heatmap_row_index(res, mm, fdr, bmean, log2fc, sval)
  if (is.null(idx) || length(idx) == 0L) integer(0) else idx
}
