readData <- function(path) {
  # pathways
  ontologies <<- list("fibrosis" = paste0("pathways/custom_fibrosis.txt"),
                      "vascular" = paste0("pathways/custom_vascular.txt"),
                      "cardio" = paste0("pathways/custom_cardio.txt"),
                      "goBP" = paste0("pathways/GO_Biological_Process_2025.txt"),
                      "goCC" = paste0("pathways/GO_Cellular_Component_2025.txt"),
                      "goMF" = paste0("pathways/GO_Molecular_Function_2025.txt"),
                      "Hallmark" = paste0("pathways/MSigDB_Hallmark_2020.txt"),
                      "Reactome" = paste0("pathways/Reactome_Pathways_2024.txt"),
                      "MPOntology" = paste0("pathways/MGI_Mammalian_Phenotype_Level_4_2024.txt"))  

  # DEGs
  pe_lfc1 <<- read.table(paste0(path, "/degs/pe_lfc1.tsv"), sep="\t", header=T, check.names = F, comment.char = "")
  no_lfc058 <<- read.table(paste0(path, "/degs/wt_lfc058.tsv"), sep="\t", header=T, check.names = F, comment.char = "")
  
  # Colormap
  mcolors <- read.table(paste0(path, "/metadata/colors.tsv"), sep="\t", header=T, check.names = F, comment.char = "")
  local_mcols <- setNames(mcolors$color, mcolors$variable)
  assign("mcols", local_mcols, envir = .GlobalEnv)
  
  # Symbol to EnsemblID mapping
  t2g <<- read.table(paste0(path, "/t2g.tsv"), sep="\t", header = T)
  
  empty <- which(t2g$symbol=="")
  t2g$symbol[empty] <- t2g$msymbol[empty]

  t2gh <<- read.table(paste0(path, "/t2gh.tsv"), sep="\t", header = T)  
  agt <<- t2g[which(t2g$symbol == "Agt"),]

  # cem object
  load(paste0(path, "/degs/cem_PE.RData"))
  cem_PE <<- cem
  graphPE <<- prepareGraph(cem_PE, t2g)

  load(paste0(path, "/degs/cem_norm.RData"))
  cem_norm <<- cem_norm
  graphNorm <<- prepareGraph(cem_norm, t2g)
  
  # Metadata
  mdata <<- read.table(paste0(path, "/metadata/metadata_h.tsv"), sep="\t", header = T)
  comps <<- read.table(paste0(path, "/metadata/comps.tsv"), sep="\t", header = T)
  comps_order_no <<- c("RDP", "RAP", "ltEP")
  #mdata <- mdata[which(mdata$Region == "LV"), ]
  #mdata <<- mdata[-which(mdata$Pheno == "np"),] 

  # Expression
  local_exprs <- read.table(paste0(path, "/exprs/lcounts_heart.tsv"), sep="\t", header = T, check.names = F)
  assign("exprs", local_exprs, envir = .GlobalEnv)

  # DEGs
  types = c("_short.tsv", ".tsv")
  type = types[2]

  regions <<- c("LV", "RV", "Sept", "Ap")  

  PEpreg_WTpreg <<- list()
  PEpost_WTpost <<- list()
  PEpreg_PEpost <<- list()
  WTpreg_WTpost <<- list()
  WTpost_np <<- list()
  WTpreg_np <<- list()
    
  for (region in regions) {
    print(region)
      
    fname = paste(region, "PEd21", region, "SDd21", sep="_")
    PEpreg_WTpreg[[region]]  <<- read.table(paste0(path, "/degs/degs_", fname, type), sep="\t", header = T)
    
    fname = paste(region, "PEpp", region, "SDpp", sep="_")
    PEpost_WTpost[[region]] <<- read.table(paste0(path, "/degs/degs_", fname, type), sep="\t", header = T)
    
    fname = paste(region, "PEd21", region, "PEpp", sep="_")
    PEpreg_PEpost[[region]] <<- read.table(paste0(path, "/degs/degs_", fname, type), sep="\t", header = T)
    PEpreg_PEpost[[region]]$log2FoldChange <<- -(PEpreg_PEpost[[region]]$log2FoldChange)
    
    fname = paste(region, "SDd21", region, "SDpp", sep="_")
    WTpreg_WTpost[[region]] <<- read.table(paste0(path, "/degs/degs_", fname, type), sep="\t", header = T)
    WTpreg_WTpost[[region]]$log2FoldChange <<- -(WTpreg_WTpost[[region]]$log2FoldChange)
    
    fname = paste(region, "SDpp", region, "np", sep="_")
    WTpost_np[[region]] <<- read.table(paste0(path, "/degs/degs_", fname, type), sep="\t", header = T)

    fname = paste(region, "SDd21", region, "np", sep="_")
    WTpreg_np[[region]] <<- read.table(paste0(path, "/degs/degs_", fname, type), sep="\t", header = T)
    
    PEpreg_WTpreg[[region]] <<- PEpreg_WTpreg[[region]][-which(PEpreg_WTpreg[[region]]$ens_gene==agt$ens_gene),]
    PEpost_WTpost[[region]] <<- PEpost_WTpost[[region]][-which(PEpost_WTpost[[region]]$ens_gene==agt$ens_gene),]  
  }

  degs_list <- list(PEpreg_WTpreg=PEpreg_WTpreg, PEpost_WTpost=PEpost_WTpost,
                    PEpreg_PEpost=PEpreg_PEpost, WTpreg_WTpost=WTpreg_WTpost,
                    WTpost_np=WTpost_np, WTpreg_np=WTpreg_np)
  gdf <<- data.frame(ens_gene = character(), symbol = character(), padj = double(), log2FC = double(),
                    group1 = character(), group2 = character(), region = character(), comp = character())
  
  for (i in 1:length(degs_list)) {
    for (region in regions) {
      group1 <- comps %>% dplyr::filter(joint == names(degs_list[i])) %>% dplyr::pull(group1)
      group2 <- comps %>% dplyr::filter(joint == names(degs_list[i])) %>% dplyr::pull(group2)
      temp <- degs_list[[i]][[region]] %>% dplyr::select(ens_gene, symbol, padj, log2FoldChange)
      temp <- temp %>% dplyr::rename(log2FC = log2FoldChange)
      temp$group1 <- group1
      temp$group2 <- group2
      temp$region <- region
      temp$comp <- comps %>% dplyr::filter(joint == names(degs_list[i])) %>% dplyr::pull(name)
      gdf <<- rbind(gdf, temp)      
    }
  }  
  
}

plotGene <- function(gene, exprs, groups, degs_list, mdata, mcols, gdf, id_type = "symbol") {
  # TODO: implement ens_gene entry mode
  if (id_type=="symbol") {
    title = gene
    ens_gene = unique(t2g$ens_gene[which(t2g$symbol == gene)])
    if (length(ens_gene)>1) {
      ens_gene = ens_gene[1]
      title = paste0(title, ", multiple genes found, showing first entry only")
    }
  } else {
    ens_gene = gene
    title = unique(t2g$symbol[which(t2g$ens_gene == gene)])
    if (length(title)==1) {
      if (title=="") {
        title = unique(t2g$msymbol[which(t2g$ens_gene== gene)])
      }
    } else {
      title = title[1]
    }
  }
  
  md <- mdata[mdata$PhenoNames %in% groups, c("SampleNumber", "Region", "PhenoNames")]
  
  gexprs <- exprs[ens_gene,]
  gexprs <- t(gexprs)
  df <- merge(md, gexprs, by.x = "SampleNumber", by.y = "row.names")
  colnames(df)[4] <- "Counts"

  egene <- ens_gene
  pval_df <- gdf %>% dplyr::filter(ens_gene == egene) %>% dplyr::filter(group1 %in% groups) %>% dplyr::filter(group2 %in% groups)
  pval_df$color <- mcols[pval_df$comp]

  pval_df$y.position = 0
  maxy <- max(df$Counts)
  for (c in levels(factor(pval_df$comp))) {
    pval_df[which(pval_df$comp==c),]$y.position <- maxy+0.5
    maxy = maxy + 0.8
  }
  
  pval_df <- pval_df %>% dplyr::select(group1, group2, padj, y.position, color, region)
  pval_df$color <- as.character(pval_df$color)
  pval_df <- pval_df %>% dplyr::rename(Region = region)
  pval_df <- pval_df %>% dplyr::filter(padj <= 0.05)
  pval_df$padj <- formatC(pval_df$padj, format = "e", digits = 2)

  ggboxplot(data=df, x="PhenoNames", y="Counts", color="PhenoNames", add="jitter", title=title, ggtheme=theme_gray()) +
    scale_color_manual(values = mcols) +
    new_scale_color() +
    stat_pvalue_manual(pval_df, label = "padj", tip.length = 0.03, color = "color") +
    scale_color_identity() +
    scale_y_continuous(expand = c(0.1, 0)) +
    facet_wrap(~Region, ncol=4)  

}

plotPathways <- function(ontologies, ontology, t2gh, ns, container_width, path, view){
  print(paste0(path, "/pathways/", view, "/", ontology, ".RData"))
  load(paste0(path, "/pathways/", view, "/", ontology, ".RData"))
  myora <- myora %>% dplyr::filter(p.adjust < 0.05)
  myora$geneRatio <- sapply(myora$GeneRatio, function(x) eval(parse(text=x)))
  myora <- myora[which(myora$Count>5),] 
  myora %>% dplyr::group_by(Module) %>% dplyr::slice_max(order_by = geneRatio, n = 30) -> myora
  
  # Ugly hardcoded value for selection of colors
  strip <- strip_themed(text_y = elem_list_text(color = rep("white", 5), size=rep(12, 5), face=rep("bold", 5)),
                        background_y = elem_list_rect(fill = mcols[12:17]))
  
  plot_data <- myora %>%
    dplyr::arrange(Module, geneRatio) %>%
    dplyr::mutate(rank = row_number())
  
  plot_data$y <- 1:nrow(plot_data)
  plot_data$y <- factor(plot_data$y)

  plot_data$ID_wrapped <- str_wrap(plot_data$ID, width = 70)
  plot_data$id <- paste(plot_data$Module, plot_data$ID, sep="_")

  idnames <- plot_data$ID_wrapped
  names(idnames) <- plot_data$y
  
  max_x <- max(plot_data$geneRatio, na.rm = TRUE)
  dynamic_nudge <- -(max_x*0.01)^(0.8)
  
  p <- ggplot(plot_data, aes(geneRatio, y)) + 
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(1, 7)) +
    geom_segment(aes(xend = 0, yend = y)) +
    theme_minimal() +
    xlab("Gene Ratio") +
    ylab(NULL) + 
    facet_grid2(Module ~ ., scales = "free_y", strip = strip, 
                drop = TRUE, axes = "margins", space = "free_y") +
    # Remove the default y-axis text
    theme(
      axis.text.y = element_text(lineheight = 1.2, color = NA),
      # Increase left margin to allow full label display
      plot.margin = unit(c(2, 1, 2, 1), "lines"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "lines")
    ) +
    # Add interactive text for y-axis labels
    geom_label_interactive(
      aes(
        x = 0, 
        y = y, 
        label = ID_wrapped,
        data_id = id,
        onclick = sprintf("Shiny.setInputValue(\"%s\", this.getAttribute(\"data-id\"), {priority: \"event\"})", ns("row_click"))
      ),
      nudge_x = dynamic_nudge,
      hover_css = "fill:none;stroke:red",
      hjust = 1, 
      vjust = 0.5,
      size = 3.4
    ) +
    # Prevent clipping of the text outside the plotting area
    coord_cartesian(clip = "off") +
    scale_y_discrete(labels = idnames,
                     expand = expansion(add = c(0.7, 0.7)))
  
  return(list(plot = p, n = nrow(myora)))
}

plotPathwayHeatmap <- function(ontology, pathway, gdf, groups, path, view) {
  load(paste0(path, "/pathways/", view, "/", ontology, ".RData"))
  curmod <- sapply(strsplit(pathway,"_"), `[`, 1) 
  pathway <- sapply(strsplit(pathway,"_"), `[`, 2) 
  myora <- myora %>% dplyr::filter(p.adjust < 0.05)  
  myora <- myora %>% dplyr::filter(Module == curmod)  
  myora <- myora %>% dplyr::filter(ID == pathway)

  geneID <- unlist(strsplit(myora$geneID, "/"))
  df <- gdf %>% dplyr::filter(ens_gene %in% geneID) %>% dplyr::filter(group1 %in% groups) %>% dplyr::filter(group2 %in% groups)

  # Create a list to store a heatmap for each region
  hts <- list()
  
  
  min_ltm <- min(df$log2FC)
  max_ltm <- max(df$log2FC)
  
  if (min_ltm>0) {
    col_fun = colorRamp2(c(0, max_ltm), c("white", "red"))
  } else if (max_ltm<0) {
    col_fun = colorRamp2(c(min_ltm, 0), c("blue","white"))
  } else {
    col_fun = colorRamp2(c(min_ltm, 0, max_ltm), c("blue", "white", "red"))
  }
  
  lgd = Legend(col_fun = col_fun, title = "log2FC", direction="vertical")  

  for (reg in regions) {
    df_wide <- df %>% dplyr::filter(region==reg)
    df_wide <- df_wide %>% dplyr::select(ens_gene, log2FC, comp, region) %>%
      pivot_wider(names_from = comp, values_from = log2FC)
    
    df_wide <- as.data.frame(df_wide)
    rownames(df_wide) <- df_wide$ens_gene
    df_wide <- df_wide[,3:ncol(df_wide)]
    df_wide <- round(df_wide, 2)
    df_wide <- df_wide[order(rownames(df_wide)),]

    df_adj <- df %>% dplyr::filter(region==reg)
    df_adj <- df_adj %>% dplyr::select(ens_gene, padj, comp, region) %>%
      pivot_wider(names_from = comp, values_from = padj)
    df_adj <- as.data.frame(df_adj)
    rownames(df_adj) <- df_adj$ens_gene
    df_adj <- df_adj[,3:ncol(df_adj)]
    df_adj <- df_adj[order(rownames(df_adj)),]
    
    df_adj[is.na(df_adj)] <- 999
    
    rowlabels <- c()
    rowlabels$ens_gene <- rownames(df_wide) 
    rowlabels <- merge(rowlabels, df[,c("ens_gene", "symbol")], by="ens_gene")
    rowlabels <- rowlabels[!(duplicated(rowlabels$ens_gene)),]
    
    csplit = c("condition", "condition", "time", "time")
    if (view=="norma") {
      csplit = c("time", "time", "time")
      df_wide <- df_wide[,comps_order_no]
      df_adj <- df_adj[,comps_order_no]
    }
    #print(df_wide)
    #print(df_adj)
    hts[[reg]] <- Heatmap(as.matrix(df_wide), name = reg, cluster_rows = F, cluster_columns = F, col=col_fun,
                  column_split = csplit, column_title = reg, border = T,
                  row_labels = rowlabels$symbol, show_heatmap_legend = F,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(df_adj[i, j] < 0.05)
                      grid.text(sprintf("%.1f", df_wide[i, j]), x, y, gp = gpar(fontsize = 10))
                  }
    )  

  }

  n_rows <- nrow(hts[[1]]@matrix)

  ht_combined <- hts[[1]] + hts[[2]] + hts[[3]] + hts[[4]]
  
  ht_opt$TITLE_PADDING = unit(c(5, 3), "mm")
  ht_opt$ANNOTATION_LEGEND_PADDING = unit(5, "mm")
  ht_opt
  ht <- draw(ht_combined, padding = unit(c(7, 7, 7, 7), "mm"), column_title=str_wrap(pathway, width = 70),
       column_title_gp=grid::gpar(fontsize=12, fontface="bold"),
       annotation_legend_list = list(lgd), legend_grouping = "original",
       ht_gap = unit(7, "mm"))
  
  return(list(ht = ht, n_rows = n_rows))
}

plotWGCNA <- function(mcols, exprs, mdata, groups, view) {
  region <- "LV"
  md <- mdata[,c("SampleNumber","Region", "PhenoNames")]
  md <- md[md$PhenoNames %in% groups,]
  md <- md[md$Region == region,]
  
  if (view=="PE") {
    m <- exprs[pe_lfc1$ens_gene, md$SampleNumber]
    
    preg = PEpreg_WTpreg[[region]]
    post = PEpost_WTpost[[region]]
    deltaPE = PEpreg_PEpost[[region]]
    RAP = WTpreg_WTpost[[region]]
    
    dlists <- list(preg=preg[which(abs(preg$log2FoldChange)>1),], post=post[which(abs(post$log2FoldChange)>1),],
                   deltaPE=deltaPE[which(abs(deltaPE$log2FoldChange)>1),], RAP=RAP[which(abs(RAP$log2FoldChange)>1),])
    cols <- list(preg=mcols["preg"], post=mcols["post"], deltaPE=mcols["deltaPE"], RAP=mcols["RAP"])
    acols <- list(Group = mcols)
    modcols <- mcols[12:17]
    
    rspl <- data.frame(modules = cem_PE@module$modules)
    rspl$ens_gene <- rownames(m)
    rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
    rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "M3", "M4", "M5", "NC"))
    
    glist <- pe_lfc1
    glist$ens_gene <- pe_lfc1$ens_gene[order(match(pe_lfc1$ens_gene, rspl$ens_gene))]
    
    order <- c("PEpreg", "PEpost", "WTpreg", "WTpost")
    pheno <- order
    ht <- compHeatmap(pheno, region, order, mdata, exprs, glist, cc=T, csplit=T, dlists, cols, acols, clegend=T, rsplit=rspl$modules, modcols)
  } else if (view=="norma") {
    m <- exprs[no_lfc058$ens_gene, md$SampleNumber]
    
    RAP = WTpreg_WTpost[[region]]
    ltEF = WTpost_np[[region]]
    RDP = WTpreg_np[[region]]
    
    dlists <- list(ltEF=ltEF[which(abs(ltEF$log2FoldChange)>0.58),], RDP=RDP[which(abs(RDP$log2FoldChange)>0.58),],
                   RAP=RAP[which(abs(RAP$log2FoldChange)>0.58),])
    cols <- list(ltEF=mcols["ltEP"], RDP=mcols["RDP"], RAP=mcols["RAP"])
    acols <- list(Group = mcols)
    modcols <- mcols[c(12:13,17)]
    
    rspl <- data.frame(modules = cem_norm@module$modules)
    rspl$ens_gene <- rownames(m)
    rspl$modules[which(rspl$modules=="Not.Correlated")] <- "NC"
    rspl$modules <- factor(rspl$modules, levels = c("M1", "M2", "NC"))
    
    glist <- no_lfc058
    glist$ens_gene <- no_lfc058$ens_gene[order(match(no_lfc058$ens_gene, rspl$ens_gene))]
    
    order <- c("np", "WTpreg", "WTpost")
    pheno <- order
    ht <- compHeatmap(pheno, region, order, mdata, exprs, glist, cc=T, csplit=T, dlists, cols, acols, clegend=T, rsplit=rspl$modules, modcols)
  }
  return(ht)
}

compHeatmap <- function(pheno, region, order=NA, metadata, exprs, glist, cc = F, csplit=T, dlists, cols, annoCol, clegend=F, rsplit=NULL, modcols){
  meta <- metadata[which(metadata$Region == region), ]
  meta <- meta[which(meta$PhenoNames %in% pheno), ]
  print(meta$PhenoNames)
  m <- exprs[glist$ens_gene, which(colnames(exprs) %in% meta$SampleNumber)]
  
  annoRow <- list()
  for (i in 1:length(dlists)) {
    genes <- dlists[[i]]$ens_gene
    anno <- rep(NA, nrow(m))
    names(anno) <- rownames(m)
    anno[which(names(anno) %in% genes)] <- cols[[names(dlists)[i]]]
    anno <- list(anno) 
    names(anno) <- names(dlists)[i]
    annoRow <- append(annoRow, anno)
  }
  
  if (!is.null(order)) {
    meta$PhenoNames <- factor(meta$PhenoNames, levels = order)
    meta <- meta[order(meta$PhenoNames),]
    m <- m[,match(meta$SampleNumber, colnames(m))]
  }
  
  if (csplit) {
    cspl <- meta$PhenoNames
  } else {
    cspl = NULL
  }
  
  anno_df = data.frame(matrix(NA, nrow = nrow(m), ncol = length(names(dlists))))
  for (i in 1:ncol(anno_df)) {
    anno_df[,i] <- rownames(m)
  }
  colnames(anno_df) <- names(dlists)
  rha = rowAnnotation(df=anno_df, col=annoRow, show_legend = F)
  cha = HeatmapAnnotation(Group = meta[,c("PhenoNames")], col=annoCol, show_legend = clegend, show_annotation_name = F,
                          annotation_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 13, fontface = "bold"),
                                                         grid_height = unit(7, "mm"), grid_width = unit(7, "mm"),
                                                         labels = order))
  
  labels = names(modcols)
  labels[which(labels=="NC")] <- ""
  
  lha = rowAnnotation(Module = anno_block(gp = gpar(fill = modcols), labels = labels))
  if (cc) {
    hc <- hclust(dist(t(m)), method="average")
    dd <- as.dendrogram(hc)
    dd <- dendextend::rotate(dd, as.character(meta$SampleNumber))
    #dd <- reorder(dd, meta$SampleNumber)
    cspl = length(pheno)
    ctitle = NULL
  } else {
    dd = cc
    ctitle = order
  }
  
  ht = Heatmap(t(scale(t(m))), show_row_names = F, show_row_dend = F, show_column_names = T, cluster_columns = dd,
               top_annotation = cha, right_annotation = rha, left_annotation = lha,
               name = "expr",
               column_split = cspl, row_split = rsplit, cluster_row_slices = F, cluster_column_slices = T,
               row_names_gp = gpar(fontsize = 14),
               column_names_gp = gpar(fontsize = 14),
               column_dend_height=unit(20, "mm"),
               heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 13),
                                           legend_direction = "horizontal", title_position = "topcenter"),
               row_title=NULL, column_title=ctitle
  )
  return(ht)
}

prepareGraph <- function(cem, t2g) {
  df <- data.frame(ens_gene = cem@module$genes, modules = cem@module$modules)
  df <- merge(df, t2g, by="ens_gene")
  
  df$symbol[which(df$symbol=="")]=df$msymbol[which(df$symbol=="")]
  rat_ids <- df[-which(duplicated(df$ens_gene)),]
  str_ids <- rba_string_map_ids(ids = rat_ids$symbol, species = 10116)
  rat_ids <- merge(rat_ids, str_ids, by.x="symbol", by.y="queryItem")
  
  edges <- rba_string_interactions_network(str_ids$stringId, species = 10116, required_score=700)
  g <- igraph::graph_from_data_frame(edges[3:ncol(edges)], directed=FALSE)
  g_tidy <- as_tbl_graph(g)
  
  rat_ids_unique <- rat_ids %>% 
    dplyr::distinct(preferredName, .keep_all = TRUE)
  
  g_tidy_merged <- g_tidy %>%
    dplyr::mutate(name = as.character(name)) %>%
    dplyr::left_join(rat_ids_unique, by = c("name" = "preferredName"))
  
  fixed_layout <- create_layout(g_tidy_merged, layout = "fr")

  return(list(rat_ids, fixed_layout))  
}

plotNetwork <- function(ontologies, ontology, pathway, path, view){
  if (view=="PE") {
    rat_ids <- graphPE[[1]]
    fixed_layout <- graphPE[[2]]
    nodcols <- mcols[12:17]
    names(nodcols) <- c("M1", "M2", "M3", "M4", "M5", "NA")
  } else if (view=="norma") {
    rat_ids <- graphNorm[[1]]
    fixed_layout <- graphNorm[[2]]
    nodcols <- mcols[c(12:13,17)]
    names(nodcols) <- c("M1", "M2", "NA")
  }  

  if (!is.null(ontology)) {
    load(paste0(path, "/pathways/", view, "/", ontology, ".RData"))
    curmod <- sapply(strsplit(pathway,"_"), `[`, 1) 
    pathway <- sapply(strsplit(pathway,"_"), `[`, 2) 
    myora <- myora %>% dplyr::filter(p.adjust < 0.05)  
    myora <- myora %>% dplyr::filter(Module == curmod)  
    myora <- myora %>% dplyr::filter(ID == pathway)
     
    geneID <- unlist(strsplit(myora$geneID, "/"))
    fun_ids <- geneID
    ll <- list(ens_gene = fun_ids)
    ids <- merge(ll, rat_ids, by="ens_gene")
  }

  p1 <- ggraph(fixed_layout) +
    geom_edge_link(alpha = 0.5, edge_colour="black") +
    geom_node_point(aes(color = modules), size=3) +
    scale_color_manual(values = nodcols) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_graph() +
    theme(legend.position="bottom", plot.title = element_text(size=11)) +
    ggtitle("Modules")
  
  p2 <- p1
  if (!is.null(pathway)) {
    p2 <- ggraph(fixed_layout) +
      geom_edge_link(alpha = 0.5, edge_colour="black") +
      geom_node_point(aes(color = name %in% ids$symbol), size = 3) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
      theme_graph() +
      theme(legend.position="none", plot.title = element_text(size=11)) +
      ggtitle(pathway)
  }
  
  combined_plot <- p1 / p2
  combined_plot
}