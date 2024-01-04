
plot_detection <- function(df_detection, subgroup, min_unique, min_spat, max_sd, color, dot_alpha, num_prot_labeled=3) {
  labeled_tmp = df_detection[order(df_detection$mean_log2_int, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$uniprot[1:num_prot_labeled]
  labeled = labeled_tmp
  
  title = paste("Proteins of ", subgroup, " in function of % of detection with\nSPAT>=", min_spat, ", #unique>=", min_unique, ", sd<=", max_sd, sep="")
  gg = ggplot(df_detection, aes_string(x="detection_p", y="mean_log2_int", color="SPAT_score"))
  gg = gg + geom_point(alpha=dot_alpha)
  gg = gg + geom_label_repel(data=filter(df_detection, uniprot %in% labeled), aes(label = Gene.names), min.segment.length=0)
  gg = gg + scale_colour_gradientn(colours=c("black", color), limits=c(0, 10))
  gg = gg + labs(y = "Mean log2(intensity + 1)", x = "Percentage of detection")
  gg = gg + theme_bw()
  gg = gg + ggtitle(title)
  gg = gg + theme(
    text = element_text(size=16),
    plot.title = element_text(hjust = 0.5)
  )
  return(gg)
}


plot_single <- function (dep, tab_uni2ens_uniq_uniprot, proteins, type = c("contrast", "centered"), plot = TRUE, df_norm=NA)
{
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(proteins), is.character(type), is.logical(plot),
                          length(plot) == 1)
  type <- match.arg(type)
  row_data <- rowData(dep, use.names = FALSE)
  row_data = merge(row_data, tab_uni2ens_uniq_uniprot, by.x="name", by.y="uniprot", all.x=TRUE)
  if (any(!c("label", "condition", "replicate") %in% colnames(colData(dep)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(dep)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)), "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  if (!"name" %in% colnames(row_data)) {
    stop("'name' column not present in '", deparse(substitute(dep)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (all(!proteins %in% row_data$name)) {
    if (length(proteins) == 1) {
      rows <- grep(substr(proteins, 1, nchar(proteins) -
                            1), row_data$name)
      possibilities <- row_data$name[rows]
    }
    else {
      rows <- lapply(proteins, function(x) grep(substr(x,
                                                       1, nchar(x) - 1), row_data$name))
      possibilities <- row_data$name[unlist(rows)]
    }
    if (length(possibilities) > 0) {
      possibilities_msg <- paste0("Do you mean: '", paste0(possibilities,
                                                           collapse = "', '"), "'")
    }
    else {
      possibilities_msg <- NULL
    }
    stop("please run `plot_single()` with a valid protein names in the 'proteins' argument\n",
         possibilities_msg, call. = FALSE)
  }
  if (any(!proteins %in% row_data$name)) {
    proteins <- proteins[proteins %in% row_data$name]
    warning("Only used the following protein(s): '", paste0(proteins,
                                                            collapse = "', '"), "'")
  }
  subset <- dep[proteins]
  select_row = which(row_data$name == proteins)
  str_title =  paste(row_data$Gene.names[select_row], row_data$name[select_row], sep=" ")
  if (type == "centered") {
    means <- rowMeans(assay(subset), na.rm = TRUE)
    df_reps <- data.frame(assay(subset) - means) %>% tibble::rownames_to_column() %>%
      tidyr::gather(ID, val, -rowname) %>% dplyr::left_join(., data.frame(colData(subset)),
                                                            by = "ID")
    
    #added to display imputed data
    df_reps$imputed <- as.vector(t(is.na(assay(df_norm[proteins]))))
    df_reps$replicate <- as.factor(df_reps$replicate)
    df_tmp <- df_reps %>% dplyr::group_by(condition, rowname) %>%  summarise(num_imp=sum(imputed))
    df_tmp$condition_search = paste("^", df_tmp$condition, sep="")
    df_tmp$condition = paste(df_tmp$condition, " (", df_tmp$num_imp, ")", sep="")
    df_reps$condition <- stri_replace_all_regex(df_reps$condition, df_tmp$condition_search, df_tmp$condition, vectorize_all=FALSE)
    
    
    df <- df_reps %>% dplyr::group_by(condition, rowname) %>% summarize(mean = mean(val,
                                                                                    na.rm = TRUE), sd = sd(val, na.rm = TRUE), n = n()) %>%
      dplyr::mutate(error = qnorm(0.975) * sd/sqrt(n), CI.L = mean -
                      error, CI.R = mean + error) %>% as.data.frame()
    df$rowname <- readr::parse_factor(df$rowname, levels = proteins)
    
    p <- ggplot(df, aes(condition, mean)) + geom_hline(yintercept = 0) +
      geom_col(colour = "black", fill = "grey") + geom_point(data = df_reps, alpha=0.6,
                                                             #aes(condition, val, col = replicate), shape = 18,
                                                             aes(condition, val, shape=imputed), #shape = 18,
                                                             size = 5, position = position_dodge(width = 0.3)) +
      geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
      labs(x = "Baits", y = expression(log[2] ~ "Centered intensity" ~
                                         "(±95% CI)"), col = "Rep") +
      ggtitle(str_title) +
      theme_DEP2() + theme(legend.position="bottom")
    #facet_wrap(~rowname) +
  }
  if (type == "contrast") {
    df <- rowData(subset, use.names = FALSE) %>% data.frame() %>%
      dplyr::select(name, dplyr::ends_with("_diff"), dplyr::ends_with("_CI.L"),
                    ends_with("_CI.R")) %>% tidyr::gather(var, val, -name) %>%
      dplyr::mutate(contrast = gsub("_diff|_CI.L|_CI.R", "", var),
                    var = gsub(".*_", "", var)) %>% tidyr::spread(var, val)
    df$name <- readr::parse_factor(df$name, levels = proteins)
    suffix <- get_suffix(df$contrast)
    if (length(suffix)) {
      df$contrast <- delete_suffix(df$contrast)
    }
    p <- ggplot(df, aes(contrast, diff)) + geom_hline(yintercept = 0) +
      geom_col(colour = "black", fill = "grey") + geom_errorbar(aes(ymin = CI.L,
                                                                    ymax = CI.R), width = 0.3) + labs(x = suffix, y = expression(log[2] ~
                                                                                                                                   "Fold change" ~ "(±95% CI)")) +
      ggtitle(str_title) +
      theme_DEP2()
    #facet_wrap(~name) +
  }
  if (plot) {
    return(p)
  }
  else {
    if (type == "centered") {
      df_reps <- df_reps %>% dplyr::select(condition, bclq_id, val, imputed)
      colnames(df_reps) <- c("condition", "sample", "log2_Centred_intensity", "imputed")
      df_reps$condition = factor(df_reps$condition, levels=df$condition)
      
      return(df_reps)
    }
    if (type == "contrast") {
      df <- df %>% dplyr::select(name, contrast, diff, CI.L, CI.R) %>%
        dplyr::mutate(contrast = paste0(contrast, suffix))
      colnames(df) <- c("protein", "contrast", "log2_fold_change",
                        "CI.L", "CI.R")
    }
    return(df)
  }
}



plot_prot_intensity <- function(df_int) {
  gg = ggplot(df_int, aes(x=condition, y=.data[['log2(Normalized Imputed Intensity + 1)']])) 
  gg = gg + geom_point(aes(shape=imputed), alpha=0.5, position=position_jitter(width=0.25, height=0))
  gg = gg + stat_summary(fun = median, geom = "crossbar", width = 0.5)
  gg = gg + theme_bw()
  gg = gg + ggtitle(df_int$Gene.names[1])
  gg = gg + theme(text = element_text(size=16), 
                  plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle=45, hjust=1))
  gg = gg + scale_x_discrete()
  
  return(gg)
}



format_vol <- function(gg_vol, title, dot_color, dot_alpha, vline=0) {
  gg_vol = gg_vol + ggtitle(title)
  gg_vol = gg_vol + geom_point(alpha=dot_alpha)
  gg_vol = gg_vol + geom_vline(xintercept=vline)
  gg_vol = gg_vol + scale_colour_gradientn(colours=c("black", dot_color))
  gg_vol = gg_vol + theme_bw()
  gg_vol = gg_vol + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5))
  gg_vol = gg_vol + scale_colour_gradientn(colours=c("black", dot_color), limits=c(0, 10))
  return(gg_vol)
}

dep_volcano <- function (dep, df_anno, contrast, num_labeled, num_min_pep, min_spat, genes_highligth, dot_color, dot_alpha, plot=TRUE, adjusted=TRUE)
{
  row_data <- rowData(dep, use.names = FALSE)
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)), "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }
  if (length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) ==
      0) {
    return(ggplot())
    valid_cntrsts <- row_data %>% data.frame() %>% dplyr::select(ends_with("_diff")) %>%
      colnames(.) %>% gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"), "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg, call. = FALSE)
  }
  diff <- grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
  p_values_adj <- grep(paste(contrast, "_p.adj", sep = ""),
                       colnames(row_data))
  p_values <- grep(paste(contrast, "_p.val", sep = ""),
                   colnames(row_data))
  
  signif <- grep(paste(contrast, "_significant", sep = ""), colnames(row_data))
  unique <- grep("#Unique", colnames(row_data))
  
  df <- data.frame(log2FC = row_data[, diff],
                   pvadj = -log10(row_data[, p_values_adj]),
                   pvalues = -log10(row_data[, p_values]),
                   significant = row_data[, signif],
                   unique = row_data[, unique], #TODO add 2 unique pep for merge data
                   uniprot = row_data$name) %>%
    filter(!is.na(significant)) %>%
    arrange(significant)
  
  df = merge(df, df_anno, by="uniprot", all.x=TRUE)
  
  if (min_spat != 0) {
    df = df[which(df$SPAT_score >= min_spat), ]
  }
  
  df = df[which(df$unique >= num_min_pep), ]
  
  if (plot == FALSE) {
    return(df)
  }
  
  labeled_tmp = df[which(df$log2FC >= 0 ),]
  labeled_tmp = labeled_tmp[order(labeled_tmp$pvadj, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$Gene.names[1:num_labeled]
  labeled = labeled_tmp
  
  labeled_tmp = df[which(df$log2FC <= 0 ),]
  labeled_tmp = labeled_tmp[order(labeled_tmp$pvadj, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$Gene.names[1:num_labeled]
  labeled = c(labeled, labeled_tmp)
  
  labeled_tmp = df[order(df$log2FC, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$Gene.names[1:num_labeled]
  labeled = c(labeled, labeled_tmp)
  
  labeled_tmp = df[order(df$log2FC, decreasing=FALSE),]
  labeled_tmp = labeled_tmp$Gene.names[1:num_labeled]
  labeled = c(labeled, labeled_tmp)
  
  if (num_labeled == 0) {
    labeled = c()
  }
  
  genes_highligth = strsplit(genes_highligth, ",")[[1]]
  if(length(genes_highligth)!= 0) {
    labeled = c(labeled, genes_highligth)
  }
  
  title = paste("Volcano of ", contrast, " with #unique>=", num_min_pep, ", SPAT>=", min_spat, sep="")
  
  gg_vol =  ggplot(df, aes_string(x="log2FC", y="pvadj", color="SPAT_score"))
  if (length(labeled) != 0) {
    gg_vol = gg_vol + geom_label_repel(data=filter(df, Gene.names %in% labeled), aes(label = Gene.names), min.segment.length=0)
  }
  gg_vol = gg_vol + ylab("-log10(pv.adj)")
  
  return(format_vol(gg_vol, title, dot_color, dot_alpha))
}


plot_limma_voom <- function(limma_res, surface_score, num_gene_labeled, min_SPAT, dot_color, dot_alpha, genes_highligth=c(), table=FALSE) {
  
  if(is.null(limma_res)){
    if(table) {
      return(data.frame())
    }
    gg = ggplot()
    gg = gg + ggtitle("Data not available, call the portal's admin!")
    return(gg)
  }
  limma_res$ID = unlist(lapply(limma_res$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  limma_res$logFC = -limma_res$logFC
  limma_res$log10_adj.P.Val = -log10(limma_res$adj.P.Val)
  limma_res$log10_P.Value = -log10(limma_res$P.Value)
  
  limma_res = merge(limma_res, surface_score, by.x="ID", by.y="ENSEMBL_id", all.x=TRUE)
  if (min_SPAT >= 0){
    limma_res = limma_res[which(limma_res$SPAT_score>=min_SPAT),]
  }
  if(table) {
    return(limma_res)
  }
  
  
  
  labeled_tmp = limma_res[which(limma_res$logFC >= 0 ),]
  labeled_tmp = labeled_tmp[order(labeled_tmp$log10_adj.P.Val, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$ID[1:num_gene_labeled]
  labeled = labeled_tmp
  
  labeled_tmp = limma_res[which(limma_res$logFC <= 0 ),]
  labeled_tmp = labeled_tmp[order(labeled_tmp$log10_adj.P.Val, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$ID[1:num_gene_labeled]
  labeled = c(labeled, labeled_tmp)
  
  labeled_tmp = limma_res[order(limma_res$logFC, decreasing=TRUE),]
  labeled_tmp = labeled_tmp$ID[1:num_gene_labeled]
  labeled = c(labeled, labeled_tmp)
  
  labeled_tmp = limma_res[order(limma_res$logFC, decreasing=FALSE),]
  labeled_tmp = labeled_tmp$ID[1:num_gene_labeled]
  labeled = c(labeled, labeled_tmp)
  
  if (num_gene_labeled == 0) {
    labeled = c()
  }
  
  genes_highligth = strsplit(genes_highligth, ",")[[1]]
  if(length(genes_highligth)!= 0) {
    genes_highligth = limma_res$ID[which(limma_res$ID %in% genes_highligth | limma_res$Gene.names %in% genes_highligth)]
    labeled = c(labeled, genes_highligth)
  }
  labeled = unique(labeled)
  
  
  title = "Limma voom volcano plot with SPAT score"
  
  gg_vol = ggplot(limma_res, aes_string(x="logFC", y="log10_adj.P.Val", color="SPAT_score"))
  if (length(labeled) != 0) {
    gg_vol = gg_vol + geom_label_repel(data=filter(limma_res, ID %in% labeled), aes(label = Gene.names), min.segment.length=0)
  }
  gg_vol = gg_vol + ylab("-log10(adj.P.Val)")
  gg_vol = gg_vol + xlab("log2FC")
  
  return(format_vol(gg_vol, title, dot_color, dot_alpha))
}


plot_trans_gene <- function(trans_subgoups, gene_exp, tab_design, gene_name, graphic_type, sel_exp) {
  if (is.null(gene_exp)) return(ggplot() + ggtitle("No expression data!"))
  if (gene_name == "") return(ggplot() + ggtitle("Choose a gene!"))
  
  gene_exp = get_gene_exp(trans_subgoups, gene_exp, tab_design, gene_name, graphic_type, sel_exp)
  
  if (sel_exp == "readcounts.xls") {
    gg_exp = plot_raw_data(gene_exp, "log2_CPM", gene_name, graphic_type=graphic_type, labels=FALSE, alpha=0.5, dodge=TRUE)
  } else {
    gg_exp = plot_raw_data(gene_exp, "log2_TPM", gene_name, graphic_type=graphic_type, labels=FALSE, alpha=0.5, dodge=TRUE)
  }
  
  return(gg_exp)
}


plot_raw_data <- function(data, y_lab, title, graphic_type="scatterplot", labels=TRUE, alpha=1, dodge=FALSE) {
  gg = ggplot(data, aes_string(x="group", y=y_lab))
  if (graphic_type == "scatterplot") {
    gg = gg + stat_summary(fun = median, geom = "crossbar", width = 0.5)
    if (dodge) {gg = gg + geom_point(shape=20, alpha=alpha, position=position_jitter(width=0.25, height=0))}
    else       {gg = gg + geom_point(shape=20, alpha=alpha)}
  } else {
    gg = gg + geom_boxplot()
  }
  gg = gg + ggtitle(title)
  if (labels) gg = gg + geom_label_repel(aes(label = sample), min.segment.length=0)
  gg = gg + theme_bw()
  gg = gg + theme(
    text = element_text(size=16), 
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust = 0.5)
  )
  
  return(gg)
}


plot_sc_umap = function(type, cellTypeGrouped=T) {
  df_umap <- read.csv(file.path(path_sc, 'umap_coordinates.csv'))
  
  annotation <- scLoom[[paste("col_attrs",type, sep="/")]][]
  if (type == "Classification") {
    if (cellTypeGrouped) {
      annotation = format_classification(annotation)
    }
  }
  
  palette_choosed = c()
  switch(type,
    Classification={
      if (cellTypeGrouped) {
        palette_choosed = color_palette_grouped
      } else {
        palette_choosed = color_palette
      }
    },
    sample={
      palette_choosed = sample_palette
    },
    AML_subgroup={
      palette_choosed = aml_palette
    }
  )
  
  df_umap <- cbind(df_umap, annotation)
  colnames(df_umap)[4] = type
  
  if (type == "Classification") {
    tmpLevels = levelCellType
    if (cellTypeGrouped) {
      tmpLevels = levelCellTypeGrouped
    }
    
    df_umap$Classification = factor(df_umap$Classification, levels=tmpLevels) 
  } else {
    df_umap[,4] <- factor(df_umap[,4], levels=unique(df_umap[,4]))
  }
  
  gg = ggplot(df_umap, aes_string(x="umap_1", y="umap_2", colour=type)) 
  gg = gg + geom_point(size=0.8, aes(colour = factor(df_umap[,4]))) 
  # gg = gg + ggtitle("UMAP")
  gg = gg + theme(panel.background = element_blank(),
    plot.title = element_text(size=20,hjust = 0.5), 
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    legend.text = element_text(size = 10)
  ) 
  gg = gg + scale_colour_manual(values = palette_choosed)
  gg = gg + guides(colour = guide_legend(override.aes = list(size=10)))
  
  return(gg)
}


plot_sc_track_genes_exp = function(genes, type, cellTypeGrouped=T) {
  df_umap <- read.csv(file.path(path_sc, 'umap_coordinates.csv'))
  
  geneExps <- c()
  for (gene in genes){
    geneExps <- cbind(geneExps, as.vector(scLoom[["matrix"]][, scLoom[["row_attrs/Gene"]][] == gene]))
  }
  colnames(geneExps) <- genes
  
  palette_choosed = c()
  switch(type,
    Classification={
      if (cellTypeGrouped) {
        palette_choosed = color_palette_grouped
      } else {
        palette_choosed = color_palette
      }
    },
    sample={
      palette_choosed = sample_palette
    },
    AML_subgroup={
      palette_choosed = aml_palette
    }
  )
  
  df_show <- data.frame(
    cells = scLoom[["col_attrs/CellID"]][],
    Classification = scLoom[["col_attrs/Classification"]][],
    Annotation = scLoom[[paste("col_attrs", type, sep="/")]][],
    geneExp = geneExps
  )
  
  if (cellTypeGrouped) {
    df_show$Classification = format_classification(df_show$Classification)
    if (type == "Classification") {
      df_show$Annotation = df_show$Classification
    }
  }
  
  colnames(df_show) <- c('cellIDs','Classification', 'Annotation', genes)
  
  df_all =c()
  dfhsc = c()
  set.seed(42)
  for (celltype in unique(df_show$Classification)){
    if (dim(df_show[df_show$Classification == celltype,])[1] >= 300){
      if (celltype == 'CD34+ HSC'){
        dfhsc <- sample_n(df_show[df_show$Classification == 'CD34+ HSC',], 1000)        
      }
      else{
        df_all <- rbind(df_all, sample_n(df_show[df_show$Classification == celltype,], 300))
      }
    }
  }
  df_show <- rbind(dfhsc,df_all)

  tmpLevels = levelCellType
  if (cellTypeGrouped) {
    tmpLevels = levelCellTypeGrouped
  }
  
  df_show$Classification = factor(df_show$Classification, levels=tmpLevels) 
  if (type == "Classification") {
    df_show$Annotation = df_show$Classification
  } else {
    df_show$Annotation <- factor(df_show$Annotation, levels=unique(df_show$Annotation))
  }
  df_show <- df_show[order(df_show$Annotation),]
  df_show <- df_show[ , -which(names(df_show) %in% c("Classification"))]
  
  econdatalong <- gather(df_show, key="Gene", value="Expression",all_of(genes))
  
  gg = ggplot(econdatalong, aes(x=factor(cellIDs, levels = unique(cellIDs)), y=Expression, fill=Annotation)) 
  gg = gg + geom_bar(position = 'dodge', stat='identity', alpha=1, width = 5) 
  gg = gg + scale_fill_manual(name=type, values=palette_choosed)
  gg = gg + theme(
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y =  element_text(size=15),
    legend.position="right",
    strip.text = element_text(size=12,hjust = 0.5, face = 'bold'), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
  gg = gg + ylim(0,4) 
  gg = gg + facet_wrap(~Gene, ncol=1)
  
  return(gg)
}


plot_sc_umap_genes_exp = function(genes) {
  df_umap <- read.csv(file.path(path_sc, 'umap_coordinates.csv'))
  
  geneExps <- c()
  for (gene in genes){
    geneExps <- cbind(geneExps, as.vector(scLoom[["matrix"]][, scLoom[["row_attrs/Gene"]][] == gene]))
  }
  colnames(geneExps) <- genes
  df_umap <- cbind(df_umap,geneExps)
  econdatalong <- gather(df_umap, key="Gene", value="Expression",all_of(genes))
  
  gg = ggplot(econdatalong %>% arrange(Expression), aes(x=umap_1, y=umap_2, colour=Expression)) 
  gg = gg + geom_point(size=0.4)
  gg = gg + theme(
    panel.background = element_blank(), legend.position = "left",
    strip.text = element_text(size=12,hjust = 0.5, face = 'bold'), 
    axis.text=element_blank(),axis.ticks = element_blank(),
    axis.title=element_blank()
  )
  gg = gg + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  gg = gg + facet_wrap(~Gene)
  
  return(gg)
}



plot_sc_umap_gene_exp = function(gene) {
  df_umap <- read.csv(file.path(path_sc, 'umap_coordinates.csv'))
  
  geneExp <- as.vector(scLoom[["matrix"]][, scLoom[["row_attrs/Gene"]][] == gene])
  
  df_umap <- cbind(df_umap,geneExp)
  colnames(df_umap)[4] <- gene
  df_umap = df_umap[order(df_umap[gene]),]
  
  gg = ggplot(df_umap, aes_string(x="umap_1", y="umap_2", colour=gene)) 
  gg = gg + geom_point(size=0.8)
  gg = gg + theme(
    panel.background = element_blank(), legend.position = "left",
    strip.text = element_text(size=12,hjust = 0.5, face = 'bold'), 
    axis.text=element_blank(),axis.ticks = element_blank(),
    axis.title=element_blank()
  )
  gg = gg + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  
  return(gg)
}




format_subgroup <- function(AML_subg) {
  AML_subg[which(AML_subg == "Complex karyotype")] = "CK"
  AML_subg[which(AML_subg == "Deletion 5q")] = "Del(5q)"
  AML_subg[which(AML_subg == "KMT2A rearranged")] = "KMT2A-r"
  AML_subg[which(AML_subg == "NK triple mut")] = "NKt"
  AML_subg[which(AML_subg == "NPM1 mutated")] = "NPM1-mut"
  AML_subg[which(AML_subg == "RUNX1 mutated")] = "RUNX1-mut"
  AML_subg[which(AML_subg == "Monosomy 7")] = "Mono_7"
  AML_subg[which(AML_subg == "inv(16)")] = "inv(16)"
  return(AML_subg)
}



format_classification <- function(Classification) {
  Classification[which(Classification == "Granulocytic-UNK")] = "Mature myeloid lineage"
  Classification[which(Classification == "Immature-Neutrophil")] = "Mature myeloid lineage"
  Classification[which(Classification == "Committed Neutrophil")] = "Mature myeloid lineage"
  Classification[which(Classification == "Neutrophil")] = "Mature myeloid lineage"
  
  Classification[which(Classification == "Pre-Dendritic")] = "Myeloid dendritic cells"
  Classification[which(Classification == "Dendritic Cell_Lympho")] = "Myeloid dendritic cells"
  Classification[which(Classification == "Dendritic Cell_Myelo")] = "Myeloid dendritic cells"
  
  Classification[which(Classification == "CD34+ MEP")] = "Erythrocytic lineage"
  Classification[which(Classification == "CD34+ ERP")] = "Erythrocytic lineage"
  Classification[which(Classification == "CD34+ ERP-Early")] = "Erythrocytic lineage"
  Classification[which(Classification == "Early-Erythroblast")] = "Erythrocytic lineage"
  Classification[which(Classification == "Erythroblast")] = "Erythrocytic lineage"
  
  Classification[which(Classification == "CD34+ MKP")] = "Megakaryocytic lineage"
  Classification[which(Classification == "Platelet")] = "Megakaryocytic lineage"
  
  Classification[which(Classification == "CD34+ pre-PC")] = "B cell lineage"
  Classification[which(Classification == "Pro-B")] = "B cell lineage"
  Classification[which(Classification == "Follicular B cell")] = "B cell lineage"
  
  Classification[which(Classification == "CD34+ pre-T")] = "T/NK cell lineage"
  Classification[which(Classification == "Naive T-cell")] = "T/NK cell lineage"
  Classification[which(Classification == "CD8 T-cell")] = "T/NK cell lineage"
  Classification[which(Classification == "NK cells")] = "T/NK cell lineage"
  
  Classification[which(Classification == "CD14+ Monocyte")] = "Monocyte"
  Classification[which(Classification == "CD16+ Monocyte")] = "Monocyte"
  
  Classification[which(Classification == "cDC")] = "Dendritic cells"
  Classification[which(Classification == "pDC")] = "Dendritic cells"
  Classification[which(Classification == "Myeloid dendritic cells")] = "Dendritic cells"
  
  return(Classification)
}


factor_classification <- function(df_res, cellTypeGrouped) {
  if (cellTypeGrouped) {
    df_res$CellType = factor(df_res$CellType, levels=c(
      "CD34+ HSC", "CD34+ MultiLin", "CD34+ LMPP", "CD34+ Gran",
      "Mature myeloid lineage", "CD34+ Eo/B/Mast", "CD34+ MDP",
      "Dendritic cells", "Monocyte", "Erythrocytic lineage",
      "Megakaryocytic lineage", "CD34+ CLP", "B cell lineage",
      "Plasma Cell", "T/NK cell lineage", "Stromal"
    ))
  } else {
    df_res$CellType = factor(df_res$CellType, levels=c(
      "CD34+ HSC", "CD34+ MultiLin", "CD34+ LMPP", "CD34+ Gran",
      "Granulocytic-UNK", "CD34+ Eo/B/Mast", "CD34+ MDP", "CD14+ Monocyte",
      "CD16+ Monocyte", "Pre-Dendritic", "pDC", "cDC", "CD34+ MEP",
      "CD34+ ERP-Early", "CD34+ ERP", "Early-Erythroblast", "Erythroblast",
      "CD34+ MKP", "Platelet", "CD34+ CLP", "CD34+ pre-PC", "Pro-B",
      "Follicular B cell", "Plasma Cell", "CD34+ pre-T", "Naive T-cell",
      "CD8 T-cell", "NK cells", "Stromal"
    ))
  }
  
  return(df_res)
}

get_zscore <- function(expr) {
  return ((expr - mean(expr))/sd(expr))
}


plot_sc_plot <- function(gene.list, cellType, zscore, zscoreByGene=TRUE, cellTypeGrouped=TRUE, palette="YlOrRd") {
  
  AML_subg = scLoom[["col_attrs/AML_subgroup"]][]
  AML_subg = format_subgroup(AML_subg)
  
  if (cellType != "All") {
    Classification = scLoom[["col_attrs/Classification"]][]
    if (cellTypeGrouped) {
      Classification = format_classification(Classification)
    }
    
    index_cell = which(Classification %in% cellType)
    AML_subg = AML_subg[index_cell]
  }
  unique_AML_subg = unique(AML_subg)
  
  subAMLExp <- c()
  geneExps <- c()
  for (gene in gene.list){
    expr = as.vector(scLoom[["matrix"]][, scLoom[["row_attrs/Gene"]][] == gene])
    if (cellType != "All") {
      expr = expr[index_cell]
    }
    geneExps <- cbind(geneExps, expr)
    
    meanExprBySubgAML <- c()
    for (subg in unique_AML_subg) {
      index_subg = which(AML_subg == subg)
      meanExprBySubgAML = rbind(meanExprBySubgAML, mean(expr[index_subg]))
    }
    subAMLExp <- cbind(subAMLExp, meanExprBySubgAML)
  }
  
  colnames(subAMLExp) <- gene.list
  colnames(geneExps) <- gene.list
  
  build_prct <- function(subg) {
    row_selected = which(AML_subg == subg)
    exprSelected = geneExps[row_selected, ]
    if (is.null(dim(exprSelected))) {
      prct = ((length(exprSelected)-sum(exprSelected==0)) / length(exprSelected)) * 100
    } else {
      prct = ((length(exprSelected[,1])-colSums(exprSelected==0)) / length(exprSelected[,1])) * 100
    }
    
    return(data.frame(gene=gene.list, subgroup=subg, prct=prct))
  }
  
  row_or_col = 1 # 1: aplly by row (cellule)
  if(zscoreByGene) {
    row_or_col = 2 #2: aplly by col (gene)
  }
  if(zscore) {
    res = apply(subAMLExp, row_or_col, get_zscore)
  } else {
    res = subAMLExp
  }
  
  if(zscoreByGene) {
    colnames(res) = gene.list
    rownames(res) = unique_AML_subg
  } else {
    colnames(res) = unique_AML_subg
    rownames(res) = gene.list
  }
  df_res <- as.data.frame(res)
  df_res = df_res %>% 
    rownames_to_column("id") %>%
    gather(conditions, values, -id)
  
  if(zscoreByGene) {
    colnames(df_res) = c("subgroup", "gene", "avg")
  } else {
    colnames(df_res) = c("gene", "subgroup", "avg")
  }
  
  res_prct = lapply(unique_AML_subg, build_prct)
  res_prct = bind_rows(res_prct)
  
  res = merge(df_res, res_prct, by=c("gene", "subgroup"))
  res$gene <- factor(res$gene, levels = rev(gene.list))
  res$subgroup <- factor(res$subgroup, levels = c(
    "KMT2A-r", "inv(16)", "Mono_7", "Del(5q)", "NKt", "RUNX1-mut", "NPM1-mut", "CK"
  ))

  gg = ggplot(data = res)
  if (zscore) {
    gg = gg + geom_point(aes(y = gene, x = subgroup, fill = avg, size = prct), pch = 21)     
    gg = gg + labs(
      fill = "zscore"
    )
  } else {
    gg = gg + geom_point(aes(y = gene, x = subgroup, fill = avg, size = prct), pch = 21)   
    gg = gg + labs(
      fill = "average normalized\nexpression"
    )
  }
  if (zscore) {
    max_zscore = max(abs(min(res$avg)), max(res$avg))
    gg = gg + scale_fill_paletteer_c(
      paste("grDevices::", palette, sep=""),
      direction = -1,
      limits = c(-max_zscore,max_zscore)
    )
  } else {
    gg = gg + scale_fill_paletteer_c(
      paste("grDevices::", palette, sep=""),
      direction = -1,
    )
  }
  gg = gg + scale_radius(range=c(0,10), limits = c(1, 100), breaks=c(20,40,60,80,100))
  gg = gg + theme_bw()
  gg = gg + labs(
    size = "percent of positive cells"
  )
  if (cellType == "All") {
    gg = gg + ggtitle("All cell types")
  } else {
    gg = gg + ggtitle(paste(cellType, " (", length(index_cell), " cells)", sep=""))  
  }
  
  gg = gg + xlab("AML subgroups")
  # gg = gg + scale_fill_gradientn(colours=c("black", dot_color))
  gg = gg + theme(
    text = element_text(size=16),
    axis.text.x = element_text(angle=45,hjust=1),
    plot.title = element_text(hjust = 0.5)
  )
  
  if(zscore) {
    if (zscoreByGene) {
      gg = gg + theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line(size=.1, color="black") 
      )
    } else {
      gg = gg + theme(
        # remove the vertical grid lines
        panel.grid.major.x =  element_line(size=.1, color="black"),
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_blank() 
      )
    }
  } else {
    gg = gg + theme(
      panel.grid = element_blank(),
    )
  }
  
  return(gg)
}




plot_sc_heatmap = function(gene.list, sc_heatmap_subgroup, zscore, cellTypeGrouped=T, palette="YlOrRd") {
  Classification = scLoom[["col_attrs/Classification"]][]
  if (cellTypeGrouped) {
    Classification = format_classification(Classification)
  }
  
  AML_subg = scLoom[["col_attrs/AML_subgroup"]][]
  AML_subg = format_subgroup(AML_subg)
  if (sc_heatmap_subgroup != "All") {
    index_subg = which(AML_subg == sc_heatmap_subgroup)
    AML_subg = AML_subg[index_subg]
    Classification = Classification[index_subg]
  }
  uniqueCellType = unique(Classification)
  
  geneExps <- c()
  cellTypeExps <- c()
  for (gene in gene.list){
    expr = as.vector(scLoom[["matrix"]][, scLoom[["row_attrs/Gene"]][] == gene])
    if (sc_heatmap_subgroup != "All") {
      expr = expr[index_subg]
    }
    geneExps <- cbind(geneExps, expr) 
    
    meanExprCellType <- c()
    for (cellType in uniqueCellType) {
      index_cell = which(Classification == cellType)
      meanExprCellType = rbind(meanExprCellType, mean(expr[index_cell]))
    }
    cellTypeExps <- cbind(cellTypeExps, meanExprCellType)
  }
  
  colnames(cellTypeExps) <- gene.list
  rownames(cellTypeExps) <- uniqueCellType
  
  if(zscore) {
    res = apply(cellTypeExps, 2, get_zscore)
  } else {
    res = cellTypeExps
  }
  colnames(res) = gene.list
  rownames(res) = uniqueCellType
  
  
  df_res <- as.data.frame(res)
  df_res = df_res %>% 
    rownames_to_column("id") %>%
    gather(conditions, values, -id)
  colnames(df_res) = c("CellType", "gene", "avg")
  
  df_res$gene <- factor(df_res$gene, levels = rev(gene.list))
  df_res = factor_classification(df_res, cellTypeGrouped) 
  
  gg = ggplot(data=df_res)
  if (zscore) {
    gg = gg + geom_raster(aes(y=gene,x=CellType, fill=avg))
    gg = gg + labs(
      fill = "zscore",
    )
  } else {
    gg = gg + geom_raster(aes(y=gene,x=CellType, fill=avg))
    gg = gg + labs(
      fill = "average normalized\nexpression",
    )
  }
  gg = gg + theme_bw()
  if (sc_heatmap_subgroup == "All") {
    gg = gg + ggtitle("All AML subgroups")
  } else {
    gg = gg + ggtitle(paste(sc_heatmap_subgroup, " (", length(index_subg), " cells)", sep=""))  
  }
  gg = gg + theme(
    text = element_text(size=16),
    axis.text.x = element_text(angle=45,hjust=1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5))
  
  if (zscore) {
    max_zscore = max(abs(min(df_res$avg)), max(df_res$avg))
    gg = gg + scale_fill_paletteer_c(
      paste("grDevices::", palette, sep=""),
      direction = -1,
      limits = c(-max_zscore,max_zscore)
    )
  } else {
    gg = gg + scale_fill_paletteer_c(
      paste("grDevices::", palette, sep=""),
      direction = -1,
    )
  }
  
  
  return(gg)
}




