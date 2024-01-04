
compute_detection <- function(samples_ids, df_log_int) {
  detection = 1 - rowSums(df_log_int == 0) / length(samples_ids)
  
  return(detection)
}


test_diff <- function (se, type = c("control", "all", "manual"), control = NULL,
                       test = NULL, design_formula = formula(~0 + condition))
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type), class(design_formula) == "formula")
  type <- match.arg(type)
  col_data <- colData(se)
  raw <- assay(se)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)),
            "'")
  }
  if (!is.null(control)) {
    assertthat::assert_that(is.character(control), length(control) ==
                              1)
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"),
           "'", call. = FALSE)
    }
  }
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>%
    as.character() %>% .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control,
                                                    "- ", sep = " "), "", .) %>% paste(" - ", control,
                                                                                       sep = "")
      }
    }
  }
  if (type == "control") {
    if (is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"), "', for example '",
           paste0(conditions[1], "_vs_", conditions[2]),
           "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  }
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_", cntrst),
                                      collapse = ", "))
  fit <- limma::lmFit(raw, design = design)
  made_contrasts <- limma::makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- limma::contrasts.fit(fit, made_contrasts)
  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i, levels = design[,
                                                                      covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates],
                                           single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[,1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[,1]
    }
  }
  eB_fit <- limma::eBayes(contrast_fit)
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- limma::topTable(fit, sort.by = "t", coef = comp, number = Inf,
                           confint = TRUE)
    res <- res[!is.na(res$t), ]
    #TODO tempory fix for https://github.com/arnesmits/DEP/issues/7
    res$qval <- res$adj.P.Val
    #fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    #res$qval <- fdr_res$qval
    #res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- tibble::rownames_to_column(res)
    return(res)
  }
  limma_res <- purrr::map_df(cntrst, retrieve_fun)
  table <- limma_res %>% dplyr::select(rowname, logFC, CI.L, CI.R,
                                       P.Value, qval, comparison) %>% dplyr::mutate(comparison = gsub(" - ",
                                                                                                      "_vs_", comparison)) %>% tidyr::gather(variable, value, -c(rowname,
                                                                                                                                                                 comparison)) %>% dplyr::mutate(variable = recode(variable, logFC = "diff",
                                                                                                                                                                                                                  P.Value = "p.val", qval = "p.adj")) %>% tidyr::unite(temp, comparison,
                                                                                                                                                                                                                                                                       variable) %>% tidyr::spread(temp, value)
  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE)
  return(se)
}

delete_suffix <- function(words) {
  # Get prefix
  suffix <- get_suffix(words)
  # Delete prefix from words
  gsub(paste0(suffix, "$"), "", words)
}



get_detection <- function(df_ms, exp_design, col_as_condition, subgroup, overview_max_sd, overview_min_pep_uniq, overview_min_spat) {
  df_ms = df_ms[which(df_ms["#Unique"] >= overview_min_pep_uniq), ]
  df_ms = df_ms[which(df_ms["SPAT_score"] >= overview_min_spat), ]
  
  samples_ids = exp_design$bclq_id[which(exp_design[col_as_condition] != 'None')]
  if(subgroup != "All") {
    samples_ids = which(exp_design[col_as_condition] == subgroup)
    samples_ids = exp_design$bclq_id[samples_ids]
  }
  samples_ids = which(colnames(df_ms) %in% samples_ids)
  
  df_log_int = log2(df_ms[, samples_ids]+1)
  detection = compute_detection(samples_ids, df_log_int)
  
  # We want calculate sd and mean on 'positive' sample
  df_log_int[df_log_int == 0] = NA
  sd_prot_intensity = apply(df_log_int, 1, function(x) sd(x, na.rm = TRUE))
  mean_prot_intensity = apply(df_log_int, 1, function(x) mean(x, na.rm = TRUE))
  max_prot_intensity = apply(df_log_int, 1, function(x) max(x, na.rm = TRUE))
  
  # mean and sd function return NA whent the protein is detected in only one sample  
  sd_prot_intensity[which(is.na(sd_prot_intensity))] = 0
  mean_prot_intensity[which(is.na(mean_prot_intensity))] = max_prot_intensity[which(is.na(mean_prot_intensity))]
  mean_prot_intensity[which(mean_prot_intensity == -Inf)] = 0
  
  df_ms$mean_log2_int = mean_prot_intensity
  df_ms$sd_log2_int = sd_prot_intensity
  df_ms$detection_p = detection * 100
  
  
  if(subgroup != "All") {
    other_ids = which(!(subgroup %in% exp_design[col_as_condition]))
    other_bclq_id = exp_design$bclq_id[other_ids]
    other_ids = which(colnames(df_ms) %in% other_bclq_id)
    df_ms = df_ms[, -other_ids]
  }
  
  df_ms = df_ms[which(df_ms$sd_log2_int <= overview_max_sd),]
  
  return(df_ms)
}


fill_replicate_surface <- function(condition, exp_design) {
  res = exp_design[which(exp_design$condition == condition),]
  res$replicate = c(1:length(res$replicate))
  return(res)
}

run_DEP <- function(df_ms, exp_design, dep_select_condition, dep_select_comparison, 
                    imputation, min_peptide, df_anno, alpha = 0.05,
                    lfc = 0, num_missval=0) {
  
  if (dep_select_comparison != 'Compare all condition pairs') {
    dep_select_comparison = strsplit(dep_select_comparison, ' ')[[1]][1]
  }
  
  param <- data.frame(alpha, lfc)
  
  exp_design = df_design
  index_condition = which(colnames(exp_design) == dep_select_condition)
  colnames(exp_design)[index_condition] = "condition"
  if (dep_select_comparison != "Compare all condition pairs") {
    exp_design$condition[which(exp_design$condition != dep_select_comparison)] = "Other"
  }
  
  conditions = unique(exp_design$condition)
  exp_design = lapply(conditions, fill_replicate_surface, exp_design)
  exp_design = rbindlist(exp_design)
  setDF(exp_design)
  exp_design = exp_design[order(exp_design$track),]
  
  exp_design$label = exp_design$bclq_id
  
  index_unique = which(colnames(df_ms) == "#Unique")
  df_ms = df_ms[which(df_ms[,index_unique] >= min_peptide), ]
  
  col_samples = which(colnames(df_ms) %in% df_design$bclq_id)
  del_samples = col_samples[-exp_design$track]
  if (length(del_samples) != 0) {
    df_ms = df_ms[,-del_samples]
  }
  col_samples = which(colnames(df_ms) %in% exp_design$bclq_id)
  
  colnames(df_ms)[col_samples] = paste("Intensity", exp_design$bclq_id, sep=".")
  
  df_raw_int = df_ms
  col_samples = which(colnames(df_raw_int) %in% paste("Intensity", exp_design$bclq_id, sep="."))
  colnames(df_raw_int)[col_samples] = exp_design$bclq_id

  # Data preparation: make unique ID
  df_ms_uniq = make_unique(df_ms, "uniprot", "uniprot", delim = ";")
  
  # Generate a SummarizedExperiment object
  intensity_columns = grep("Intensity.", colnames(df_ms_uniq))
  df_ms_uniq[,intensity_columns] = sapply(df_ms_uniq[,intensity_columns], as.numeric)
  
  df_ms_se <- make_se(df_ms_uniq, intensity_columns, exp_design)
  
  # Filter for proteins that are identified in all replicates of at least one condition
  df_ms_filt <- filter_missval(df_ms_se, thr = num_missval)
  print(head(df_ms_filt))
  
  # Normalize the data
  df_ms_norm <- normalize_vsn(df_ms_filt)
  
  # All possible imputation methods are printed in an error, if an invalid function name is given.
  # impute(df_ms_norm, fun = "")
  
  # Impute missing data
  # TODO q = 0.01 for some imputation functions
  df_ms_imputed = NULL
  tryCatch(
    expr = {
      df_ms_imputed = impute(df_ms_norm, fun=imputation)
    },
    error = function(e) {
      print('!!! ERROR: Impute')
      print(e)
    }
  )
  imputed_error = F
  if (is.null(df_ms_imputed)) {
    print("Error during imputation")
    df_ms_imputed = impute(df_ms_norm, fun='zero')
    imputed_error = T
  }
  
  # Test every sample versus control
  if (dep_select_comparison != "Compare all condition pairs") {
    cond = unique(exp_design$condition[which(exp_design$condition != "Other")])
    df_ms_diff <- test_diff(df_ms_imputed, type = "manual", test=paste(cond, "vs", "Other", sep="_"))
  } else {
    df_ms_diff <- test_diff(df_ms_imputed, type = "all")
  }
  
  # Denote significant proteins based on user defined cutoffs
  dep <- add_rejections(df_ms_diff, alpha=alpha, lfc=lfc)
  
  df_dep = dep %>% rowData() %>% as.data.frame() %>% dplyr::select(ends_with("_diff"))
  df_dep$uniprot = rownames(df_dep)
  df_dep = merge(df_dep, tab_uni2ens, by="uniprot")
  proteo_biotype = unique(df_dep$gene_biotype)
  
  df_raw_int = merge(df_raw_int, tab_uni2ens, by="uniprot")
  
  df_raw_int$id = paste(df_raw_int$Gene.names, " (", df_raw_int$uniprot, ")", sep="")
  
  df_dep_res = list(dep=dep, df_ms=df_ms, df_raw_int=df_raw_int,
                   filt=df_ms_filt, norm=df_ms_norm, imputed=df_ms_imputed,
                   param=param, samples=df_design, design=exp_design,
                   biotype=proteo_biotype, imputed_error=imputed_error)
  
  return(df_dep_res)
}

get_dep_table <- function(df_dep_res, df_anno, type, significant=FALSE) {
  res = NA
  if (type == "filtred") res = as.data.frame(assay(df_dep_res$filt))
  if (type == "normalized") res = as.data.frame(assay(df_dep_res$norm))
  if (type == "imputed") res = as.data.frame(assay(df_dep_res$imputed))
  
  colnames(res) = df_dep_res$design$label
  res$uniprot = rownames(res)
  
  if(significant) {
    df_dep = df_dep_res$dep %>% .[rowData(.)$significant, ] %>% rowData() %>% data.frame()
    df_dep$uniprot = rownames(df_dep)
    res = res[which(res$uniprot %in% df_dep$uniprot),]
  }
  
  df_ms = df_dep_res$df_ms
  col_name = c("uniprot", "Coverage (%)", "#Peptides",	"#Unique")
  col_sel = unlist(lapply(col_name, function(x) which(colnames(df_ms) == x)))
  if (length(col_sel) > 1) {
    df_ms = df_ms[,col_sel]
    colnames(df_ms) = c("uniprot", "coverage", "#pep", "#uniq")
    res = merge(res, df_ms, by="uniprot", all.x=TRUE)
  } else {
    res = res[which(res$uniprot %in% df_ms$uniprot), ]
  }
  
  res = merge(res, df_anno, by="uniprot", all.x=TRUE)
  
  return(res)
  
}

get_gene_exp <- function(trans_subgoups, gene_exp, tab_design, gene_name, graphic_type, sel_exp) {
  
  if (is.null(gene_exp)) return(data.frame())
  if (gene_name == "") return(data.frame())
  
  row_gene_id = which(gene_exp$Gene.names == gene_name | gene_exp$ID2 == gene_name)
  gene_name = gene_exp$Gene.names[row_gene_id]
  
  gene_exp = gene_exp[row_gene_id, ]
  gene_exp = reshape2::melt(gene_exp, id.vars = c("ID", "ID2", "Gene.names"), stringsAsFactors=FALSE)
  
  colnames(gene_exp)[4] = "sample"
  if (sel_exp == "readcounts.xls") {
    colnames(gene_exp)[5] = "log2_CPM"
    gene_exp$log2_CPM = log2(as.numeric(gene_exp$log2_CPM) + 1)
  } else {
    colnames(gene_exp)[5] = "log2_TPM"
    gene_exp$log2_TPM = log2(as.numeric(gene_exp$log2_TPM) + 1)
  }
  gene_exp = merge(gene_exp, tab_design, by="sample")
  
  gene_exp$group[which(gene_exp$group == "Query")] = trans_subgoups[1]
  if(length(trans_subgoups) == 2) {
    gene_exp$group[which(gene_exp$group == "Ref")] = trans_subgoups[2]
    gene_exp$group = factor(gene_exp$group, levels= c(trans_subgoups[1], trans_subgoups[2]))
  } else {
    gene_exp$group[which(gene_exp$group == "Ref")] = "other"
    gene_exp$group = factor(gene_exp$group, levels= c(trans_subgoups[1], "other"))
  }
  
  return(gene_exp)
}

update_all_gene_exp <- function(all_gene_exp) {
  
  count_only = all_gene_exp[,-which(colnames(all_gene_exp) == "ID")]
  count_only = apply(count_only, 2,
                     function(x) {
                       total_count = sum(x)
                       count_norm = (x / total_count) * 1e6
                     }
  )
  df_cpm_anno = cbind(data.frame(ID=all_gene_exp$ID, stringsAsFactors=F), count_only)
  
  df_cpm_anno$ID2 = unlist(lapply(df_cpm_anno$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  df_cpm_anno = merge(df_cpm_anno, tab_ens_gene, by.x="ID2", by.y="ens_gene", all.x=TRUE)
  return(df_cpm_anno)
}

melt_df_imputed <-function(df_imputed, df_dep_res, dep_select_gene, tab_list_genes_dep) {
  
  gene_name = dep_select_gene
  row_gene_id = which(tab_list_genes_dep$Gene.names == gene_name | tab_list_genes_dep$uniprot == gene_name | tab_list_genes_dep$ens_gene == gene_name | tab_list_genes_dep$ens_trans == gene_name)
  gene_name = tab_list_genes_dep$uniprot[row_gene_id]
  
  df_imputed = df_imputed[which(df_imputed$uniprot == gene_name),]
  
  design = df_dep_res$design
  col_2_melt = which(colnames(df_imputed) %in% design$bclq_id)
  df_imputed = reshape2::melt(df_imputed, id.vars = colnames(df_imputed)[-col_2_melt])
  df_imputed <- plyr::rename(df_imputed,c('variable'='bclq_id'))
  df_imputed <- plyr::rename(df_imputed,c('value'='log2(Normalized Imputed Intensity + 1)'))
  df_imputed$imputed = as.vector(t(is.na(assay(df_dep_res$norm[gene_name]))))
  df_imputed = merge(df_imputed, design, by="bclq_id")
  
  return(df_imputed)
}

## TO DO, figure out how to print a message if 
## k is too big
get_scKNN = function(gene_name,k){ 
  #kNN_error = FALSE
  if (k > ncol(scGenesNN)){
    #return (paste("Only ",ncol," neighbors are computed",sep=""))
    k = ncol(scGenesNN)
    #kNN_error = TRUE
  }
    return (rownames(scGenesNN)[unlist(scGenesNN[gene_name,1:k])])
}

get_corGenes = function(gene_name,k){ 
  return (names(sort(scGeneCorrelation[gene_name,],decreasing = T)[1:k]))
}

order_genes = function(geneSet,geneExp){
  return(hclust(dist(geneExp[geneSet,])))
}


