#Rscript --vanilla path_design path_counts
#ex:
# - [path to design]
# - [path to count]
# - [path to save output file]


library(limma)
library(edgeR)

args = commandArgs(trailingOnly=TRUE)

path_design=args[1]
path_counts=args[2]
path_output=args[3]

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

file_out = file.path(path_output, "limma_voom_genes.xls")

design_file = file.path(path_design, "design.tsv")
design = read.table(design_file, header = TRUE, stringsAsFactors=FALSE)

file_counts = file.path(path_counts, "readcounts.xls")
counts = read.delim(file_counts, header=TRUE, row.names=1)
colnames(counts) = gsub("^X", "", colnames(counts))
design$sample = gsub("-", ".", design$sample)
counts = counts[, which(colnames(counts) %in% design$sample)]
counts = trunc(counts)
counts = counts[apply(counts, 1, function(x) !all(x==0)),]
design = design[which(design$sample %in% colnames(counts)), ]
if (length(unique(design$group)) >= 2) {
  sampleTable = data.frame(condition = design$group)

  #rownames(counts) = unlist(lapply(rownames(counts), function(x) unlist(strsplit(x, "[.]"))[1]))

  counts = counts[, match(design$sample, colnames(counts))]

  dge <- DGEList(counts=counts, group=design$group)

  #keep <- filterByExpr(dge)
  #dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

  design = model.matrix(~condition, data = sampleTable)
  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  #To give more weight to FC in the ranking
  #fit <- treat(fit, lfc=log2(1.2))
  #topTreat(fit, coef=ncol(design))
  res = topTable(fit, n=Inf, sort.by="p")
  res = cbind(data.frame(ID=row.names(res), stringsAsFactors=FALSE), res)
  write.table(res, file_out, quote=FALSE, row.names=FALSE, sep="\t")
}
