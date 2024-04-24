library(data.table)

path_data = file.path("..")
path_biomart = file.path(path_data, "biomart")
file_biomart = file.path(path_biomart, "human", "ensembl_gene_trans.txt")

file_uni2ens = file.path("uni2ens.txt")
file_uni2ens_unique = file.path("uni2ens_unique.txt")
if(!file.exists(file_uni2ens)) {
  tab_ens_gene = fread(file_biomart, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
  file_tmp = file.path("HUMAN_9606_idmapping.dat")
  tab_uni2ens = fread(file_tmp, header=FALSE, stringsAsFactors=FALSE, sep="\t", quote = "")
  colnames(tab_uni2ens) = c("uniprot", "type", "trans_id")
  tab_uni2ens = tab_uni2ens[which(tab_uni2ens$type == "Ensembl_TRS" | tab_uni2ens$type == "Gene_Name"),]
  tab_uni2ens$uniprot = unlist(lapply(tab_uni2ens$uniprot, function(x) unlist(strsplit(x, "-"))[1]))
  tab_uni2ens = tab_uni2ens[!duplicated(tab_uni2ens[, c("uniprot", "type")]),]
  tab_uni2ens = dcast(tab_uni2ens, uniprot ~ type, value.var=c("trans_id"))
  colnames(tab_uni2ens) = c("uniprot", "ens_trans", "Gene.names")

  tab_uni2ens = merge(tab_uni2ens, tab_ens_gene, by=c("ens_trans", "Gene.names"), all.x=TRUE)
  write.table(tab_uni2ens, file=file_uni2ens, quote=FALSE, row.names=FALSE, sep="\t")
  tab_uni2ens_uniq_uniprot = tab_uni2ens[!duplicated(tab_uni2ens[, c("uniprot")]),]
  write.table(tab_uni2ens_uniq_uniprot, file=file_uni2ens_unique, quote=FALSE, row.names=FALSE, sep="\t")

}


