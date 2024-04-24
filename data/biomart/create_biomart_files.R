
file_biomart = file.path("ensembl_gene_trans.txt")
if(!file.exists(file_biomart)) {
  mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  tab_ens_gene = biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "gene_biotype"), mart = mart)
  tab_ens_gene = dplyr::rename(tab_ens_gene, ens_gene = ensembl_gene_id, ens_trans = ensembl_transcript_id, Gene.names=external_gene_name)
  write.table(tab_ens_gene, file=file_biomart, quote=FALSE, row.names=FALSE, sep="\t")
}

