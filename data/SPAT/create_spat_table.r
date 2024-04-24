
# version = "2021-10-21" # using SPAT 3.3
version = "2023-02-02"
path_data = "./"

file = file.path(path_data, version, "SPAT_annotation.tsv")
SPAT_uniprot <- read.delim(file)
SPAT_uniprot$uniprot_id = vapply(strsplit(SPAT_uniprot$uniprot_id,"\\-"), `[`, 1, FUN.VALUE=character(1))
file = file.path(path_data, version, "SPAT_annotation.tmp.tsv")
write.table(SPAT_uniprot, file=file, row.names = F, sep=",", quote=F)



file = file.path(path_data, version, "SPAT_annotation.tmp.tsv")
SPAT <- read.delim(file, sep=",")
file = file.path(path_data, "all_uniprot_ids.txt")
list_uniprot <- read.delim(file, header=F)

#Need data in SPAT git repo (the link will be available when SPAT is published)
# TODO update this path according to SPAT version
file = "./SPAT/SPAT3.4/database/translator.txt"
trans <- read.delim(file, header = F)

res1 = merge(SPAT, trans, by.x="uniprot_id", by.y="V2")
res1 = res1[, -c(1,3,4)]
res1 = res1[!duplicated(res1), ]

res2 = merge(list_uniprot, trans, by.x="V1", by.y="V2")
res2 = res2[, -c(2,3)]

res = merge(res1, res2, by.x="V4", by.y="V4", all=T)
res[,1] = res[,3]
res = res[,-3]
colnames(res) = c("uniprot_id", "SPAT_score")

res = res[-which(is.na(res$SPAT_score)), ]
res = res[!duplicated(res), ]
file = file.path(path_data, version, "SPAT_score_uniprot_id.txt")
write.table(res, file=file, quote=F, sep="\t", row.names=F)




res1 = merge(SPAT, trans, by.x="uniprot_id", by.y="V2")
res1 = res1[, -c(1,3,4)]
res1 = res1[!duplicated(res1), ]

res2 = merge(list_uniprot, trans, by.x="V1", by.y="V2")
res2 = res2[, -c(2,3)]

res = merge(res1, res2, by.x="V4", by.y="V4", all=T)
res[,1] = res[,3]
res = res[,-3]
colnames(res) = c("uniprot_id", "SPAT_score")

res = merge(res, trans, by.x="uniprot_id", by.y="V2")
res[,1] = res[,4]
res = res[,-c(3, 4, 5)]
res = res[!duplicated(res), ]
colnames(res)[1] = "ENSEMBL_id"

file = file.path(path_data, version, "SPAT_scores.ensembl_id.txt")
write.table(res, file=file, quote=F, sep="\t", row.names=F)
