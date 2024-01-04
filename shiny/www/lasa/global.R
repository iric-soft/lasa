library(ggrepel)
require(ggplot2)
library(cowplot)
library(ggbiplot)
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library(DEP)
library(shiny)
library(shinydashboard)
library(data.table)
library(shinyWidgets)
library(shinyjs)
library(here)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(stringi)
library(tidyr)
library(Laurae)
library(umap)
library(tableHTML)
library(plotly)
library(loomR)
library(Dict)
library(paletteer)

# On big screen
#plot_height = c(800, 590, 750, 490)
# On laptop
plot_height = c("700px", "500px", "600px", "300px", "2000px")
script_dir = here::here()
spat_version = "2021-10-21"

path_data = file.path(script_dir, "data")
path_surface = file.path(path_data, "surfaceome")
path_globale = file.path(path_data, "globaleome")
path_spat = file.path(path_data, "SPAT", spat_version)
path_uniprot = file.path(path_data, "uniprot")
path_biomart = file.path(path_data, "biomart")
path_diff = file.path(path_data, "transcriptome")
path_sc = file.path(path_data,"singleCell")
path_project = path_surface

path_logo = file.path(script_dir, "shiny", "www", "lasa", "LASA_logo.jpg")

tab_samples = NA

#file_biomart = file.path(path_biomart, "ensembl_gene_trans.txt")
#tab_ens_gene_trans = fread(file_biomart, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

file_uni2ens=file.path(path_uniprot,"uni2ens.tsv")
tab_uni2ens = fread(file_uni2ens, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
tab_uni2ens = as.data.frame(tab_uni2ens)

file_uni2ens_uniq=file.path(path_uniprot,"uni2ens_uniq.tsv")
tab_uni2ens_uniq_uniprot <<- fread(file_uni2ens_uniq, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
tab_uni2ens_uniq_uniprot <<- as.data.frame(tab_uni2ens_uniq_uniprot)

file_in = file.path(path_spat, "SPAT_score_uniprot_id.txt")
surface_score_uniprot = fread(file_in, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
df_anno = merge(tab_uni2ens_uniq_uniprot, surface_score_uniprot, by.x="uniprot", by.y="uniprot_id", all.x=TRUE)

custom.settings = umap.defaults
custom.settings$random_state = 42
custom.settings$n_neighbors = 10
custom.settings$n_epochs = 1000
custom.settings$min_dist = 0.2

file_design = file.path(path_project, 'lasa_paper', "samples.tsv")
df_design <<- read.table(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

file_ms = file.path(path_project, 'lasa_paper', "Reviewed", "proteins.csv")
df_ms <<- read.csv(file_ms, header = TRUE, stringsAsFactors=FALSE, na.strings = c("NA","-"), check.names = F)

list_condition <<- c("batch", "subgroup", "FAB")
list_subgroup <<- unique(df_design$subgroup)

cond = unique(df_design$subgroup)
cond = paste(cond, ' vs other', sep='')
list_comparison <<- c("All pairs of comparison", cond)

list_by <<- colnames(df_design)
remove_name = c("bclq_id", "track", "batch", "replicate", "mg")
list_by <<- list_by[which(!list_by %in% remove_name)]
cond = paste0(df_design$subgroup, ' vs other')
list_on <<- c("All pairs of comparison", cond)

df_ms = merge(df_ms, df_anno, by="uniprot", all.x=TRUE)

trans_list_comparison = list.dirs(file.path(path_diff, "design", 'leucegene2_430'), full.names=FALSE, recursive=F)

file_biomart = file.path(path_biomart, "ensembl_gene.txt")
tab_ens_gene = fread(file_biomart, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

file_in = file.path(path_spat, "SPAT_scores.ensembl_id.txt")
trans_surface_score = fread(file_in, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

trans_surface_score = merge(trans_surface_score, tab_ens_gene, by.x="ENSEMBL_id", by.y="ens_gene", all=T)

# single cell

# file_in = file.path(path_sc, "intensityByGroup.csv")
# df_sc_intensity <<- read.csv(file_in, header = TRUE, stringsAsFactors=FALSE, na.strings = c("NA","-"), check.names = F)

file_in = file.path(path_sc, "sc_samples.loom")
scLoom <<- connect(filename=file_in, mode = "r", skip.validate = TRUE)


levelCellTypeGrouped = c(
  "CD34+ HSC", "CD34+ MultiLin", "CD34+ LMPP", "CD34+ Gran",
  "Mature myeloid lineage", "CD34+ Eo/B/Mast", "CD34+ MDP",
  "Dendritic cells", "Monocyte", "Erythrocytic lineage",
  "Megakaryocytic lineage", "CD34+ CLP", "B cell lineage",
  "Plasma Cell", "T/NK cell lineage", "Stromal"
)

color_list_grouped  <<- Dict$new(
  'CD34+ HSC' = '#5A3AAB', 'CD34+ MultiLin' = '#BF7EED',   
  'CD34+ Gran' = '#08306B', 'Mature myeloid lineage' = '#08519C', 
  'Monocyte' = '#26A657',
  'Dendritic cells' = '#89F5D1', 'CD34+ MDP' = '#71E3BD', 
  'CD34+ LMPP' = "#C83DE0",  'CD34+ CLP' = "#CF5BE3", 
  'B cell lineage' = '#FFC2CD', 
  'CD34+ Eo/B/Mast' = '#575f88',
  'Erythrocytic lineage' = "#DC143C", 'Megakaryocytic lineage' = '#ed9393',
  'T/NK cell lineage' = '#FDAE6B',  
  'Plasma Cell' = '#C71585', 'Stromal' = "#ffd000"
)

color_palette_grouped <<- as.list(color_list_grouped$values)
names(color_palette_grouped) <- color_list_grouped$keys



levelCellType = c(
  "CD34+ HSC", "CD34+ MultiLin", "CD34+ LMPP", "CD34+ Gran",
  "Granulocytic-UNK", "CD34+ Eo/B/Mast", "CD34+ MDP", "CD14+ Monocyte",
  "CD16+ Monocyte", "Pre-Dendritic", "pDC", "cDC", "CD34+ MEP",
  "CD34+ ERP-Early", "CD34+ ERP", "Early-Erythroblast", "Erythroblast",
  "CD34+ MKP", "Platelet", "CD34+ CLP", "CD34+ pre-PC", "Pro-B",
  "Follicular B cell", "Plasma Cell", "CD34+ pre-T", "Naive T-cell",
  "CD8 T-cell", "NK cells", "Stromal"
)

color_list  <<- Dict$new(
  'CD34+ HSC'  = '#5A3AAB',  'CD34+ MultiLin'  = '#BF7EED',   
  'CD34+ Gran'  = '#08306B',  'Granulocytic-UNK'  = '#08519C', 
  'CD14+ Monocyte'  = '#26A657', 'CD16+ Monocyte'  = '#1C7A3F',
  'Pre-Dendritic'  = '#11BA82',  'cDC'  = '#4ED9AB',  'pDC'  = '#89F5D1',  'CD34+ MDP'  = '#71E3BD', 
  'CD34+ LMPP'  = "#C83DE0",  'CD34+ CLP'  = "#CF5BE3",  'CD34+ pre-PC'  = '#DB309C', 
  'Pro-B'  = '#FFC2CD',  'Follicular B cell'  =  '#FF69B4', 
  'CD34+ MEP'  = '#78222E', 'CD34+ Eo/B/Mast'  = '#575f88',
  'CD34+ ERP-Early'  = '#6E0202', 'CD34+ ERP'  = '#A11212', 
  'Early-Erythroblast'  = '#FF0000',  'Erythroblast'  = "#DC143C", 'CD34+ MKP'  = '#db7676',  'Platelet'  = '#ed9393',
  'CD34+ pre-T'  = '#D94801','Naive T-cell'  = '#F16913',  'CD8 T-cell'  =  '#FD8D3C',  'NK cells'  = '#FDAE6B',  
  'Plasma Cell'  = '#C71585',  'Stromal'  = "#ffd000"
)

color_palette <<- as.list(color_list$values)
names(color_palette) <- color_list$keys





sample_colors <<- Dict$new(
  '06H088' ='#c5fe74', #KMT2A-r
  '05H066' ='#85ea46', #KMT2A-r
  '02H017' ='#3cd672', #KMT2A-r
  '09H032' ='#0fec59', #KMT2A-r
  '09H010' ='#04b135', #KMT2A-r
  
  '18H047' ='#0474b1', #CK
  '12H138' ='#098fd7', #CK
  '12H106' ='#7fc9f1', #CK
  '11H103' ='#3ac5e8', #CK
  
  '14H007' ='#a53ae8', #NK triple-mut
  '12H010' ='#c47ff2', #NK triple-mut
  '07H134' ='#ed88d3', #NK triple-mut
  '09H070' ='#f50cba', #NPM1-mut
  
  '11H097' ='#585554', #del(5q)
  '17H065' ='#aa9f9b', #del(5q)
  '05H034' ='#eecdc2', #del(5q)
  '04H096' ='#d9d5d4', #mono 7
  '05H193' ='#c19587', #mono 7
  
  '08H087' ='#f3960c', #RUNX1-mut
  '03H112' ='#f46e31'  #inv(16)
)
sample_palette <<- as.list(sample_colors$values)
names(sample_palette) <- sample_colors$keys

AML_colors <<- Dict$new(
  'KMT2A rearranged'  ='#3cd672', #KMT2A-r '06H088','05H066','02H017','09H032','09H010'
  'Complex karyotype' ='#098fd7', #CK '18H047','12H138','12H106','11H103'
  'NK triple mut'     ='#c47ff2', #NK triple-mut  '14H007','12H010','07H134'
  'NPM1 mutated'      ='#f50cba', #NPM1-mut '09H070'
  'Deletion 5q'       ='#aa9f9b', #del(5q) '11H097','17H065','05H034'
  'Monosomy 7'        ='#c19587', #mono 7 '04H096','05H193'
  'RUNX1 mutated'     ='#f3960c', #RUNX1-mut
  'inv(16)'           ='#f46e31' #inv(16)
)
aml_palette <<- as.list(AML_colors$values)
names(aml_palette) <- AML_colors$keys

