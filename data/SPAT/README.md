Small documentation to generate these data

```
# In this directory

# 1. Download HUMAN_9606_idmapping.dat from uniprot 

# 2. Generare all_uniprot_ids.txt file
cut -f1 ./HUMAN_9606_idmapping.dat | sort | uniq > all_uniprot_ids.txt

# 3. Generate SPAT_annotation_uniprot_id.txt file on SPAT shiny web portal (https://spat.leucegene.ca/), using all_uniprot_ids.txt

# 4. Edit SPAT_annotation_uniprot_id.txt to keep uniport id and SPAT score columns only

# 5. Rename these two columsn in: uniprot_id and SPAT_score

# 6. Edit create_spat_table.r to change version

# 7. Create 2 files used by LASA (shiny server)
Rscript ./create_spat_table.r

# 8. In root directory, edit shiny/www/lasa/global.R script to update spat_version variable
```

