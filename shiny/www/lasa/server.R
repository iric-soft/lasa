
source('global.R')
source('functions.R')
source('fig_functions.R')

shinyServer(function(input, output, session) {
  disable("downloadDepTable")
  
  updateSelectizeInput(
    session, 'sc_gene_names', 
    choices = scLoom[["row_attrs/Gene"]][], 
    selected = c("PTPRC", "CD38", "LAIR1", "CD47", "CD33", "CD37", "IL3RA", "FLT3", "CLEC12A", "ITGA4", "VSIR", "CD74"),
    server = TRUE
  )
  
  updateSelectizeInput(
    session, 'sc_track_gene_names', 
    choices = scLoom[["row_attrs/Gene"]][], 
    selected = c("CD34", 'CD33','IL3RA'),
    server = TRUE
  )
  
  updateSelectizeInput(
    session, 'sc_umap_gene_names', 
    choices = scLoom[["row_attrs/Gene"]][], 
    selected = c("CD34", 'CD33','IL3RA'),
    server = TRUE
  )
  
  updateSelectizeInput(
    session, 'sc_umap_gene_selected', 
    choices = scLoom[["row_attrs/Gene"]][], 
    server = TRUE
  )
  
  updateSelectizeInput(
    session, 'sc_gene_names_cell_type', 
    choices = scLoom[["row_attrs/Gene"]][], 
    selected = c("CD34", 'CD33','IL3RA'),
    server = TRUE
  )
  

  output$logo <- renderImage({
    list(src = path_logo,
         alt = "LASA logo")
  }, deleteFile = FALSE)

  output$downloadMsRaw <- downloadHandler(
    file = "leucegene_surfaceome_maxquant_raw_output.tsv",
    content = function(file) {
      withProgress(message = 'Prepare data...', value = 0.80, {
        write.table(df_ms, file, row.names = FALSE, sep="\t", quote=FALSE)
      })
    }
  )

  output$downloadSurfaceomeDesign <- downloadHandler(
    file = "leucegene_surfaceome_design.tsv",
    content = function(file) {
      withProgress(message = 'Prepare data...', value = 0.80, {
        write.table(df_design[, -c(2,5,6)], file, row.names = FALSE, sep="\t", quote=FALSE)
      })
    }
  )
  
  output$downloadGlobaleome <- downloadHandler(
    file = "globaleome.zip",
    content = function(file) {
      withProgress(message = 'Prepare data...', value = 0.80, {
        utils::zip(file, files=list_globaleome, extras = '-j')
      })
    }
  )

  output$downloadReadcounts <- downloadHandler(
    file = "leucegene_transcriptome_readcounts.tsv",
    content = function(file) {
      withProgress(message = 'Prepare data...', value = 0.80, {
        file_exp = file.path(path_diff, "STAR_RSEM_unstranded_v32", 'leucegene2_430', "readcounts.xls")
        file.copy(file_exp, file)
      })
    }
  )

  output$downloadTanscriptomeDesign <- downloadHandler(
    file = "leucegene_transcriptome_design.tsv",
    content = function(file) {
      withProgress(message = 'Prepare data...', value = 0.80, {
        design_dir = file.path(path_diff, "design", 'leucegene2_430')
        file_design = file.path(design_dir, "all_design.csv")
        file.copy(file_design, file)
      })
    }
  )
  
  #output$downloadSingleCellDesign <- downloadHandler(
  #  file = "leucegene_singleCell_design.tsv",
  #  content = function(file) {
  #    withProgress(message = 'Prepare data...', value = 0.80, {
  #      file_design = file.path(path_sc, "single_cell_design.tsv")
  #      file.copy(file_design, file)
  #    })
  #  }
  #)
  
  output$loomObject <- renderUI({
    url <- a("Single cell Loom file", href="https://lasa.leucegene.ca/sc/sc_samples.loom")
    tagList("Click on link to download loom sc file: ", url)
  })
  
  update_list_comparison <- function(input) {
    cond = unique(df_design$subgroup)
    cond = paste(cond, ' vs other', sep='')
    list_comparison <<- c("Compare all condition pairs", cond)
    cond = unique(df_design[,input$dep_select_condition])
    cond = paste(cond, ' vs other', sep='')
    list_comparison <<- c("Compare all condition pairs", cond)
    selected = ifelse(length(list_comparison) > 0, list_comparison[1], "")
    updateSelectInput(session, "dep_select_comparison", choices=list_comparison, selected=selected)
  }

  update_list_subgroup <- function(input) {
    uniq_subgroup = unique(df_design[input$overview_select_condition])
    uniq_subgroup = uniq_subgroup[uniq_subgroup != 'None']
    list_subgroup <<- c("All", uniq_subgroup)
    selected = ifelse(length(list_subgroup) > 0, list_subgroup[1], "")
    updateSelectInput(session, "overview_select_subgroup", choices=list_subgroup, selected=selected)
  }

  observeEvent(input$dep_select_condition, {
    withProgress(message = 'Update list of comparaison...', value = 0.80, {
      update_list_comparison(input)
    })
  })

  observeEvent(input$overview_select_condition, {
    withProgress(message = 'Update list of condition...', value = 0.80, {
      update_list_subgroup(input)
    })
  })


  overview_detection <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      if (input$overview_select_subgroup == ""){
        ggplot()
      } else {
        df_detection = get_detection(df_ms, df_design,
                               input$overview_select_condition,
                               input$overview_select_subgroup,
                               input$overview_max_sd,
                               input$overview_min_pep_uniq,
                               input$overview_min_spat
        )
        plot_detection(df_detection, input$overview_select_subgroup,
                    input$overview_min_pep_uniq,
                    input$overview_min_spat,
                    input$overview_max_sd,
                    input$overview_select_color,
                    input$overview_select_alpha
        )
      }
    })
  })

  output$overview_detection <- renderPlot({
    overview_detection()
  })

  output$downloadOverviewFig <- downloadHandler(
    filename = "percentage_presence.pdf",
    content = function(file) {
      pdf(file, useDingbats=FALSE)
      print(overview_detection())
      dev.off()
    }
  )
  
  
  sc_plot <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      if (length(input$sc_gene_names) == 0){
        ggplot()
      } else {
#        if (input$scPlotZscore) {
#          zscoreByGene <- switch(
#            input$scPlotZscoreGene,
#            gene = TRUE,
#            subg = FALSE
#          )
#          zscoreOnMeanSubg = T
          #zscoreOnMeanSubg <- switch(
          #  input$scPlotZscoreMean,
          #  mean = TRUE,
          #  raw = FALSE
          #)
#        } else {
#          zscoreByGene = T 
#          zscoreOnMeanSubg = T 
#        }
        
        plot_sc_plot(
          input$sc_gene_names,
          input$scCellTypeSelected,
          input$scPlotZscore,
#          zscoreByGene,
#          zscoreOnMeanSubg
          palette = input$scPaletteColor
        )
      }
    })
  })
  
  output$sc_plot <- renderPlot({
    sc_plot()
  })
  
  output$downloadScPlotFig <- downloadHandler(
    filename = "singleCellBySubGroup.pdf",
    content = function(file) {
      pdf(file, width=10, height=(1+length(input$sc_gene_names)*0.5), useDingbats=FALSE)
      print(sc_plot())
      dev.off()
    }
  )
  
  
  sc_heatmap <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      if (length(input$sc_gene_names_cell_type) == 0){
        ggplot()
      } else {
        plot_sc_heatmap(
          input$sc_gene_names_cell_type, input$sc_heatmap_subgroup,
          input$scHeatmapZscore,
          input$scHeatmapCellType,
          palette = input$scHeatmapPaletteColor
        )
      }
    })
  })
  
  output$sc_heatmap <- renderPlot({
    sc_heatmap()
  })
  
  output$downloadScHeatmapFig <- downloadHandler(
    filename = "singleCellByCellType.pdf",
    content = function(file) {
      pdf(file, width=10, height=max(4,length(input$sc_gene_names_cell_type)*0.5), useDingbats=FALSE)
      print(sc_heatmap())
      dev.off()
    }
  )
  
  
  
  sc_umap <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      plot_sc_umap(
        input$sc_umap_type,
        input$scUmapCellTypeGrouped
      )
    })
  })
  
  output$sc_umap <- renderPlot({
    sc_umap()
  })
  
  output$downloadScUmapFig <- downloadHandler(
    filename = "umap.pdf",
    content = function(file) {
      pdf(file, width=9, height=7, useDingbats=FALSE)
      print(sc_umap())
      dev.off()
    }
  )
  
  
  
  sc_track_genes_exp <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      if (length(input$sc_track_gene_names) == 0){
        ggplot()
      } else {
        plot_sc_track_genes_exp(
          input$sc_track_gene_names,
          input$sc_track_type,
          input$scTrackCellTypeGrouped
        )
      }
    })
  })
  
  output$sc_track_genes_exp <- renderPlot({
    sc_track_genes_exp()
  })
  
  output$downloadScTrackGenesExpFig <- downloadHandler(
    filename = "track_genes_exp.pdf",
    content = function(file) {
      pdf(file, width=10, height=max(5,length(input$sc_track_gene_names)*0.5), useDingbats=FALSE)
      print(sc_track_genes_exp())
      dev.off()
    }
  )
  
  
  
  sc_umap_genes_exp <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      if (length(input$sc_umap_gene_names) == 0){
        ggplot()
      } else {
        plot_sc_umap_genes_exp(
          input$sc_umap_gene_names
        )
      }
    })
  })
  
  output$sc_umap_genes_exp <- renderPlot({
    sc_umap_genes_exp()
  })
  
  output$downloadScUmapGenesExpFig <- downloadHandler(
    filename = "umap_genes_exp.pdf",
    content = function(file) {
      pdf(file, width=7, height=max(4,length(input$sc_umap_gene_names)*0.7), useDingbats=FALSE)
      print(sc_umap_genes_exp())
      dev.off()
    }
  )
  
  
  
  sc_umap_gene_exp <- reactive({
    withProgress(message = 'Plotting', value = 0.80, {
      if (length(input$sc_umap_gene_selected) == 0){
        ggplot()
      } else {
        plot_sc_umap_gene_exp(
          input$sc_umap_gene_selected
        )
      }
    })
  })
  
  output$sc_umap_gene_exp <- renderPlot({
    sc_umap_gene_exp()
  })
  
  output$downloadScUmapGeneExpFig <- downloadHandler(
    filename = "umap_gene_exp.pdf",
    content = function(file) {
      pdf(file, width=8, height=7, useDingbats=FALSE)
      print(sc_umap_gene_exp())
      dev.off()
    }
  )
  
  
  
  
  

  get_overview_table <- function() {
    df_detection = get_detection(df_ms, df_design,
                           input$overview_select_condition,
                           input$overview_select_subgroup,
                           input$overview_max_sd,
                           input$overview_min_pep_uniq,
                           input$overview_min_spat
    )
    df_selected = df_detection[which(df_detection$mean_log2_int>=input$overview_detection_brush$ymin),]
    df_selected = df_selected[which(df_selected$mean_log2_int<=input$overview_detection_brush$ymax),]
    df_selected = df_selected[which(df_selected$detection>=input$overview_detection_brush$xmin),]
    df_selected = df_selected[which(df_selected$detection<=input$overview_detection_brush$xmax),]

    del_col = c(df_design$bclq_id, "ens_gene", "ens_trans")
    del_col = unlist(lapply(del_col, function(x) grep(x, colnames(df_selected))))
    df_selected = df_selected[,-del_col]

    gene_col = "Gene.names"
    gene_col = which(colnames(df_selected) == gene_col)
    uniprot_col = "uniprot"
    uniprot_col = which(colnames(df_selected) == uniprot_col)
    tmp = df_selected[, -c(gene_col, uniprot_col)]
    df_selected = cbind(df_selected[, c(gene_col, uniprot_col)], tmp)

    return(df_selected)
  }

  output$overview_detection_table <- renderDataTable({
    if (is.null(input$overview_detection_brush$xmin)) return(data.frame())
    else {
      return(get_overview_table())
    }
  })

  output$overview_detection_table_dl <- downloadHandler(
    file = "percenage_presence_table.tsv",
    content = function(file) {
      df_selected = get_overview_table()
      write.table(df_selected, file, row.names = FALSE, sep="\t", quote=FALSE)
    }
  )

  observeEvent(input$analyze, {
    withProgress(message = 'Computing...', value = 0.80, {
      file_log = file.path(path_project, "use.log")
      write(format(Sys.time(), "%a %b %d %X %Y"),file=file_log,append=TRUE)

      df_dep_res = run_DEP(
        df_ms, df_design,
        input$dep_select_condition, input$dep_select_comparison,
        input$dep_imputation, input$dep_min_peptide,
        df_anno, alpha=0.01, lfc=1, input$dep_num_missval
      )

      if(df_dep_res$imputed_error) {
        showModal(modalDialog(
          title = "Warning",
          "There are too many missing values to impute the dataset with knn.
           Missing values were replaced by 0. Please reduce the number of
           missing values per condition or use the zero or knn imputation
           method.",
          easyClose = TRUE,
          footer = NULL
        ))
      }


      if( !is.null(df_dep_res) ) {
        list_genes_dep = df_dep_res$dep %>% rowData() %>% data.frame() %>% select(ends_with("name"))
        tab_list_genes_dep = tab_uni2ens_uniq_uniprot[which(tab_uni2ens_uniq_uniprot$uniprot %in% list_genes_dep$name),]
        tmp = list_genes_dep[which(!list_genes_dep$name %in% tab_uni2ens_uniq_uniprot$uniprot),]
        dep_list_genes = split(as.matrix(tab_list_genes_dep[,c("Gene.names", "uniprot", "ens_trans", "ens_gene")]), tab_list_genes_dep$uniprot)

        updateSelectizeInput(session, "dep_select_gene", choices=dep_list_genes, selected=dep_list_genes[[1]][2], server = TRUE)

        list_contrast = colnames(df_dep_res$dep %>% rowData())
        list_contrast = list_contrast[grep("_diff", list_contrast)]
        list_contrast = gsub("_diff", "", list_contrast)
        updateSelectInput(session, "dep_select_contrast_volcano", choices=list_contrast, selected=list_contrast[1])


      } else {
        updateSelectizeInput(session, "dep_select_gene", choices=c(), selected="", server = TRUE)
        updateSelectInput(session, "dep_select_contrast_volcano", choices=c(), selected="")
      }

    })

    ####################################################################
    # reactive
    ####################################################################
    dep_signif_heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        df = df_dep_res$dep %>% .[rowData(.)$significant, ] %>% rowData() %>% data.frame()

        if (length(df[,1]) == 0 ) {
          return(ggplot() + ggtitle("No significant genes."))
        }

        res = get_dep_table(df_dep_res, df_anno, "imputed")
        res = res[which(res$SPAT_score>= input$heatmap_dep_spat),]
        res = res[which(res$uniprot %in% df$uniprot),]

        col_ids = which(colnames(res) %in% df_dep_res$design$label)
        df = res[,col_ids]

        annotation_col = data.frame(
                            condition = df_dep_res$design$condition,
                            batch = as.character(df_dep_res$design$batch)
                        )
        rownames(annotation_col) = df_dep_res$design$label

        batch_colors = colorRampPalette(brewer.pal(n=9, name="Greens"))(length(unique(df_dep_res$design$batch)))
        names(batch_colors) = unique(df_dep_res$design$batch)
        condition_colors = colorRampPalette(brewer.pal(n=8, name="Set1"))(length(unique(df_dep_res$design$condition)))
        names(condition_colors) = unique(df_dep_res$design$condition)
        annotation_colors = list(
          batch=batch_colors,
          condition=condition_colors
        )

        pheatmap(
          mat=df,
          color = inferno(10),
          cutree_cols = input$heatmap_dep_num_col_cluster,
          cutree_rows = input$heatmap_dep_num_row_cluster,
          cluster_rows=TRUE,
          cluster_cols=TRUE,
          show_rownames=FALSE,
          fontsize_col = 12,
          annotation_col = annotation_col,
          annotation_colors = annotation_colors
        )
      })
    })

    dep_var_umap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        df = df_dep_res$dep %>% .[rowData(.)$significant, ] %>% rowData() %>% data.frame()

        if (length(df[,1]) == 0 ) {
          return(ggplot() + ggtitle("No significant genes."))
        }

        res = get_dep_table(df_dep_res, df_anno, "imputed")
        res = res[which(res$SPAT_score>= input$umap_dep_spat),]

        col_ids = which(colnames(res) %in% df_dep_res$design$label)
        int = res[,col_ids]

        RowVar <- function(x, ...) {
          rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
        }
        var_value = RowVar(int)
        int["var"] = var_value
        int = int[order(int$var, decreasing = TRUE), ]
        max_prot = length(df[,1])
        top_x = min(max_prot, input$num_gene_umap)
        most_var = int[1:top_x, 1:length(col_ids)]

        df_umap = umap(t(most_var), config=custom.settings)
        df_res = as.data.frame(df_umap$layout)
        colnames(df_res) = c("x", "y")
        df_res$condition = df_dep_res$design$condition
        df_res$batch = as.character(df_dep_res$design$batch)
        df_res$sample = as.character(df_dep_res$design$label)

        gg_umap = ggplot(df_res, aes_string(x="x", y="y", shape=input$sel_umap_col, color=input$sel_umap_col))
        gg_umap = gg_umap + geom_point(size = 5)
        gg_umap = gg_umap + scale_shape_manual(values = rep(c(0:2,5,6,15:18), 5))
        gg_umap = gg_umap + theme_bw()
        gg_umap = gg_umap + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5))
        gg_umap = gg_umap + ggtitle(paste0("UMAP build on ", top_x, " most var proteins"))

        return(gg_umap)
      })
    })

    dep_var_heatmap_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        res = get_dep_table(df_dep_res, df_anno, "imputed")

        res = res[which(res$SPAT_score>= input$heatmap_dep_spat),]

        col_ids = which(colnames(res) %in% df_dep_res$design$label)
        df = res[,col_ids]

        RowVar <- function(x, ...) {
          rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
        }


        r_var = RowVar(df)
        df$var = r_var
        df = df[order(df$var, decreasing = TRUE), ]
        top_x = min(length(r_var),input$num_gene_heatmap )
        df = df[1:top_x, 1:length(col_ids)]

        annotation_col = data.frame(
                            condition = df_dep_res$design$condition,
                            batch = df_dep_res$design$batch
                        )
        rownames(annotation_col) = df_dep_res$design$label

        batch_colors = colorRampPalette(brewer.pal(n=9, name="Greens"))(length(unique(df_dep_res$design$batch)))
        names(batch_colors) = unique(df_dep_res$design$batch)
        condition_colors = colorRampPalette(brewer.pal(n=8, name="Set1"))(length(unique(df_dep_res$design$condition)))
        names(condition_colors) = unique(df_dep_res$design$condition)
        annotation_colors = list(
          batch=batch_colors,
          condition=condition_colors
        )

        pheatmap(
          mat=df,
          color = inferno(10),
          cutree_cols = input$heatmap_dep_num_col_cluster,
          cutree_rows = input$heatmap_dep_num_row_cluster,
          cluster_rows=TRUE,
          cluster_cols=TRUE,
          show_rownames=FALSE,
          fontsize_col = 12,
          annotation_col = annotation_col,
          annotation_colors = annotation_colors
        )
      })
    })

    volcano_input <- reactive({
      if (input$dep_select_contrast_volcano != "") {
        withProgress(message = 'Plotting', value = 0.80, {
          dep_volcano(
            df_dep_res$dep, df_anno,
            input$dep_select_contrast_volcano,
            input$num_labeled_volcano, input$pep_dep_volcano,
            input$spat_dep_volcano, input$genes_volcano,
            input$volcano_select_color,
            input$volcano_select_alpha
          )
        })
      }
    })

    numbers_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        plot_numbers(df_dep_res$filt)
      })
    })

    coverage_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        plot_coverage(df_dep_res$filt)
      })
    })

    norm_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        plot_normalization(df_dep_res$filt, df_dep_res$norm)
      })
    })

    detect_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        plot_detect(df_dep_res$filt)
      })
    })

    missval_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        plot_missval(df_dep_res$filt)
      })
    })

    imputation_input <- reactive({
      withProgress(message = 'Plotting', value = 0.80, {
        plot_imputation(df_dep_res$norm, df_dep_res$imputed)
      })
    })

    dep_centred_intensity_input <- reactive({
      if (input$dep_select_gene != "") {
        withProgress(message = 'Plotting', value = 0.80, {
          gene_name = input$dep_select_gene
          row_gene_id = which(tab_list_genes_dep$Gene.names == gene_name | tab_list_genes_dep$uniprot == gene_name | tab_list_genes_dep$ens_gene == gene_name | tab_list_genes_dep$ens_trans == gene_name)
          gene_name = tab_list_genes_dep$uniprot[row_gene_id]
          plot_single(df_dep_res$dep, tab_uni2ens_uniq_uniprot, gene_name, type = "centered", df_norm=df_dep_res$norm)
        })
      }
    })

    dep_intensity_input <- reactive({
      if (input$dep_select_gene != "") {
        withProgress(message = 'Plotting', value = 0.80, {
          df_imputed = get_dep_table(df_dep_res, df_anno, 'imputed')
          df_imputed = melt_df_imputed(df_imputed, df_dep_res, input$dep_select_gene, tab_list_genes_dep)

          return(plot_prot_intensity(df_imputed))
        })
      }
    })
    

    ####################################################################
    # Output functions
    ####################################################################
    
    output$dep_var_umap <- renderPlot({
      dep_var_umap_input()
    })

    output$dep_var_heatmap <- renderPlot({
      dep_var_heatmap_input()
    })

    output$volcano <- renderPlot({
      volcano_input()
    })

    output$norm <- renderPlot({
      norm_input()
    })

    output$missval <- renderPlot({
      missval_input()
    })

    output$detect <- renderPlot({
      detect_input()
    })

    output$imputation <- renderPlot({
      imputation_input()
    })

    output$numbers <- renderPlot({
      numbers_input()
    })

    output$coverage <- renderPlot({
      coverage_input()
    })

    output$dep_centred_intensity <- renderPlot({
      dep_centred_intensity_input()
    })

    output$dep_intensity <- renderPlot({
      dep_intensity_input()
    })
    
    build_table_dep_volcano <- function() {
      df = dep_volcano(
        df_dep_res$dep, df_anno, input$dep_select_contrast_volcano,
        input$num_labeled_volcano, input$pep_dep_volcano,
        input$spat_dep_volcano, input$genes_volcano,
        input$volcano_select_color,
        input$volcano_select_alpha,
        plot=FALSE
      )
      df = df[which(df$log2FC >= input$dep_volcano_brush$xmin & df$log2FC <= input$dep_volcano_brush$xmax),]
      df = df[which(df$pvadj >= input$dep_volcano_brush$ymin & df$pvadj <= input$dep_volcano_brush$ymax),]

      sel_col_name = c("Gene.names", "uniprot", "log2FC", "pvalues", "pvadj", "unique", "SPAT_score", "gene_biotype")
      df = df[,which(colnames(df) %in% sel_col_name)]

      colnames(df)[3]="-log10(pv.adj)"
      colnames(df)[4]="-log10(pvalues)"
      colnames(df)[5]="# unique peptide"

      df = df[,c(1,6,7,2,3,4,5,8)]
      return(df)
    }
    ####################################################################
    # click functions
    ####################################################################
    output$table_dep_volcano <- renderDataTable({
      if (is.null(input$dep_volcano_brush$xmin)) return(NULL)
      else {
        return(build_table_dep_volcano())
      }
    })


    output$table_dep_int_centred <- renderTable({
      if (is.null(input$dep_centred_int_brush$xmin)) return(data.frame())
      else {
        gene_name = input$dep_select_gene
        row_gene_id = which(tab_list_genes_dep$Gene.names == gene_name | tab_list_genes_dep$uniprot == gene_name | tab_list_genes_dep$ens_gene == gene_name | tab_list_genes_dep$ens_trans == gene_name)
        gene_name = tab_list_genes_dep$uniprot[row_gene_id]
        res = plot_single(df_dep_res$dep, tab_uni2ens_uniq_uniprot, gene_name, type = "centered", plot=FALSE, df_norm=df_dep_res$norm)

        id_cond = c(1:length(levels(res$condition)))
        id_cond = which(id_cond >= input$dep_centred_int_brush$xmin & id_cond <= input$dep_centred_int_brush$xmax)
        x_cond = levels(res$condition)[id_cond]
        res = res[which(res$condition %in% x_cond),]
        res = res[which(res$log2_Centred_intensity>=input$dep_centred_int_brush$ymin),]
        res = res[which(res$log2_Centred_intensity<=input$dep_centred_int_brush$ymax),]

        return(res)
      }
    })


    output$table_dep_int <- renderTable({
      if (is.null(input$dep_int_brush$xmin)) return(data.frame())
      else {
        df_imputed = get_dep_table(df_dep_res, df_anno, 'imputed')
        df_imputed = melt_df_imputed(df_imputed, df_dep_res, input$dep_select_gene, tab_list_genes_dep)
        id_cond = c(1:length(unique(df_imputed$condition)))
        id_cond = which(id_cond >= input$dep_int_brush$xmin & id_cond <= input$dep_int_brush$xmax)
        x_cond = levels(as.factor(df_imputed$condition))[id_cond]
        df_imputed = df_imputed[which(df_imputed$condition %in% x_cond),]
        df_imputed = df_imputed[which(df_imputed['log2(Normalized Imputed Intensity + 1)']>=input$dep_int_brush$ymin),]
        df_imputed = df_imputed[which(df_imputed['log2(Normalized Imputed Intensity + 1)']<=input$dep_int_brush$ymax),]
        print(colnames(df_imputed))
        df_imputed = df_imputed[, c(16, 1, 12, 13)]

        return(df_imputed)
      }
    })

    ####################################################################
    # download figure functions
    ####################################################################
    output$downloadDepVarUmap <- downloadHandler(
      filename = 'umap.pdf',
      content = function(file) {
        pdf(file, height = 10, width=14, paper = "a4r", useDingbats=FALSE)
        print(dep_var_umap_input())
        dev.off()
      }
    )

    output$downloadDepVarHeatmap <- downloadHandler(
      filename = 'Heatmap.pdf',
      content = function(file) {
        pdf(file, height = 10, paper = "a4", useDingbats=FALSE)
        print(dep_var_heatmap_input())
        dev.off()
      }
    )

    output$downloadVolcano <- downloadHandler(
      filename = function() {
        paste0("Volcano.pdf")
      },
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(volcano_input())
        dev.off()
      }
    )

    output$downloadNorm <- downloadHandler(
      filename = "normalization.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(norm_input())
        dev.off()
      }
    )

    output$downloadMissval <- downloadHandler(
      filename = "missing_values_heatmap.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(missval_input())
        dev.off()
      }
    )

    output$downloadDetect <- downloadHandler(
      filename = "missing_values_quant.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        gridExtra::grid.arrange(detect_input())
        dev.off()
      }
    )

    output$downloadImputation <- downloadHandler(
      filename = "imputation.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(imputation_input())
        dev.off()
      }
    )

    output$downloadNumbers <- downloadHandler(
      filename = "numbers.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE, width=14, height=7)
        print(numbers_input())
        dev.off()
      }
    )

    output$downloadCoverage <- downloadHandler(
      filename = "coverage.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(coverage_input())
        dev.off()
      }
    )

    output$downloadDepCentredIntensity <- downloadHandler(
      filename = "dep_centred_intensity.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(dep_centred_intensity_input())
        dev.off()
      }
    )

    output$downloadDepIntensity <- downloadHandler(
      filename = "dep_centred_intensity.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(dep_intensity_input())
        dev.off()
      }
    )

    ####################################################################
    # download table functions
    ####################################################################
    output$table_dep_volcano_dw <- downloadHandler(
      file = "volcano selected_proteins.tsv",
      content = function(file) {
        df = data.frame()
        if (!is.null(input$dep_volcano_brush$xmin)) df = build_table_dep_volcano()
        write.table(df, file, row.names = FALSE, sep="\t", quote=FALSE)
      }
    )

  })


  ########################### Transcriptome differential analysis ##########################

  observeEvent(input$trans_analyze, {
    withProgress(message = 'Computing...', value = 0.80, {

      file_log = file.path(path_diff, "use.log")
      write(format(Sys.time(), "%a %b %d %X %Y"),file=file_log,append=TRUE)

      trans_comparaison = input$transcripome_selected_comparison
      trans_comparaison = strsplit(trans_comparaison, ' ')[[1]][1]

      file_exp = file.path(path_diff, "STAR_RSEM_unstranded_v32", 'leucegene2_430', "readcounts.xls")
      all_gene_exp = fread(file_exp, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
      setDF(all_gene_exp)
      df_cpm_anno <<- update_all_gene_exp(all_gene_exp)

      analysis_dir = file.path(path_diff, "design", 'leucegene2_430', trans_comparaison, "STAR_RSEM_unstranded_v32")

      id_col = which(colnames(df_cpm_anno) %in% c("Gene.names", "ID2"))
      list_genes = split(as.matrix(df_cpm_anno[,id_col]), df_cpm_anno$ID2)
      old_gene = input$trans_selected_gene
      if (old_gene == "" | ! old_gene %in% list_genes) old_gene = "MYH11"
      if (! old_gene %in% list_genes) old_gene = list_genes[[1]][2]
      updateSelectizeInput(session, "trans_selected_gene", choices=list_genes, selected=old_gene, server = TRUE)

      design_dir <<- file.path(path_diff, "design", 'leucegene2_430', input$transcripome_selected_comparison)
      file_design = file.path(design_dir, "design.tsv")
      tab_design <<- fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
      setDF(tab_design)

      dir_res = file.path(analysis_dir, "counts")
      file_in = file.path(dir_res, "limma_voom_genes.xls")
      limma_res <<- NULL
      if(file.exists(file_in)){
        limma_res <<- read.table(file_in, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
      }

      trans_subgoups <<- unlist(strsplit(trans_comparaison, "_vs_"))

    })


    limma_voom_fig <- reactive({
      withProgress(message = 'Plotting...', value = 0.80, {
        plot_limma_voom(
          limma_res, trans_surface_score,
          input$trans_num_gene_labeled, input$tans_min_spat,
          input$trans_volcano_select_color, input$trans_volcano_select_alpha,
          genes_highligth=input$trans_limma_highlight_genes
        )
      })
    })

    limma_voom_tab <- reactive({
      withProgress(message = 'Plotting...', value = 0.80, {
        plot_limma_voom(
          limma_res, trans_surface_score,
          input$trans_num_gene_labeled, input$tans_min_spat,
          input$trans_volcano_select_color, input$trans_volcano_select_alpha,
          table=TRUE
        )
      })
    })

    build_table_limma_volcano <- function () {
      df = plot_limma_voom(
        limma_res, trans_surface_score,
        input$trans_num_gene_labeled, input$tans_min_spat,
        input$trans_volcano_select_color, input$trans_volcano_select_alpha,
        table=TRUE
      )
      df = df[which(df$logFC >= input$limma_volcano_brush$xmin & df$logFC <= input$limma_volcano_brush$xmax),]
      df = df[which(df$log10_adj.P.Val >= input$limma_volcano_brush$ymin & df$log10_adj.P.Val <= input$limma_volcano_brush$ymax),]

      sel_col_name = c("ID", "Gene.names", "logFC", "log10_adj.P.Val", "log10_P.Value", "SPAT_score", "gene_biotype")
      df = df[,which(colnames(df) %in% sel_col_name)]

      colnames(df)[1]="Ensembl"
      colnames(df)[2]="log2FC"
      colnames(df)[3]="-log10(pv.adj)"
      colnames(df)[4]="-log10(pvalues)"

      df = df[,c(1,6,7,2,3,4,5)]
      return(df)
    }

    output$table_limma_volcano <- renderDataTable({
      if (is.null(input$limma_volcano_brush$xmin)) return(NULL)
      else {
        return(build_table_limma_volcano())
      }
    })


    output$trans_limma_voom <- renderPlot({
      limma_voom_fig()
    })

    output$downloadTransLimmaVoom <- downloadHandler(
      filename = 'trans_diff_limma_voom.pdf',
      content = function(file) {
        pdf(file, height = 10, paper = "a4", useDingbats=FALSE)
        print(limma_voom_fig())
        dev.off()
      }
    )

    output$table_limma_volcano_dw <- downloadHandler(
      filename = 'limma_voom_surface.tsv',
      content = function(file) {
        df = data.frame()
        if (!is.null(input$limma_volcano_brush$xmin)) df = build_table_limma_volcano()
        write.table(df, file, row.names = FALSE, sep="\t", quote=FALSE)
      }
    )


    trans_gene_fig <- reactive({
      withProgress(message = 'Plotting...', value = 0.80, {
        plot_trans_gene(trans_subgoups, df_cpm_anno, tab_design, input$trans_selected_gene, input$trans_selected_graphic, 'readcounts.xls')
      })
    })

    output$trans_gene <- renderPlot({
      trans_gene_fig()
    })

    output$downloadTransGene <- downloadHandler(
      filename = "gene_readcounts.pdf",
      content = function(file) {
        pdf(file, useDingbats=FALSE)
        print(trans_gene_fig())
        dev.off()
      }
    )

    output$trans_selected_sample_table <- renderTable({
      if (is.null(input$trans_gene_brush$xmin)) return(data.frame())
      else {
        exp_data = get_gene_exp(trans_subgoups, df_cpm_anno, tab_design, input$trans_selected_gene, input$trans_selected_graphic, 'readcounts.xls')

        id_group = c(1:length(levels(exp_data$group)))
        id_group = which(id_group >= input$trans_gene_brush$xmin & id_group <= input$trans_gene_brush$xmax)
        x_group = levels(exp_data$group)[id_group]

        exp_data = exp_data[which(exp_data$group %in% x_group),]

        exp_data = exp_data[which(exp_data$log2_CPM>=input$trans_gene_brush$ymin),]
        exp_data = exp_data[which(exp_data$log2_CPM<=input$trans_gene_brush$ymax),]
        selected_col = which(colnames(exp_data) %in% c("sample", "log2_CPM"))

        return(exp_data[,selected_col])
      }
    })
  }) ##observeEvent
})
