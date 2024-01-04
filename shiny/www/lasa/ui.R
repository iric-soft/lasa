side_width = 300

shinyUI(
  navbarPage(
    "Leucegene: LASA",
    header = tagList(
      useShinydashboard()
    ),
    tabPanel("Home",
      fluidRow(
        column(4),
        column(4,
          imageOutput("logo", height = 200),
        ),
        column(4)
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            HTML('<h4> This web portal was developed as part of the
                  <a href="https://leucegene.ca/">Leucegene</a>
                  project and reports: </h4>')
          )
        ),
        column(1),
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            HTML(
              '<ul>
                <li>Mass spectrometry data from proteins identified in the surface proteome of 100 primary human AML specimens.</li>
                <li>RNA sequencing data from 430 primary human AML specimens.</li>
                <li>Single-cell RNA sequencing data from 20 primary human AML specimens.</li>
              </ul>'
            )
          )
        ),
        column(1)
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            h4("More details on data and analyses are available in the publication: Immunotherapeutic targeting of surfaceome heterogeneity in AML [ref].")
          )
        ),
        column(1)
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            h4("Please use the bar on top of this web page to navigate between sections:")
          )
        ),
        column(1)
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            HTML(
              '<ul>
                <li>Surfaceome - Dataset overview: To explore surface proteins according to their percentage of detection in samples of the Leucegene AML surfaceome cohort.</li>
                <li>Surfaceome - Differential analysis: To perform differential analyses, using DEP.</li>
                <li>Transcriptome - Differential analysis: To perform differential analyses, using limma voom.</li>
                <li>Single Cell RNA sequencing - Data analysis </li>
                <li>Download: To download quantitative datasets and design files used in this portal.</li>
              </ul>'
            )
          )
        ),
        column(1)
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            h4("References:")
          )
        ),
        column(1)
      ),
      fluidRow(
        column(1),
        column(10,
          div(
            HTML(
              '<ul>
                <li><a href="https://www.maxquant.org/">Maxquant</a>: Cox J and Mann M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. Nat Biotechnol, 2008 p.1367-72. </li>
                <li><a href="https://bioconductor.org/packages/release/bioc/html/DEP.html">DEP</a>: Zhang X, Smits A, van Tilburg G, Ovaa H, Huber W and Vermeulen M. Proteome-wide identification of ubiquitin interactions using UbIA-MS. Nature Protocols, 2018, p.530–550. </li>
                <li><a href="https://bioconductor.org/packages/release/bioc/html/limma.html">Limma voom</a>: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 2015 e47. </li>
                <li><a href="https://spat.leucegene.ca/">SPAT</a>: Spinella, J.-F., Theret, L., Aubert, L., Audemard, E., Boucher, G., Pfammatter, S., Bonneil, É., Bordeleau, M.-E., Thibault, P., Hébert, J., Roux, P.P., Sauvageau G. (2023). SPAT: Surface Protein Annotation Tool. bioRxiv, doi: 10.1101/2023.07.07.547075. </li>
              </ul>'
            )
          )
        ),
        column(1)
      )
    ),
    tabPanel("Surfaceome - Dataset overview",
      fluidRow(
        column(
          width = 12,
          box(
            title = "Explore proteins according to their percentage of detection",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              column(
                width = 2,
                selectInput(
                  "overview_select_condition", "Select a condition:",
                  choices=list_condition, selected="subgroup"
                )
              ),
              column(
                width = 2,
                selectInput(
                  "overview_select_subgroup", "Select an annotation:",
                  choices=list_subgroup, selected="All"
                )
              ),
              column(
                width = 2,
                numericInput(
                  "overview_max_sd",
                  "Max standard deviation (on protein intensity):",
                  min=0, max=100, value=10
                )
              ),
              column(
                width = 2,
                numericInput(
                  "overview_min_pep_uniq",
                  "Min # of unique peptide:",
                  min=0, max=1000, value=2
                )
              ),
              column(
                width = 2,
                numericInput(
                  "overview_min_spat",
                  "Min SPAT:",
                  min=0, max=10, value=8
                )
              )
            ),
            fluidRow(
              column(
                width = 2,
                selectInput(
                  "overview_select_color", "Select legend color:",
                  choices=colors(), selected="red"
                )
              ),
              column(
                width = 2,
                numericInput(
                  "overview_select_alpha",
                  "Dot opacity:",
                  min=0, max=1, value=0.5
                )
              ),
              column(
                style = "margin-top: 25px;",
                width = 2,
                downloadButton('downloadOverviewFig', 'Figure')
              )
            ),
            fluidRow(
              column(
                width = 12,
                plotOutput("overview_detection",
                  height = plot_height[4],
                  brush="overview_detection_brush"
                )
              )
            ),
            fluidRow(
              box(
                title = "Table of selected proteins",
                status = "primary",
                solidHeader = FALSE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE,
                dataTableOutput('overview_detection_table'),
                downloadButton('overview_detection_table_dl',"Download table")
              )
            )
          )
        )
      )
    ),
    tabPanel("Surfaceome - Differential analysis",
      fluidRow(
        column(
          width = 12,
          box(
            title = "DEP differential analysis parameters",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "dep_select_condition", "Select a comparison design:",
                  choices=list_condition, selected="subgroup"
                )
              ),
              column(
                width = 4,
                selectInput(
                  "dep_select_comparison", "Select a comparison method:",
                  choices=list_comparison, selected="Compare all condition pairs"
                )
              )
            ),
            fluidRow(
              column(
                width = 3,
                numericInput(
                  "dep_min_peptide", "Min. # unique peptide:",
                  min=0, max=100, value=2
                )
              ),
              column(
                width = 3,
                numericInput(
                  "dep_num_missval", "Max. # of missing value by condition:",
                  min=0, max=100, value=0
                ),
                p(a("Need help ?",
                    href = "https://rdrr.io/bioc/DEP/man/filter_missval.html",
                    target="_blank")
                )
              ),
            ),
            fluidRow(
              column(
                width = 4,
                radioButtons(
                  "dep_imputation", "Imputation type",
                  choices = c("knn", 'min', 'zero'),
                  selected = "zero"
                ),
                p(a("Need help to select imputation ?",
                    href = "https://www.rdocumentation.org/packages/MSnbase/versions/1.20.7/topics/impute-methods",
                    target="_blank")
                )
              )
            ),
            fluidRow(
              column(
                width = 4,
                actionButton("analyze", "Run analysis")
              )
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          box(
            title = "Analysis results",
            status = "success",
            solidHeader = TRUE,
            width = 8,
            fluidRow(
              tabBox(
                title = " ", width = 12,
                tabPanel(
                  title = "DEP QC", width = 12,
                  fluidRow(
                    tabBox(
                      title = "Plots", width = 12,
                      tabPanel(
                        title = "Protein #",
                        plotOutput("numbers", height = plot_height[2]),
                        downloadButton('downloadNumbers', 'Figure')
                      ),
                      tabPanel(
                        title = "Sample cov.",
                        plotOutput("coverage", height = plot_height[2]),
                        downloadButton('downloadCoverage', 'Figure')
                      ),
                      tabPanel(
                        title = "Normalization",
                        plotOutput("norm", height = plot_height[5]),
                        downloadButton('downloadNorm', 'Figure')
                      ),
                      tabPanel(
                        title = "Missing values",
                        plotOutput("detect", height = plot_height[2]),
                        downloadButton('downloadDetect', 'Figure')
                      ),
                      tabPanel(
                        title = "Missing value heatmap",
                        plotOutput("missval", height = plot_height[2]),
                        downloadButton('downloadMissval', 'Figure')
                      ),
                      tabPanel(
                        title = "Imputation",
                        plotOutput("imputation", height = plot_height[2]),
                        downloadButton('downloadImputation', 'Figure')
                      )
                    )
                  )
                ),
                tabPanel(
                  title = "Other QC", width = 12,
                  fluidRow(
                    tabBox(
                      title = "Plots", width = 12,
                      tabPanel(
                        title = "Heatmap",
                        fluidRow(
                          column(
                            width = 2,
                            numericInput(
                              "heatmap_dep_spat",
                              "Min SPAT:",
                              min=0, max=10, value=0
                            )
                          ),
                          column(
                            width = 2,
                            numericInput(
                              "heatmap_dep_num_row_cluster",
                              "Num row cluster:",
                              min=1, max=20, value=1
                            )
                          ),
                          column(
                            width = 2,
                            numericInput(
                              "heatmap_dep_num_col_cluster",
                              "Num. col. cluster:",
                              min=1, max=20, value=6
                            )
                          ),
                          column(
                            width = 3,
                            numericInput(
                              "num_gene_heatmap",
                              "Top x most variant genes:",
                              min=1, max=10000, value=100
                            )
                          ),
                          column(
                            style = "margin-top: 25px;",
                            width = 2,
                            downloadButton('downloadDepVarHeatmap', 'Figure')
                          )
                        ),
                        fluidRow(
                          column(
                            width = 12,
                            plotOutput("dep_var_heatmap", height = plot_height[3])
                          )
                        )
                      ),
                      tabPanel(
                        title = "UMAP",
                        fluidRow(
                          column(
                            width = 2,
                            numericInput(
                              "umap_dep_spat",
                              "Min SPAT:",
                              min=0, max=10, value=0
                            )
                          ),
                          column(
                            width = 2,
                            selectInput(
                              "sel_umap_col",
                              "Colored by:",
                              choices=c("condition", "batch"),
                              selected="condition"
                            )
                          ),
                          column(
                            width = 3,
                            numericInput(
                              "num_gene_umap",
                              "Top x most variant genes:",
                              min=1, max=10000, value=100
                            )
                          ),
                          column(
                            style = "margin-top: 25px;",
                            width = 2,
                            downloadButton('downloadDepVarUmap', 'Figure')
                          )
                        ),
                        fluidRow(
                          column(
                            width = 12,
                            plotOutput("dep_var_umap", height = plot_height[3])
                          )
                        )
                      )
                    )
                  )
                ),
                tabPanel(
                  title = "Volcano",
                  fluidRow(
                    column(
                      width = 3,
                      selectInput(
                        "dep_select_contrast_volcano",
                        "Choose a comparison:",
                        choices=c()
                      )
                    ),
                    column(
                      width = 3,
                      numericInput(
                        "num_labeled_volcano",
                        "Number of labeled proteins:",
                        min=1, max=50, value=10
                      )
                    ),
                    column(
                      width = 2,
                      numericInput(
                        "pep_dep_volcano",
                        "Min # of unique peptide:",
                        min=0, max=1000, value=2
                      )
                    ),
                    column(
                      width = 2,
                      numericInput(
                        "spat_dep_volcano",
                        "Min SPAT:",
                        min=0, max=10, value=0
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 2,
                      selectInput(
                        "volcano_select_color", "Select legend color:",
                        choices=colors(), selected="red"
                      )
                    ),
                    column(
                      width = 2,
                      numericInput(
                        "volcano_select_alpha",
                        "Dot opacity:",
                        min=0, max=1, value=0.5
                      )
                    ),
                    column(
                      style = "margin-top: 25px;",
                      width = 2,
                      downloadButton('downloadVolcano', 'Figure')
                    )
                  ),
                  fluidRow(
                    plotOutput(
                      "volcano",
                      height = plot_height[2],
                      brush="dep_volcano_brush"
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      textInput("genes_volcano", label ="Enter genes's name to highlight separated by ',' (without space)", value = "")
                    )
                  ),
                  fluidRow(
                    box(
                      title = "Selected proteins table",
                      status = "primary",
                      solidHeader = FALSE,
                      width = 12,
                      collapsible = TRUE,
                      collapsed = FALSE,
                      dataTableOutput('table_dep_volcano'),
                      downloadButton('table_dep_volcano_dw',"Download table")
                    )
                  )
                )
              )
            )
          ),
          box(
            title = "Protein view",
            status = "info",
            solidHeader = TRUE,
            width = 4,
            fluidRow(
              column(
                width = 12,
                selectizeInput(
                  "dep_select_gene",
                  "Choose a gene/protein:",
                  choices=c(),
                  multiple = FALSE
                )
              )
            ),
            fluidRow(
              tabBox(
                title = "", width = 12,
                tabPanel(
                  title = "Intensity", width = 12,
                  fluidRow(
                    column(
                      width = 12,
                      plotOutput("dep_intensity",
                                 height = plot_height[3],
                                 brush="dep_int_brush"
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      downloadButton('downloadDepIntensity', 'Figure')
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      h4("Selected samples:"),
                      tableOutput('table_dep_int')
                    )
                  )
                ),
                tabPanel(
                  title = "Centred Intensity", width = 12,
                  fluidRow(
                    column(
                      width = 12,
                      plotOutput("dep_centred_intensity",
                        height = plot_height[3],
                        brush="dep_centred_int_brush"
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      downloadButton('downloadDepCentredIntensity', 'Figure')
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      h4("Selected samples:"),
                      tableOutput('table_dep_int_centred')
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
    tabPanel("Transcriptome - Differential analysis",
      fluidRow(
        column(
          width = 12,
          box(
            title = "Limma voom differential analysis parameters",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "transcripome_selected_comparison", "Select a comparison design:",
                  choices=trans_list_comparison, selected="68_CK"
                )
              )
            ),
            fluidRow(
              column(
                width = 4,
                actionButton("trans_analyze", "Run analysis")
              )
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          box(
            title = "Analysis results",
            status = "success",
            solidHeader = TRUE,
            width = 8,
            fluidRow(
              column(
                width = 2,
                numericInput(
                  "tans_min_spat",
                  "Min SPAT:",
                  min=0, max=10, value=0
                )
              ),
              column(
                width = 2,
                numericInput(
                  "trans_num_gene_labeled",
                  "Number of labeled genes:",
                  min=0, max=50, value=5
                )
              )
            ),
            fluidRow(
              column(
                width = 2,
                selectInput(
                  "trans_volcano_select_color", "Select legend color:",
                  choices=colors(), selected="red"
                )
              ),
              column(
                width = 2,
                numericInput(
                  "trans_volcano_select_alpha",
                  "Dot opacity:",
                  min=0, max=1, value=0.5
                )
              ),
              column(
                width = 2,
                style = "margin-top: 25px;",
                downloadButton('downloadTransLimmaVoom', 'Figure')
              )
            ),
            fluidRow(
              column(
                width = 12,
                plotOutput(
                  "trans_limma_voom",
                  height = plot_height[2],
                  brush="limma_volcano_brush"
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                textInput("trans_limma_highlight_genes", label ="Enter genes's name to highlight separated by ',' (without space)", value = "")
              )
            ),
            fluidRow(
              box(
                title = "Selected genes table",
                status = "primary",
                solidHeader = FALSE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE,
                dataTableOutput('table_limma_volcano'),
                downloadButton('table_limma_volcano_dw',"Download table")
              )
            )
          ),
          box(
            title = "Gene view",
            status = "info",
            solidHeader = TRUE,
            width = 4,

            fluidRow(
              column(
                width = 6,
                selectizeInput(
                  "trans_selected_gene",
                  "Choose a gene/protein:",
                  choices=c(),
                  multiple = FALSE
                )
              ),
              column(
                width = 6,
                selectInput(
                  "trans_selected_graphic",
                  "Choose a graphic:",
                  choices = c("scatterplot", "boxplot"),
                  selected = "scatterplot"
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                plotOutput("trans_gene",
                           height = plot_height[3],
                           brush="trans_gene_brush"
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                downloadButton('downloadTransGene', 'Figure')
              )
            ),
            fluidRow(
              column(
                width = 12,
                h4("Selected samples:"),
                tableOutput('trans_selected_sample_table')
              )
            )
          )
        )
      )
    ),
    tabPanel("Single-cell RNA sequencing - Data analysis",
      fluidRow(
        column(
          width = 12,
          box(
            title = "Analysis results",
            status = "success",
            solidHeader = TRUE,
            width = 7,
            fluidRow(
              tabBox(
                title = "", width = 12,
                tabPanel(
                  title = "Overview",
                  width = 12,
                  fluidRow(
                    tabBox(
                      title = "", width = 12,
                      tabPanel(
                        title = "By AML subgroup",
                        width = 12,
                        fluidRow(
                          column(
                            width = 8,
                            selectInput("sc_gene_names","Select gene names:",
                                        choices=NULL, multiple = TRUE
                            )
                          ),
#                          column(
#                            width = 4,
#                            style = "margin-top: 25px;",
#                            switchInput(
#                              "scCellType",
#                              label="Cells type grouped", value = TRUE,
#                              inline=TRUE,
#                              labelWidth='100%'
#                            )
#                          ),
                          column(
                            width = 2,
                            selectInput(
                              "scCellTypeSelected", "Select a cell type:",
                              choices=c(
                                "All", "CD34+ HSC", "CD34+ MultiLin", "CD34+ LMPP", "CD34+ Gran",
                                "Mature myeloid lineage", "CD34+ Eo/B/Mast", "CD34+ MDP",
                                "Dendritic cells", "Monocyte", "Erythrocytic lineage",
                                "Megakaryocytic lineage", "CD34+ CLP", "B cell lineage",
                                "Plasma Cell", "T/NK cell lineage", "Stromal"
                              ),
                              selected="All"
                            )
                          ),
                        ),
                        fluidRow(
                          column(
                            width = 3,
                            style = "margin-top: 25px;",
                            switchInput(
                              "scPlotZscore",
                              label="zscore", value = TRUE
                            )
                          ),
                          column(
                            width = 2,
                            selectInput(
                              "scPaletteColor", "Select palette color:",
                              choices=c(
                                "YlOrRd", "RdBu",
                                "RdGy", "PiYG", "PRGn",
                                "BrBG", "Spectral"
                              ),
                              selected="YlOrRd"
                            )
                          ),
#                          conditionalPanel(
#                            condition = "input.scPlotZscore",
#                            column(
#                              width = 2,
#                              radioButtons(
#                                "scPlotZscoreGene",
#                                label="Zscore computed on:",
#                                choices=c("Gene" = "gene", "AML subg" = "subg"),
#                                selected = "gene"
#                              )
#                            ),
#                            column(
#                              width = 3,
#                              radioButtons(
#                                "scPlotZscoreMean",
#                                label="Zscore computed on:",
#                                choices=c(
#                                  "Mean  expression by AML subg" = "mean",
#                                  "Cells expression" = "raw"
#                                ),
#                                selected = "mean"
#                              )
#                            )
#                          ),
                          column(
                            style = "margin-top: 25px;",
                            width = 2,
                            downloadButton('downloadScPlotFig', 'Figure')
                          )
                        ),
                        fluidRow(
                          column(
                            width = 12,
                            plotOutput(
                              "sc_plot",
                              height = plot_height[2]
                            )
                          )
                        )
                      ),
                      tabPanel(
                        title = "By Cell type",
                        width = 12,
                        fluidRow(
                          column(
                            width = 8,
                            selectInput("sc_gene_names_cell_type","Select gene names:",
                                        choices=NULL, multiple = TRUE
                            )
                          ),
                          column(
                            width = 2,
                            selectInput(
                              "sc_heatmap_subgroup", "Select a subgroup:",
                              choices=c(
                                "All","CK", "Del(5q)", "inv(16)",
                                "KMT2A-r", "Mono_7", "NKt", "NPM1-mut",
                                "RUNX1-mut"
                              ),
                              selected="all"
                            )
                          )
                        ),
                        fluidRow(
                          column(
                            width = 2,
                            style = "margin-top: 25px;",
                            switchInput(
                              "scHeatmapZscore",
                              label="zscore", value = TRUE
                            )
                          ),
                          column(
                            width = 2,
                            selectInput(
                              "scHeatmapPaletteColor", "Select palette color:",
                              choices=c(
                                "YlOrRd", "RdBu",
                                "RdGy", "PiYG", "PRGn",
                                "BrBG", "Spectral"
                              ),
                              selected="YlOrRd"
                            )
                          ),
                          column(
                            width = 4,
                            style = "margin-top: 25px;",
                            switchInput(
                              "scHeatmapCellType",
                              label="Cells type grouped", value = TRUE,
                              inline=TRUE,
                              labelWidth='100%'
                            )
                          ),
                          column(
                            style = "margin-top: 25px;",
                            width = 2,
                            downloadButton('downloadScHeatmapFig', 'Figure')
                          )
                        ),
                        fluidRow(
                          column(
                            width = 12,
                            plotOutput(
                              "sc_heatmap",
                              height = plot_height[2]
                            )
                          )
                        )
                      )
                    )
                  )
                ),
                tabPanel(
                  title = "Genes Expression UMAPs",
                  width = 12,
                  fluidRow(
                    column(
                      width = 7,
                      selectInput("sc_umap_gene_names","Select gene names:",
                                  choices=NULL, multiple = TRUE
                      )
                    ),
                    column(
                      style = "margin-top: 25px;",
                      width = 2,
                      downloadButton('downloadScUmapGenesExpFig', 'Figure')
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      plotOutput(
                        "sc_umap_genes_exp",
                        height = plot_height[3]
                      )
                    )
                  )
                ),
                tabPanel(
                  title = "Track",
                  width = 12,
                  fluidRow(
                    column(
                      width = 8,
                      selectInput("sc_track_gene_names","Select gene names:",
                                  choices=NULL, multiple = TRUE
                      )
                    ),
                  ),
                  fluidRow(
                    column(
                      width = 4,
                      selectInput(
                        "sc_track_type", "Select a color legend:",
                        choices=c("Classification", "sample", "AML_subgroup"),
                        selected="Classification"
                      )
                    ),
                    conditionalPanel(
                      condition = "input.sc_track_type == 'Classification'",
                      column(
                        width = 5,
                        style = "margin-top: 25px;",
                        switchInput(
                          "scTrackCellTypeGrouped",
                          label="Cells type grouped", value = TRUE,
                          inline=TRUE,
                          labelWidth='100%'
                        )
                      ),
                    ),
                    column(
                      style = "margin-top: 25px;",
                      width = 2,
                      downloadButton('downloadScTrackGenesExpFig', 'Figure')
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      plotOutput(
                        "sc_track_genes_exp",
                        height = plot_height[3]
                      )
                    )
                  )
                ),
              )
            )
          ),
          box(
            title = "UMAP view",
            status = "info",
            solidHeader = TRUE,
            width = 5,
            fluidRow(
              tabBox(
                title = "", width = 12,
                tabPanel(
                  title = "Annotations",
                  width = 12,
                  fluidRow(
                    column(
                      width = 3,
                      selectInput(
                        "sc_umap_type", "Select a color legend:",
                        choices=c("Classification", "sample", "AML_subgroup"),
                        selected="Classification"
                      )
                    ),
                    conditionalPanel(
                      condition = "input.sc_umap_type == 'Classification'",
                      column(
                        width = 7,
                        style = "margin-top: 25px;",
                        switchInput(
                          "scUmapCellTypeGrouped",
                          label="Cells type grouped", value = TRUE,
                          inline=TRUE,
                          labelWidth='100%'
                        )
                      ),
                    ),
                    column(
                      style = "margin-top: 25px;",
                      width = 2,
                      downloadButton('downloadScUmapFig', 'Figure')
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      plotOutput(
                        "sc_umap",
                        height = plot_height[3]
                      )
                    )
                  )
                ),
                tabPanel(
                  title = "Gene Expression",
                  width = 12,
                  fluidRow(
                    column(
                      width = 10,
                      selectInput(
                        "sc_umap_gene_selected", "Select a gene:",
                        choices=c(),
                      )
                    ),
                    column(
                      style = "margin-top: 25px;",
                      width = 2,
                      downloadButton('downloadScUmapGeneExpFig', 'Figure')
                    )
                  ),
                  fluidRow(
                    column(
                      width = 12,
                      plotOutput(
                        "sc_umap_gene_exp",
                        height = plot_height[3]
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
    tabPanel("Download",
      fluidRow(
        column(
          width = 12,
          box(
            width = 12,
            title = "Download files used for analyses",
            status = "success",
            solidHeader = TRUE,
            collapsible = FALSE,
            fluidRow(
              column(1,),
              column(6,
                downloadButton('downloadMsRaw', 'Surface proteomics data from Maxquant')
              )
            ),
            fluidRow(
              column(1,),
              column(6,
                     downloadButton('downloadSurfaceomeDesign', 'Surface proteomics - analysis design')
              )
            ),
            fluidRow(
              column(1,),
              column(6,
                downloadButton('downloadReadcounts', 'Transcriptome data (readcounts from STAR+RSEM)')
              )
            ),
            fluidRow(
              column(1,),
              column(6,
                downloadButton('downloadTanscriptomeDesign', 'Transcriptome - analysis design')
              )
            ),
            fluidRow(
              column(1,),
              column(6,
                     downloadButton('downloadSingleCellDesign', 'Single-cell RNA sequencing - analysis design')
              )
            )
          )
        )
      )
    )
  )
)

