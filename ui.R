#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

ui <- navbarPage(
    "Mammary Gland Transcriptomes", 
    tabPanel(
        "Gene Plots",
        fluidPage(
            fluidRow(
                shinyjs::useShinyjs(),
                column(3, 
                       textInput("gene",
                                 "Gene name:", 
                                 value = "Zfp57"),
                       actionButton("plot", "Plot")
                ),
                column(9, 
                       textOutput("plotted_gene"),
                       # downloadButton("download_plot_data", "Download plot data"),
                       plotOutput("zfp57_expr"),
                       plotOutput("imprint_expr"),
                       plotOutput("imprint_isolde")
                )
            )
        )
    ),
    tabPanel(
        "Data Download",
        fluidPage(
            fluidRow(
                column(12, 
                       h3("Thresholds for filtering"),
                       fluidRow(
                           column(3,
                                  h4("ZFP57 data"),
                                  sliderInput("fdr_zfp57", "FDR", 
                                              min = 0, max = 1,
                                              value = 0.05, step = 0.05),
                                  sliderInput("fc_threshold", "Fold change",
                                              min = 1, max = 10,
                                              value = 2, step = 1)
                           ),
                           column(3,
                                  h3(""),
                                  h4("Hybrid data"),
                                  checkboxGroupInput("isolde_status", "Isolde status",
                                                     choices = list(
                                                         `Biased expression` = "ASE",
                                                         `Biallelic` = "BA",
                                                         `Undetermined` = "UN",
                                                         `Filtered out` = "FILT"
                                                     ),
                                                     selected = "ASE"),
                                  sliderInput("ase_bias_threshold", "Allele bias threshold",
                                              min = 0, max = 1, value = 0.7, step = 0.1)),
                           column(2,
                                  h4("Cell type"),
                                  checkboxGroupInput("cell_type", "",
                                                     choices = list(
                                                         "Adipocytes",
                                                         "Basal",
                                                         "Endothelial",
                                                         "Luminal Differentiated",
                                                         "Luminal Progenitors",
                                                         "Stromal"
                                                     ),
                                                     selected = c("Adipocytes",
                                                                  "Basal",
                                                                  "Endothelial",
                                                                  "Luminal Differentiated",
                                                                  "Luminal Progenitors",
                                                                  "Stromal"))
                                  ),
                           column(2,
                                  h4("Stage"),
                                  checkboxGroupInput("stage", "",
                                                     choices = list(
                                                         "Nulliparous",
                                                         "Gestation D5.5",
                                                         "Gestation D9.5",
                                                         "Gestation D14.5",
                                                         "Lactation D2",
                                                         "Laction D5",
                                                         "Lactation D10",
                                                         "Lactation D15",
                                                         "Involution D1",
                                                         "Involution D6",
                                                         "Involution D14"
                                                     ),
                                                     selected = c("Nulliparous",
                                                                  "Gestation D5.5",
                                                                  "Gestation D9.5",
                                                                  "Gestation D14.5",
                                                                  "Lactation D2",
                                                                  "Laction D5",
                                                                  "Lactation D10",
                                                                  "Lactation D15",
                                                                  "Involution D1",
                                                                  "Involution D6",
                                                                  "Involution D14"))),
                           column(2,
                                  h3("Genes"),
                                  textAreaInput("genes", "one per row; leave empty for all genes passing selected filters"),
                                  hr(),
                                  downloadButton("download_data", "Download data")
                       ),
                       fluidRow(
                           column(6,
                              includeHTML("functions/data_explanation.html"))
                           )
                       ))
                
            )
        )
    ),
    tabPanel(
        "PCA",
        fluidPage(
            fluidRow(
                column(2, 
                       h4("Will add download button soon")
                ),
                column(10, 
                       h3("Hybrid data"),
                       plotOutput("hybrid_pca"),
                       h3("Zfp57 data"),
                       plotOutput("zfp57_pca")
                )
            )
        )
    )
)
