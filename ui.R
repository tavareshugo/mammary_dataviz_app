#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

ui <- navbarPage(
    "Mammary ASE Explorer", 
    # allow for enter key to trigger plot button
    header = tags$head(
      tags$script(HTML(
        "$(document).on('keypress', function(e) {
          if (e.which == 13 && $('#gene').is(':focus')) {
            $('#plot').click();
          }
        });"
      ))
    ),
    tabPanel(
        "🏠 Home",
        fluidPage(
          div(
            style = "background-color:#f8f9fa; padding:5px; text-align:center; font-size:small; border-bottom:1px solid #ddd;",
            HTML("If you use data or outputs from this app, please cite: Hanin, G. et al. (2024). <i>Dynamic allelic expression in mouse mammary gland across the adult developmental cycle</i>, <i>bioRxiv</i>, <a href='https://doi.org/10.1101/2024.09.02.610775' target='_blank'>https://doi.org/10.1101/2024.09.02.610775</a>")
            ),
            fluidRow(
                # add a title
                column(12, 
                       h1("Dynamic Allelic Expression in Mouse Mammary Gland Across the Adult Developmental Cycle")
                ),
                column(12, 
                       shinycssloaders::withSpinner(plotOutput("hybrid_pca"))
                ),
                br(),
                column(12, 
                       downloadButton("download_pca_data", "Download PCA data")
                ),
                br(),
                column(6, 
                       includeMarkdown("www/homepage_intro.md")
                ),
                column(6, 
                       includeMarkdown("www/homepage_howto.md")
                )
            )
        )
    ),
    tabPanel(
        "📈 Expression Plots",
        fluidPage(
          div(
            style = "background-color:#f8f9fa; padding:5px; text-align:center; font-size:small; border-bottom:1px solid #ddd;",
            HTML("If you use data or outputs from this app, please cite: Hanin, G. et al. (2024). <i>Dynamic allelic expression in mouse mammary gland across the adult developmental cycle</i>, <i>bioRxiv</i>, <a href='https://doi.org/10.1101/2024.09.02.610775' target='_blank'>https://doi.org/10.1101/2024.09.02.610775</a>")
            ),
          br(),
            fluidRow(
                shinyjs::useShinyjs(),
                column(2, 
                       textInput("gene",
                                 "Gene name or Ensembl ID:"),
                       actionButton("plot", "Plot"),
                       br(), br(),
                       radioButtons(
                         inputId = "expr_type",
                         label = "Expression scale:",
                         choices = c("VST" = "vst",
                                     "Log2 normalised counts" = "log2counts",
                                     "Log2 TPM" = "log2tpm"),
                         selected = "vst"
                       ),
                       actionLink("show_norm_note", "📘 Learn more about normalisation methods")
                ),
                column(9, 
                      #  downloadButton("download_plot_data", "Download plot data"),
                       shinycssloaders::withSpinner(uiOutput("plotted_gene"),
                                                    type = 1),
                       br(),
                       shinycssloaders::withSpinner(uiOutput("download_plot_data_ui"), 
                                                    type = 0),
                       br(),
                       shinycssloaders::withSpinner(plotOutput("imprint_expr"),
                                                    type = 1),
                       shinycssloaders::withSpinner(plotOutput("imprint_isolde"),
                                                    type = 0)
                )
            )
        ),
        tags$footer(
  style = "position:fixed; bottom:0; width:100%; background-color:#f8f9fa; padding:5px; text-align:center; font-size:small; border-top:1px solid #ddd;",
  HTML("Please cite: Hanin, G., Costello, K. R., Tavares, H., AlSulaiti, B., Patel, S., Edwards, C. A., & Ferguson-Smith, A. C. (2024). <i>Dynamic allelic expression in mouse mammary gland across the adult developmental cycle</i>. <i>bioRxiv</i>. <a href='https://doi.org/10.1101/2024.09.02.610775' target='_blank'>https://doi.org/10.1101/2024.09.02.610775</a>")
)
    )
)
