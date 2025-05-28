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
    # https://stackoverflow.com/a/56771353
    header = tags$head(
        tags$script(
            "$(document).on('shiny:inputchanged', function(event) {
          if (event.name != 'changed') {
            Shiny.setInputValue('changed', event.name);
          }
        });"
        )
    ),
    tabPanel(
        "Gene Plots",
        fluidPage(
            fluidRow(
                shinyjs::useShinyjs(),
                column(2, 
                       textInput("gene",
                                 "Gene name or Ensembl ID:"),
                       actionButton("plot", "Plot"),
                       br(), br(),
                       downloadButton("download_plot_data", "Download Plot Data")
                ),
                column(9, 
                      #  downloadButton("download_plot_data", "Download plot data"),
                       shinycssloaders::withSpinner(uiOutput("plotted_gene"),
                                                    type = 1),
                       shinycssloaders::withSpinner(plotOutput("imprint_expr"),
                                                    type = 0),
                       shinycssloaders::withSpinner(plotOutput("imprint_isolde"),
                                                    type = 0)
                )
            )
        )
    ),
    tabPanel(
        "PCA",
        fluidPage(
            fluidRow(
                column(2, 
                       downloadButton("download_pca_data", "Download PCA data")
                ),
                column(10, 
                       shinycssloaders::withSpinner(plotOutput("hybrid_pca"))
                )
            )
        )
    )
)
