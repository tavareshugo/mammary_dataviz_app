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
        "🏠 Home",
        fluidPage(
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
            fluidRow(
                shinyjs::useShinyjs(),
                column(2, 
                       textInput("gene",
                                 "Gene name or Ensembl ID:"),
                       actionButton("plot", "Plot")
                ),
                column(9, 
                      #  downloadButton("download_plot_data", "Download plot data"),
                       shinycssloaders::withSpinner(uiOutput("plotted_gene"),
                                                    type = 1),
                       br(),
                       uiOutput("download_plot_data_ui"),
                       br(),
                       shinycssloaders::withSpinner(plotOutput("imprint_expr"),
                                                    type = 0),
                       shinycssloaders::withSpinner(plotOutput("imprint_isolde"),
                                                    type = 0)
                )
            )
        )
    )
)
