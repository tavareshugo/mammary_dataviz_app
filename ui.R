#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

ui <- fluidPage(

    # Application title
    titlePanel("Mammary Gland Transcriptomes"),

    fluidRow(
        column(3, 
               textInput("gene",
                         "Gene name:",
                         value = "Zfp57"),
               actionButton("plot", "Plot")
        ),
        column(9, 
               plotOutput("zfp57_expr"),
               plotOutput("imprint_expr"),
               plotOutput("imprint_isolde")
        )
    )

    # sidebarLayout(
    #     sidebarPanel(
    #         textInput("gene",
    #                   "Gene name:",
    #                   value = "Zfp57"),
    #         actionButton("plot", "Plot")
    #     ),
    # 
    #     mainPanel(
    #        plotOutput("zfp57_expr"),
    #        plotOutput("imprint_expr"),
    #        plotOutput("imprint_isolde")
    #     )
    # )
)
