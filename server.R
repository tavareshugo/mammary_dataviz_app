library(shiny)
library(tidyverse)
library(patchwork)
library(writexl)
theme_set(theme_minimal(base_size = 16))


# colouring for plots
stage_colours <- c("Nulliparous" = "grey",
                   "Gestation D5.5" = "#a6cee3",
                   "Gestation D9.5" = "#9ecae1",
                   "gestation D14.5" = "#3182bd",
                   "Lactation D5" = "#e5f5e0",
                   "Lactation D10" = "#a1d99b",
                   "Lactation D15" = "#31a354",
                   "Involution D1" = "#fee6ce",
                   "Involution D6" = "#fdae6b",
                   "Involution D14" = "#e6550d",
                   "Lactation D2" = "#a1d99b")

cell_colours <- RColorBrewer::brewer.pal(6, "Dark2")
names(cell_colours) <- c("Adipocytes", "Basal", "Endothelial", "Luminal Differentiated", "Luminal Progenitors", "Stromal")

cross_colours <- c("NC" = "black", "CB" = "brown")

genotype_colours <- c("WT" = "black", "KO" = "brown")


# load data
zfp57_expr <- readRDS("data/zfp57_normalised_expression.rds")
hybrid_expr <- readRDS("data/hybrid_normalised_expression.rds")
isolde <- readRDS("data/hybrid_isolde.rds")
annot <- readRDS("data/gene_annotation.rds")
sample_info <- readRDS("data/sample_info.rds")
diffexp <- readRDS("data/zfp57_differential_expression.rds")


# Define server logic required to draw a histogram
server <- function(input, output) {

  # check if things should be plotted
  target_gene <- eventReactive(input$plot, {
    annot %>% 
      filter(gene == toupper(input$gene) | toupper(name) == toupper(input$gene)) %>% 
      distinct(gene)
  })
  
  # check if gene is valid
  valid_target_gene <- function(x = target_gene()$gene){
    if(length(x) == 0){
      "No genes found."
    } else if (length(x) > 1) {
      paste(c("Multiple genes found matching that name:", x), collapse = " ")
    } else if (length(x) == 1) {
      NULL
    } else {
      "Unknown error."
    }
  }
  
  
  # Plots ----
  
  # zfp57 expression plots
  output$zfp57_expr <- renderPlot({
    validate(valid_target_gene())
    
    p1 <- target_gene() %>% 
      left_join(diffexp, by = "gene") %>% 
      ggplot(aes(stage, -log2FoldChange)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(aes(group = cell_type, colour = cell_type)) +
      geom_point(aes(shape = padj < 0.05)) +
      labs(x = "", y = "LFC(KO/WT)", colour = "Cell", shape = "FDR < 5%") +
      scale_colour_manual(values = cell_colours) +
      scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- target_gene() %>% 
      left_join(zfp57_expr, by = "gene") %>% 
      left_join(sample_info, by = "sample") %>% 
      ggplot(aes(stage, expr)) +
      ggbeeswarm::geom_quasirandom(aes(colour = genotype), 
                                   dodge.width = 0.5, size = 2) +
      facet_grid(~ cell_type) +
      labs(title = "Expression in Zfp57 WT/KO samples", 
           colour = "Genotype",
           x = "", y = "Normalised Expression") +
      scale_colour_manual(values = genotype_colours) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    wrap_plots(p1, p2) +
      plot_layout(ncol = 2, widths = c(1, 5))
  })
  
  # expression in hybrid data
  output$imprint_expr <- renderPlot({
    validate(valid_target_gene())
    
    target_gene() %>% 
      left_join(hybrid_expr) %>% 
      left_join(sample_info, by = "sample") %>% 
      ggplot(aes(stage, expr, colour = cell_type)) +
      ggbeeswarm::geom_quasirandom(dodge.width = 0.5, size = 2) +
      geom_line(stat = "summary", fun = "median", aes(group = 1),
                size = 1) +
      facet_grid(~ cell_type) +
      labs(title = "Expression in hybrid dataset", 
           colour = "Cell",
           x = "", y = "Normalised Expression") +
      scale_colour_manual(values = cell_colours) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # isolde plots
  output$imprint_isolde <- renderPlot({
    validate(valid_target_gene())
    
    target_gene() %>% 
      left_join(isolde, by = "gene") %>% 
      ggplot(aes(stage, diff_prop)) +
      geom_hline(yintercept = 0) +
      geom_line(aes(colour = cell_type, group = cell_type),
                size = 1) +
      geom_point(aes(colour = cell_type), size = 2) +
      labs(title = "Allele-specific bias in hybrid dataset", 
           colour = "Cell",
           x = "", y = "Difference (M - P)") +
      scale_y_continuous(limits = c(-1, 1)) +
      scale_colour_manual(values = cell_colours) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Download data ----
  output$download_data <- downloadHandler(
    filename = function() {
      paste(target_gene()$gene, ".xlsx", sep = "")
    },
    content = function(file) {
      out <- list(isolde = target_gene() %>% left_join(isolde, by = "gene"),
                  hybrid_expression = target_gene() %>% 
                    left_join(hybrid_expr) %>% 
                    left_join(sample_info, by = "sample"),
                  zfp57_expression = target_gene() %>% 
                    left_join(zfp57_expr, by = "gene") %>% 
                    left_join(sample_info, by = "sample"))
      write_xlsx(out, file)
    }
  )
  
  observe({
    shinyjs::toggleState("download_data", is.null(valid_target_gene()))
  })
  
}
