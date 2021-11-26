library(shiny)
library(tidyverse)
library(patchwork)
library(writexl)
library(ggfortify)
theme_set(theme_minimal(base_size = 16))

# source functions
source("functions/prep_data_for_download.R")

# colouring for plots
stage_colours <- c("Nulliparous" = "grey",
                   "Gestation D5.5" = "#c6dbef",
                   "Gestation D9.5" = "#6baed6",
                   "Gestation D14.5" = "#08519c",
                   "Lactation D2" = "#bae4b3",
                   "Lactation D5" = "#74c476",
                   "Lactation D10" = "#31a354",
                   "Lactation D15" = "#006d2c",
                   "Involution D1" = "#fdd0a2",
                   "Involution D6" = "#fd8d3c",
                   "Involution D14" = "#a63603")

cell_colours <- RColorBrewer::brewer.pal(6, "Dark2")
names(cell_colours) <- c("Adipocytes", "Basal", "Endothelial", "Luminal Differentiated", "Luminal Progenitors", "Stromal")

cross_colours <- c("BC" = "black", "CB" = "brown")

genotype_colours <- c("WT" = "black", "KO" = "brown")


# load data
zfp57_expr <- readRDS("data/zfp57_normalised_expression.rds")
hybrid_expr <- readRDS("data/hybrid_normalised_expression.rds")
isolde <- readRDS("data/hybrid_isolde.rds")
annot <- readRDS("data/gene_annotation.rds")
sample_info <- readRDS("data/sample_info.rds")
diffexp <- readRDS("data/zfp57_differential_expression.rds")
zfp57_pca <- readRDS("data/zfp57_pca.rds")
hybrid_pca <- readRDS("data/hybrid_pca.rds")
hybrid_isoform_expr <- readRDS("data/hybrid_isoform_normalised_expression.rds")
isolde_isoform <- readRDS("data/hybrid_isoform_isolde.rds")
isoform_annot <- readRDS("data/isoform_annotation.rds")


# Define server logic 
server <- function(input, output, session) {
  
  # Plot panel -----
  
  # check if things should be plotted
  # target_gene <- eventReactive(c(input$plot, input$isoform_plot), {
  #   annot %>% 
  #     filter(gene == toupper(input$gene) | toupper(name) == toupper(input$gene)) %>% 
  #     distinct(gene)      
  # })
  
  # react to both buttons being pushed
  target_gene <- reactiveValues()
  
  observeEvent(input$plot, {
    # update the other input text box
    updateTextInput(session, inputId = "isoform_gene", value = input$gene)
    target_gene$gene <- annot %>% 
      filter(gene == toupper(input$gene) | toupper(name) == toupper(input$gene)) %>% 
      distinct(gene)      
  })
  
  observeEvent(input$isoform_plot, {
    # update the other input text box
    updateTextInput(session, inputId = "gene", value = input$isoform_gene)
    target_gene$gene <- annot %>% 
      filter(gene == toupper(input$isoform_gene) | toupper(name) == toupper(input$isoform_gene)) %>% 
      distinct(gene)
  })
  
  # check if gene is valid
  valid_target_gene <- function(x = target_gene$gene){
    if(is.null(x)){
      "Choose a gene to plot."
    } else if(nrow(x) == 0){
      "No genes found."
    } else if (nrow(x) > 1) {
      paste(c("Multiple genes found matching that name:", x$gene), collapse = " ")
    } else if (nrow(x) == 1) {
      NULL
    } else {
      "Unknown error."
    }
  }
  
  # Expression plots ----
  
  # information about which gene is being plotted
  output$plotted_gene <- renderText({
    validate(valid_target_gene())
    plotted_gene <- target_gene$gene %>% 
      left_join(annot, by = "gene") %>% 
      group_by(gene) %>% 
      summarise(name = paste(name, collapse = "/")) %>% ungroup()
    paste0("Showing: ", plotted_gene$name, " (", plotted_gene$gene, ")")
  })
  
  # zfp57 expression plots
  output$zfp57_expr <- renderPlot({
    validate(valid_target_gene())
    
    p1 <- target_gene$gene %>% 
      left_join(diffexp, by = "gene") %>% 
      mutate(padj = ifelse(is.na(padj), 1, padj)) %>% 
      ggplot(aes(stage, -log2FoldChange)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(aes(group = cell_type, colour = cell_type)) +
      geom_point(aes(shape = padj < 0.05)) +
      labs(x = "", y = "LFC(KO/WT)", colour = "Cell", shape = "FDR < 5%") +
      scale_colour_manual(values = cell_colours) +
      scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- target_gene$gene %>% 
      left_join(zfp57_expr, by = "gene") %>% 
      left_join(sample_info, by = "sample") %>% 
      ggplot(aes(stage, expr)) +
      ggbeeswarm::geom_quasirandom(aes(colour = genotype), 
                                   dodge.width = 0.5, size = 2,
                                   groupOnX = TRUE) +
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
    
    target_gene$gene %>% 
      left_join(hybrid_expr) %>% 
      left_join(sample_info, by = "sample") %>% 
      ggplot(aes(stage, expr, colour = cell_type)) +
      ggbeeswarm::geom_quasirandom(dodge.width = 0.5, size = 2, groupOnX = TRUE) +
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
    
    target_gene$gene %>% 
      left_join(isolde, by = "gene") %>% 
      drop_na(diff_prop) %>% 
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
  
  # Isoform expression ----
  
  # information about which gene/isoform is being plotted
  output$plotted_isoform <- renderText({
    validate(valid_target_gene())
    plotted_gene <- target_gene$gene %>% 
      left_join(annot, by = "gene") %>% 
      group_by(gene) %>% 
      summarise(name = paste(name, collapse = "/")) %>% ungroup()
    paste0("Showing: ", plotted_gene$name, " (", plotted_gene$gene, ")")
  })
  output$ensembl_isoform_link <- renderUI({
    validate(valid_target_gene())
    url <- a("Ensembl Link", 
             href = paste0("https://nov2020.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=",
                    target_gene$gene))
    tagList(url)
  })
  
  # this is used to adjust the height of the plots
  n_isoforms <- reactive({
    if(is.null(target_gene$gene)) return(1)
    
    n_isoforms <- target_gene$gene %>%
      left_join(isoform_annot, by = "gene") %>%
      distinct(transcript) %>%
      nrow()
    return(n_isoforms)
  })

  # expression in hybrid data
  output$isoform_expression <- renderPlot({
    
    validate(valid_target_gene())
    
    # target_gene$gene %>% 
    #   left_join(isoform_annot, by = "gene") %>% 
    #   inner_join(hybrid_isoform_expr, by = "transcript") %>% 
    #   left_join(sample_info, by = "sample") %>% 
    #   ggplot(aes(stage, expr, colour = cell_type)) +
    #   ggbeeswarm::geom_quasirandom(dodge.width = 0.5, size = 2, groupOnX = TRUE) +
    #   geom_line(stat = "summary", fun = "median", aes(group = 1),
    #             size = 1) +
    #   facet_grid(transcript_name ~ cell_type) +
    #   labs(title = "Expression in hybrid dataset", 
    #        colour = "Cell",
    #        x = "", y = "Normalised Expression") +
    #   scale_colour_manual(values = cell_colours) +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
    #         panel.border = element_rect(fill = NA, colour = "black"))
    
    target_gene$gene %>% 
      left_join(isoform_annot, by = "gene") %>% 
      inner_join(hybrid_isoform_expr, by = "transcript") %>% 
      left_join(sample_info, by = "sample") %>% 
      group_by(transcript_name, stage, cell_type) %>% 
      summarise(expr = median(expr)) %>% 
      ggplot(aes(stage, transcript_name, fill = expr)) +
      geom_tile() + 
      facet_grid( ~ cell_type) +
      labs(title = "Expression in hybrid dataset", 
           subtitle = "Heatmap shows median expression across 4 replicates in each condition (cell type + stage)",
           fill = "Normalised\nExpression",
           x = "", y = "") +
      scale_fill_viridis_c() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_rect(fill = NA, colour = "black"))
    
  }, height = function() max(400, as.integer(ceiling(60 * n_isoforms()/2))))
  
  # isolde plots
  output$isoform_ase <- renderPlot({
    validate(valid_target_gene())
    
    target_gene$gene %>% 
      left_join(isoform_annot, by = "gene") %>% 
      inner_join(isolde_isoform, by = "transcript") %>% 
      drop_na(diff_prop) %>% 
      ggplot(aes(stage, diff_prop)) +
      geom_hline(yintercept = 0) +
      geom_line(aes(colour = cell_type, group = cell_type),
                size = 1) +
      geom_point(aes(colour = cell_type), size = 2) +
      facet_wrap(~ transcript_name, ncol = 4) +
      labs(title = "Allele-specific bias", 
           colour = "Cell",
           x = "", y = "Difference (M - P)") +
      scale_y_continuous(limits = c(-1, 1)) +
      scale_colour_manual(values = cell_colours) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_rect(fill = NA, colour = "black"))
    
  }, height = function() max(400, as.integer(ceiling(200 * n_isoforms()/8))))

  # Data panel ----
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("mammary_gland_expression", ".xlsx", sep = "")
    },
    content = function(file) {
      withProgress({
        out <- prep_data_for_download(annot = annot,
                                      diffexp = diffexp,
                                      zfp57_expr = zfp57_expr,
                                      isolde = isolde,
                                      hybrid_expr = hybrid_expr,
                                      sample_info = sample_info,
                                      fdr_threshold = input$fdr_zfp57,
                                      fc_threshold = input$fc_threshold,
                                      ase_status = input$isolde_status,
                                      ase_threshold = input$ase_bias_threshold,
                                      stage = input$stage,
                                      cell_type = input$cell_type,
                                      genes = input$genes)
      })
      write_xlsx(out, file)
    }
  )

  # PCA panel ----
  output$hybrid_pca <- renderPlot({
    wrap_plots(
      autoplot(hybrid_pca,
               data = sample_info[match(rownames(hybrid_pca$x), sample_info$sample), ],
               x = 1, y = 2, colour = "cell_type") +
        scale_colour_manual(values = cell_colours, name = "Cell"),
      autoplot(hybrid_pca,
               data = sample_info[match(rownames(hybrid_pca$x), sample_info$sample), ],
               x = 1, y = 2, colour = "stage") +
        scale_colour_manual(values = stage_colours, name = "Stage"),
      autoplot(hybrid_pca,
               data = sample_info[match(rownames(hybrid_pca$x), sample_info$sample), ],
               x = 1, y = 2, colour = "cross") +
        scale_colour_manual(values = cross_colours, name = "Cross"),
      nrow = 1
    ) + plot_annotation(title = "PCA on Hybrid Samples") & coord_equal()
  })
  
  output$zfp57_pca <- renderPlot({
    
    wrap_plots(
      autoplot(zfp57_pca,
               data = sample_info[match(rownames(zfp57_pca$x), sample_info$sample), ],
               x = 1, y = 2, colour = "cell_type") +
        scale_colour_manual(values = cell_colours, name = "Cell"),
      autoplot(zfp57_pca,
               data = sample_info[match(rownames(zfp57_pca$x), sample_info$sample), ],
               x = 1, y = 2, colour = "stage") +
        scale_colour_manual(values = stage_colours, name = "Stage"),
      autoplot(zfp57_pca,
               data = sample_info[match(rownames(zfp57_pca$x), sample_info$sample), ],
               x = 1, y = 2, colour = "genotype") +
        scale_colour_manual(values = genotype_colours, name = "Genotype"),
      nrow = 1
    ) + plot_annotation(title = "PCA on ZFP57 Samples") & coord_equal()
    
  })
  
  output$download_pca_data <- downloadHandler(
    filename = function() {
      paste("mammary_gland_pca", ".xlsx", sep = "")
    },
    content = function(file) {
      withProgress({
        zfp57_scores <- zfp57_pca$x %>% 
          as_tibble(rownames = "sample") %>% 
          select(sample, PC1:PC10) %>% 
          left_join(sample_info, by = "sample") %>% 
          select(cell_type, stage, genotype, animal_id, matches("PC"))
        
        zfp57_variance <- tibble(variance = zfp57_pca$sdev^2) %>% 
          mutate(PC = paste0("PC", 1:n()),
                 pct_variance = variance/sum(variance)*100) %>% 
          filter(PC %in% paste0("PC", 1:10)) %>% 
          select(PC, variance, pct_variance)
        
        hybrid_scores <- hybrid_pca$x %>% 
          as_tibble(rownames = "sample") %>% 
          select(sample, PC1:PC10) %>% 
          left_join(sample_info, by = "sample") %>% 
          select(cell_type, stage, genotype, animal_id, matches("PC"))
        
        hybrid_variance <- tibble(variance = hybrid_pca$sdev^2) %>% 
          mutate(PC = paste0("PC", 1:n()),
                 pct_variance = variance/sum(variance)*100) %>% 
          filter(PC %in% paste0("PC", 1:10)) %>% 
          select(PC, variance, pct_variance)
        
      })
      write_xlsx(
        list(`Hybrid PCA` = hybrid_scores,
             `Hybrid PCA variance` = hybrid_variance,
             `ZFP57 PCA` = zfp57_scores,
             `ZFP57 PCA variance` = zfp57_variance), 
        file
        )
    }
  )
  
  
}
