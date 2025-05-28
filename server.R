library(shiny)
library(duckplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(writexl)
library(ggfortify)
theme_set(theme_minimal(base_size = 16))

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

# load data
hybrid_expr <- read_parquet_duckdb("data/hybrid_normalised_expression.parquet")
isolde <- read_parquet_duckdb("data/hybrid_isolde.parquet")
annot <- read_parquet_duckdb("data/gene_annotation.parquet")
sample_info <- read_parquet_duckdb("data/sample_info.parquet")
hybrid_pca <- readRDS("data/hybrid_pca.rds")


# Define server logic 
server <- function(input, output, session) {
  
  # Plot panel -----
  
  target_gene <- reactiveValues()
  
  observeEvent(input$plot, {
    target_gene$gene <- annot |> 
      filter(gene == toupper(input$gene) | toupper(name) == toupper(input$gene)) |> 
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
  output$plotted_gene <- renderUI({
    validate(valid_target_gene())
    
    matches <- target_gene$gene |>
      left_join(annot, by = "gene")
    
    all_names <- unique(matches$name)
    user_input <- isolate(toupper(input$gene))

    # Determine primary name
    matched_names <- all_names[toupper(all_names) == user_input]
    is_ensembl_id <- grepl("^ENSMUSG\\d+$", user_input)

    primary <- if (is_ensembl_id || length(matched_names) == 0) {
      all_names[1]
    } else {
      matched_names[1]
    }

    aliases <- setdiff(all_names, primary)

    alias_text <- ""
    if (length(aliases) == 1) {
      alias_text <- paste0(", also known as ", aliases)
    } else if (length(aliases) > 1 && length(aliases) <= 3) {
      alias_text <- paste0(", also known as ", paste(aliases, collapse = ", "))
    } else if (length(aliases) > 3) {
      alias_text <- paste0(", also known as ",
                          paste(head(aliases, 3), collapse = ", "),
                          ", and ", length(aliases) - 3, " others")
    }

    url <- paste0(
      '<a href="https://nov2020.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=',
      target_gene$gene[1],
      '" target="_blank">Ensembl Link</a>'
    )

    HTML(paste0("<b>Showing:</b> ", primary, alias_text, " (", url, ")"))
  })

  
  # expression in hybrid data
  output$imprint_expr <- renderPlot({
    if (is.null(target_gene$gene) || nrow(target_gene$gene) != 1) return(NULL)

    target_gene$gene |> 
      left_join(hybrid_expr) |> 
      left_join(sample_info, by = "sample") |> 
      mutate(stage = factor(stage, levels = names(stage_colours))) |>
      ggplot(aes(stage, expr, colour = cell_type)) +
      ggbeeswarm::geom_quasirandom(dodge.width = 0.5, size = 2, groupOnX = TRUE) +
      geom_line(stat = "summary", fun = "median", aes(group = 1),
                size = 1) +
      facet_grid(~ cell_type) +
      labs(title = "Expression", 
           colour = "Cell",
           x = "", y = "Normalised Expression") +
      scale_colour_manual(values = cell_colours) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.position = "none")
  })
  
  # isolde plots
  output$imprint_isolde <- renderPlot({
    if (is.null(target_gene$gene) || nrow(target_gene$gene) != 1) return(NULL)
    
    target_gene$gene |> 
      left_join(isolde, by = "gene") |> 
      filter(!is.na(diff_prop)) |> 
      mutate(stage = factor(stage, levels = names(stage_colours))) |>
      ggplot(aes(stage, diff_prop)) +
      geom_hline(yintercept = 0) +
      geom_line(aes(colour = cell_type, group = cell_type),
                size = 1) +
      geom_point(aes(colour = cell_type), size = 2) +
      labs(title = "Allele-specific bias (ISoLDE)", 
           colour = "Cell",
           x = "", y = "Difference (M - P)") +
      scale_y_continuous(limits = c(-1, 1)) +
      scale_colour_manual(values = cell_colours) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$download_plot_data <- downloadHandler(
    filename = function() {
      paste0("expression_", target_gene$gene$gene[1], ".xlsx")
    },
    content = function(file) {
      withProgress(message = "Preparing plot data...", {
        
        req(target_gene$gene)
        validate(valid_target_gene())
        
        # Expression data
        expr_data <- target_gene$gene |>
          left_join(hybrid_expr) |>
          left_join(sample_info, by = "sample") |> 
          select(-genotype)
        
        # Isolde data
        isolde_data <- target_gene$gene |>
          left_join(isolde, by = "gene")
        
        # Annot info
        gene_info <- target_gene$gene |>
          left_join(annot, by = "gene")
        
        write_xlsx(
          list(
            Gene_Info = gene_info,
            Expression = expr_data,
            Isolde = isolde_data
          ),
          path = file
        )
      })
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
  
  output$download_pca_data <- downloadHandler(
    filename = function() {
      paste("mammary_gland_pca", ".xlsx", sep = "")
    },
    content = function(file) {
      withProgress({
        
        hybrid_scores <- hybrid_pca$x |> 
          as_tibble(rownames = "sample") |> 
          select(sample, PC1:PC10) |> 
          left_join(sample_info, by = "sample") |> 
          select(cell_type, stage, genotype, animal_id, matches("PC"))
        
        hybrid_variance <- tibble(variance = hybrid_pca$sdev^2) |> 
          mutate(PC = paste0("PC", 1:n()),
                 pct_variance = variance/sum(variance)*100) |> 
          filter(PC %in% paste0("PC", 1:10)) |> 
          select(PC, variance, pct_variance)
        
      })
      write_xlsx(
        list(`Hybrid PCA` = hybrid_scores,
             `Hybrid PCA variance` = hybrid_variance), 
        file
      )
    }
  )
}
