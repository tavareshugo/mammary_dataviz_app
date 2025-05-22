prep_data_for_download <- function(annot,
                                   diffexp,
                                   zfp57_expr,
                                   isolde,
                                   hybrid_expr,
                                   sample_info,
                                   fdr_threshold, 
                                   fc_threshold,
                                   ase_status, 
                                   ase_threshold,
                                   stage,
                                   cell_type,
                                   genes = NULL){
  setProgress(value = 0, message = "Preparing data...")

  # requested genes to subset ----
  if (!is.null(genes) && genes != "") {
    user_genes <- read.table(text = genes, 
                             header = FALSE, 
                             col.names = "name") |> 
      mutate(name = str_trim(name))
    user_genes <- which(tolower(annot$gene) %in% tolower(user_genes$name) | 
                          tolower(annot$name) %in% tolower(user_genes$name))
    user_genes <- annot$gene[user_genes]
  } else {
    user_genes <- annot$gene
  }

  # tidy annotation ----
  annot <- annot |> 
    collect() |>
    summarise(name = paste(name, collapse = "/"), .by = "gene")
    
  # Genes to retain ----
  gene_set_degs <- diffexp |> 
    collect() |>
    # filter
    filter(padj <= fdr_threshold & abs(log2FoldChange) >= log2(fc_threshold) &
             stage %in% stage & cell_type %in% cell_type) |> 
    distinct(gene) |> pull(gene)
  gene_set_isolde <- isolde |> 
    collect() |>
    filter(abs(diff_prop) >= ase_threshold & status %in% ase_status &
             stage %in% stage & cell_type %in% cell_type) |> 
    distinct(gene) |> pull(gene)
  gene_set <- intersect(gene_set_degs, gene_set_isolde)
  if(length(user_genes) > 0){
    gene_set <- intersect(gene_set, user_genes)
  }

  # if(!is.null(genes)){
  #   genes <- str_squish(genes)
  #   annot <- annot |> filter(gene %in% genes)
  # }
  
  # differential expression ----
  setProgress(value = 0.2, message = "Preparing data...")
  
  # filters
  diffexp_out <- diffexp |> 
    collect() |>
    # filter
    filter(gene %in% gene_set & 
             padj <= fdr_threshold & 
             abs(log2FoldChange) >= log2(fc_threshold) &
             stage %in% stage & 
             cell_type %in% cell_type) |>
    # add annotation
    inner_join(annot, by = "gene") |> 
    # calculate fold-change
    mutate(fc = 2^log2FoldChange) |> 
    # tidy columns
    select(gene, name, cell_type, stage, fc, 
           lfc = log2FoldChange, lfcSE, pvalue, padj) |> 
    arrange(gene, cell_type, stage)
  
  # zfp57 expression ----
  setProgress(value = 0.4, message = "Preparing data...")
  
  zfp57_expr_out <- zfp57_expr |> 
    collect() |>
    filter(gene %in% gene_set) |> 
    # right_join(diffexp_out |> distinct(gene), by = "gene") |> 
    left_join(sample_info, by = "sample") |> 
    filter(stage %in% stage & cell_type %in% cell_type) |> 
    inner_join(annot, by = "gene") |> 
    select(gene, name, cell_type, stage, cross, animal_id, genotype, 
           normalised_expression = expr) |> 
    arrange(gene, cell_type, stage, genotype)
  
  
  # Hybrid expression ----
  setProgress(value = 0.6, message = "Preparing data...")
  
  hybrid_expr_out <- hybrid_expr |> 
    collect() |>
    filter(gene %in% gene_set) |> 
    # right_join(diffexp_out |> distinct(gene), by = "gene") |> 
    left_join(sample_info, by = "sample") |> 
    filter(stage %in% stage & cell_type %in% cell_type) |> 
    inner_join(annot, by = "gene") |> 
    select(gene, name, cell_type, stage, cross, animal_id, 
           normalised_expression = expr) |> 
    arrange(gene, cell_type, stage)
  
  # Isolde test ----
  setProgress(value = 0.8, message = "Preparing data...")
  
  isolde_out <- isolde |> 
    collect() |>
    # retain only DEGs
    # right_join(diffexp_out |> distinct(gene), by = "gene") |> 
    # filter
    filter(gene %in% gene_set & 
             abs(diff_prop) >= ase_threshold & 
             status %in% ase_status &
             stage %in% stage & 
             cell_type %in% cell_type) |> 
    # add annotation
    inner_join(annot, by = "gene") |> 
    # tidy columns
    select(gene, name, cell_type, stage, maternal_minus_paternal = diff_prop, 
           status) |> 
    arrange(gene, cell_type, stage)
  
  # return result ----
  setProgress(value = 1, message = "Preparing data...")
  
  list(`ZFP57 differential expression` = diffexp_out |> as.data.frame(),
       `ZFP57 normalised expression` = zfp57_expr_out |> as.data.frame(),
       `Hybrid normalised expression` = hybrid_expr_out |> as.data.frame(),
       `Hybrid Isolde analysis` = isolde_out |> as.data.frame())
}
