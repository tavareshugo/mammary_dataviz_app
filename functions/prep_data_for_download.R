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
  # requested genes to subset ----
  user_genes <- read.table(text = genes, 
                           header = FALSE, 
                           col.names = "name") %>% 
    mutate(name = str_trim(name))
  user_genes <- which(tolower(annot$gene) %in% tolower(user_genes$name) | 
                        tolower(annot$name) %in% tolower(user_genes$name))
  user_genes <- annot$gene[user_genes]
  
  # tidy annotation ----
  annot <- annot %>% 
    group_by(gene) %>% 
    summarise(name = paste(name, collapse = "/")) %>% 
    ungroup()

  # Genes to retain ----
  gene_set_degs <- diffexp %>% 
    # filter
    filter(padj <= fdr_threshold & abs(log2FoldChange) >= log2(fc_threshold) &
             stage %in% stage & cell_type %in% cell_type) %>% 
    distinct(gene) %>% pull(gene)
  gene_set_isolde <- isolde %>% 
    filter(abs(diff_prop) >= ase_threshold & status %in% ase_status &
             stage %in% stage & cell_type %in% cell_type) %>% 
    distinct(gene) %>% pull(gene)
  gene_set <- intersect(gene_set_degs, gene_set_isolde)
  if(length(user_genes) > 0){
    gene_set <- intersect(gene_set, user_genes)
  }
  
  
  # if(!is.null(genes)){
  #   genes <- str_squish(genes)
  #   annot <- annot %>% filter(gene %in% genes)
  # }
  
  # differential expression ----
  # filters
  diffexp_out <- diffexp %>% 
    # filter
    filter(gene %in% gene_set & 
             padj <= fdr_threshold & 
             abs(log2FoldChange) >= log2(fc_threshold) &
             stage %in% stage & 
             cell_type %in% cell_type) %>%
    # add annotation
    inner_join(annot, by = "gene") %>% 
    # calculate fold-change
    mutate(fc = 2^log2FoldChange) %>% 
    # tidy columns
    select(gene, name, cell_type, stage, fc, 
           lfc = log2FoldChange, lfcSE, pvalue, padj) %>% 
    arrange(gene, cell_type, stage)
  
  # zfp57 expression ----
  zfp57_expr_out <- zfp57_expr %>% 
    filter(gene %in% gene_set) %>% 
    # right_join(diffexp_out %>% distinct(gene), by = "gene") %>% 
    left_join(sample_info, by = "sample") %>% 
    filter(stage %in% stage & cell_type %in% cell_type) %>% 
    inner_join(annot, by = "gene") %>% 
    select(gene, name, cell_type, stage, cross, animal_id, genotype, 
           normalised_expression = expr) %>% 
    arrange(gene, cell_type, stage, genotype)
  
  
  # Hybrid expression ----
  hybrid_expr_out <- hybrid_expr %>% 
    filter(gene %in% gene_set) %>% 
    # right_join(diffexp_out %>% distinct(gene), by = "gene") %>% 
    left_join(sample_info, by = "sample") %>% 
    filter(stage %in% stage & cell_type %in% cell_type) %>% 
    inner_join(annot, by = "gene") %>% 
    select(gene, name, cell_type, stage, cross, animal_id, 
           normalised_expression = expr) %>% 
    arrange(gene, cell_type, stage)
  
  # Isolde test ----
  isolde_out <- isolde %>% 
    # retain only DEGs
    # right_join(diffexp_out %>% distinct(gene), by = "gene") %>% 
    # filter
    filter(gene %in% gene_set & 
             abs(diff_prop) >= ase_threshold & 
             status %in% ase_status &
             stage %in% stage & 
             cell_type %in% cell_type) %>% 
    # add annotation
    inner_join(annot, by = "gene") %>% 
    # tidy columns
    select(gene, name, cell_type, stage, maternal_minus_paternal = diff_prop, 
           status) %>% 
    arrange(gene, cell_type, stage)
  
  # return result ----
  list(`ZFP57 differential expression` = diffexp_out,
       `ZFP57 normalised expression` = zfp57_expr_out,
       `Hybrid normalised expression` = hybrid_expr_out,
       `Hybrid Isolde analysis` = isolde_out)
  
}
