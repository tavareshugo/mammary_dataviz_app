library(DESeq2)
library(biomaRt)
library(tidyverse)
library(nanoparquet)

# Read dds objects -----

# hybrid dataset
hybrid <- readRDS("../mammary_cellsort_hybrid_rnaseq/results/DESeqDataSet/dds_gene_regular.rds")

# zfp57 dataset
zfp57 <- readRDS("../mammary_cellsort_zfp57_rnaseq/results/diffexp/dds.rds")

# simplify sample names
colnames(hybrid) <- paste0("hybrid_", 1:ncol(hybrid))
colnames(zfp57) <- paste0("zfp57_", 1:ncol(zfp57))



# Prepare sample information ----------------------------------------------

hybrid_samples <- hybrid |> 
  colData() |> 
  as_tibble(rownames = "sample") |> 
  distinct(sample, cell_type, stage, cross, animal_id) |> 
  mutate_all(as.character)

zfp57_samples <- zfp57 |> 
  colData() |> 
  as_tibble(rownames = "sample") |> 
  distinct(sample, cell_type, stage, genotype, animal_id) |> 
  mutate_all(as.character)

sample_info <- bind_rows(hybrid_samples, zfp57_samples)

# tidy up factor variables
sample_info <- sample_info |> 
  mutate(cell_type = case_when(cell_type == "endo" ~ "endothelial",
                               cell_type == "lumhi" ~ "luminal progenitors", 
                               cell_type == "lumlo" ~ "luminal differentiated",
                               TRUE ~ cell_type),
         stage = case_when(stage == "Gest 9.5" ~ "gestation d9.5",
                           stage == "Lac 2" ~ "lactation d2",
                           TRUE ~ stage)) |> 
  mutate(stage = str_replace(stage, "gestation", "pregnancy")) |>
  mutate(cell_type = str_to_title(cell_type), 
         stage = str_to_title(stage),
         cross = str_to_upper(cross), 
         genotype = str_to_upper(genotype)) |> 
  mutate(stage = factor(stage, 
                        levels = c("Nulliparous", "Pregnancy D5.5", 
                                   "Pregnancy D9.5", "Pregnancy D14.5", 
                                   "Lactation D2", "Lactation D5",
                                   "Lactation D10", "Lactation D15", 
                                   "Involution D1", "Involution D6", 
                                   "Involution D14")))


# Prepare expression information -------------

hybrid_expr <- map_dfr(
  c("vst", "logcounts", "tpm"),
  function(assay_name) {
    hybrid |>
      assay(assay_name) |>
      as_tibble(rownames = "gene") |>
      pivot_longer(-gene, names_to = "sample", values_to = "expr") |>
      mutate(assay = assay_name)
  }
) |>
  pivot_wider(names_from = "assay", values_from = "expr") |> 
  mutate(log2tpm = log2(tpm + 1)) |> 
  rename(log2counts = logcounts) |>
  select(-tpm)

zfp57_expr <- map_dfr(
  c("vst", "logcounts", "tpm"),
  function(assay_name) {
    zfp57 |>
      assay(assay_name) |>
      as_tibble(rownames = "gene") |>
      pivot_longer(-gene, names_to = "sample", values_to = "expr") |>
      mutate(assay = assay_name)
  }
) |>
  pivot_wider(names_from = "assay", values_from = "expr") |> 
  mutate(log2tpm = log2(tpm + 1)) |> 
  rename(log2counts = logcounts) |>
  select(-tpm)


# Differential Expression -------------------------------------------------

diffexp <- read_csv("../mammary_cellsort_zfp57_rnaseq/results/diffexp/diffexp_wt_vs_ko.csv")
diffexp <- diffexp |> 
  mutate(cell_type = case_when(cell_type == "endo" ~ "endothelial",
                               cell_type == "lumhi" ~ "luminal progenitors", 
                               cell_type == "lumlo" ~ "luminal differentiated",
                               TRUE ~ cell_type),
         stage = case_when(timepoint == "t0" ~ "nulliparous",
                           timepoint == "t2" ~ "gestation d9.5",
                           timepoint == "t4" ~ "lactation d2",
                           TRUE ~ NA_character_)) |> 
    mutate(stage = str_replace(stage, "gestation", "pregnancy")) |>
    mutate(cell_type = str_to_title(cell_type), 
           stage = str_to_title(stage)) |> 
    mutate(stage = factor(stage, 
                          levels = c("Nulliparous", "Pregnancy D5.5", 
                                     "Pregnancy D9.5", "Pregnancy D14.5", 
                                     "Lactation D2", "Lactation D5",
                                     "Lactation D10", "Lactation D15", 
                                     "Involution D1", "Involution D6", 
                                     "Involution D14")))


# Isolde data ------------

# read data
isolde <- read_csv("../mammary_cellsort_hybrid_rnaseq/results/ase_isolde/isolde_all.csv")

# tidy
isolde <- isolde |>
  rename(gene = name) |>
  mutate(cell_type = case_when(cell_type == "luminalp" ~ "luminal progenitors",
                               cell_type == "luminald" ~ "luminal differentiated",
                               TRUE ~ cell_type)) |>
  inner_join(hybrid |> colData() |> as_tibble(rownames = "sample") |>
               distinct(cell_type, timepoint, stage),
             by = c("cell_type", "timepoint")) |> 
  mutate(stage = str_replace(stage, "gestation", "pregnancy")) |>
  mutate(cell_type = str_to_title(cell_type), 
         stage = str_to_title(stage)) |> 
  mutate(stage = factor(stage, 
                        levels = c("Nulliparous", "Pregnancy D5.5", 
                                   "Pregnancy D9.5", "Pregnancy D14.5", 
                                   "Lactation D2", "Lactation D5",
                                   "Lactation D10", "Lactation D15", 
                                   "Involution D1", "Involution D6", 
                                   "Involution D14"))) |> 
  select(gene, cell_type, stage, everything())


# Allele-specific counts --------------

# allele-specific counts
dds_isolde <- readRDS("../mammary_cellsort_hybrid_rnaseq/results/DESeqDataSet/dds_gene_allele.rds")

# estimate size factors for normalisation (ISoLDE recommends doing this)
dds_isolde <- estimateSizeFactors(dds_isolde)

# normalised counts table
isolde_cts <- dds_isolde[unique(isolde$gene), ] |> 
  counts(normalized = TRUE) |> 
  as_tibble(rownames = "gene") |>
  pivot_longer(-gene, names_to = "sample", values_to = "cts")

# add sample info and sample ID
isolde_cts <- isolde_cts |> 
  left_join(dds_isolde |> colData() |> as_tibble(rownames = "sample"),
            by = "sample") |> 
  # mutate(sample = sample |> str_remove("_B$") |> str_remove("_C$")) |> 
  # test to check everything is there
  # anti_join(hybrid_samples, by = c("cell_type", "stage", "cross", "animal_id"))
  select(-sample) |> 
  full_join(hybrid_samples, by = c("cell_type", "stage", "cross", "animal_id"))

isolde_cts <- isolde_cts |> 
  select(gene, sample, allele_strain, allele_parent, cts)


# Gene Annotation -----------

# get gene annotation - this was mapped to GRCm38 assembly
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",
                   host="nov2020.archive.ensembl.org")

# gene information
gene_annot <- getBM(attributes = c("ensembl_gene_id",
                                   "external_gene_name",
                                   "external_synonym"),
                    mart = mart)

gene_annot <- gene_annot |> 
  pivot_longer(-ensembl_gene_id) |> 
  dplyr::select(gene = ensembl_gene_id, name = value) |> 
  distinct() |> 
  filter(name != "")

isoform_annot <- getBM(attributes = c("ensembl_gene_id",
                                      "ensembl_transcript_id",
                                      "external_transcript_name"),
                       mart = mart)
isoform_annot <- isoform_annot |> 
  distinct() |> 
  rename(gene = ensembl_gene_id, transcript = ensembl_transcript_id, 
         transcript_name = external_transcript_name)


# Isoform analysis ----

# deseq dataset
isoform <- readRDS("../mammary_cellsort_hybrid_rnaseq/results/DESeqDataSet/dds_isoform_regular.rds")
colnames(isoform) <- paste0("hybrid_", 1:ncol(isoform))

isoform_expr <- isoform |> 
  assay("vst") |> 
  as_tibble(rownames = "transcript") |> 
  pivot_longer(-transcript, names_to = "sample", values_to = "expr")

# isolde analysis
isolde_isoform <- read_csv("../mammary_cellsort_hybrid_rnaseq/results/ase_isolde_isoform/isolde_isoform_all.csv")

isolde_isoform <- isolde_isoform |> 
  rename(transcript = name) |>
  mutate(cell_type = case_when(cell_type == "luminalp" ~ "luminal progenitors",
                               cell_type == "luminald" ~ "luminal differentiated",
                               TRUE ~ cell_type)) |>
  inner_join(isoform |> colData() |> as_tibble(rownames = "sample") |>
               distinct(cell_type, timepoint, stage),
             by = c("cell_type", "timepoint")) |> 
  mutate(stage = str_replace(stage, "gestation", "pregnancy")) |>
  mutate(cell_type = str_to_title(cell_type), 
         stage = str_to_title(stage)) |> 
  mutate(stage = factor(stage, 
                        levels = c("Nulliparous", "Pregnancy D5.5", 
                                   "Pregnancy D9.5", "Pregnancy D14.5", 
                                   "Lactation D2", "Lactation D5",
                                   "Lactation D10", "Lactation D15", 
                                   "Involution D1", "Involution D6", 
                                   "Involution D14"))) |> 
  select(transcript, cell_type, stage, everything())



# Dimensionality reduction ------

# PCA on hybrid data
hybrid_topn <- order(rowVars(assay(hybrid, "vst")), decreasing = TRUE)[1:1000]
hybrid_pca <- hybrid[hybrid_topn, ] |>
  assay("vst") |>
  t() |> #scale() |>
  prcomp()

# PCA on zfp57 data
zfp57_topn <- order(rowVars(assay(zfp57, "vst")), decreasing = TRUE)[1:1000]
zfp57_pca <- zfp57[zfp57_topn, ] |>
  assay("vst") |>
  t() |> #scale() |>
  prcomp()


# Save data ---------

nanoparquet::write_parquet(
  sample_info, "data/sample_info.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  zfp57_expr, "data/zfp57_normalised_expression.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  hybrid_expr, "data/hybrid_normalised_expression.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  isolde_cts, "data/hybrid_allele_cts.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  isolde, "data/hybrid_isolde.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  gene_annot, "data/gene_annotation.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  diffexp, "data/zfp57_differential_expression.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  isoform_annot, "data/isoform_annotation.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  isoform_expr, "data/hybrid_isoform_normalised_expression.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
nanoparquet::write_parquet(
  isolde_isoform, "data/hybrid_isoform_isolde.parquet",
  compression = "gzip",
  options = nanoparquet::parquet_options(compression_level = 9)
)
zfp57_pca |> saveRDS("data/zfp57_pca.rds")
hybrid_pca |> saveRDS("data/hybrid_pca.rds")
