## Summary

This application provides access to allele-specific expression (ASE) data from 480 RNA-seq libraries generated from reciprocal F1 hybrid mice (C57BL/6J × CAST/EiJ) sampled across several developmental stages of the mouse mammary gland (<a href="https://doi.org/10.1101/2024.09.02.610775" target="_blank">Hanin et al. 2025</a>).

Tissues were collected at ten stages: nulliparous and 3 timepoints each at gestation, lactation and involution.
For each animal, cells were sorted into six major cell populations using FACS: basal, luminal progenitor, luminal differentiated, stromal, endothelial and adipocyte cells.

RNA was extracted from flow-sorted cell populations and sequenced using Illumina. 
Expression quantification was performed using Salmon, and allele-specific expression was assessed using a custom pipeline incorporating CAST/EiJ-specific genome annotations and the ISoLDE package. 
Normalisation and differential expression analysis was done using `DESeq2`, and PCA was applied to variance-stabilised counts.

