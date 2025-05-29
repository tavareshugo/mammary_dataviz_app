**TL;DR:**

TPM is fine for comparing genes within the same sample, but not across samples. 
For between-sample comparisons, we recommend using **variance-stabilised counts** or **DESeq2-normalised counts**.

----

**Why not TPM for comparing across samples?**

TPM (and its cousin FPKM) is a popular way to normalise RNA-seq data by transcript length and sequencing depth. 
While this works well for comparing different genes within the same sample, it performs poorly for comparing gene expression across samples, especially when transcript composition differs (e.g. due to tissue type, treatment, stage, etc.).

This is because TPM is a relative measure: if one gene's expression increases, others may appear to decrease even if they haven’t changed at all. 
This makes TPM less robust and misleading for between-sample comparisons.
Figure 1 of <a href="https://doi.org/10.1093/bib/bbx008" target="_blank">Evans et al. 2018</a> illustrates this issue clearly. 


**What should I use instead?**

RNA-seq analysis packages recommend using normalisation strategies that are robust to differences in transcript composition between samples. 
Two popular methods include:

* Median-of-ratios (used by `DESeq2`)
* TMM (trimmed mean of M-values) (used by `edgeR`)

These approaches adjust for sample-specific biases, producing more reliable comparisons across experimental conditions.


**What about variance?**

RNA-seq count data typically shows _heteroscedasticity_, whereby lowly expressed genes tend to have more variability relative to their mean. 
This can, in some cases, affect visualisation and clustering. 
To address this, `DESeq2` provides a _variance stabilising transformation (VST)_ that equalises variance across the expression range, making downstream analyses more robust.


**What’s included in this app:**

We provide three options for visualising expression:

- **VST (Variance Stabilising Transform)**: Recommended. Transforms normalised counts to reduce mean-variance dependence. Useful for PCA, clustering, and visualisation.
- **Log2 normalised counts**: Based on `DESeq2`'s median-of-ratios normalisation, followed by log transformation: log2(normcts + 1). 
- **Log2 TPM**: Log-transformed transcripts per million: log2(TPM + 1). Included due to their popularity, but interpret between-sample comparisons with caution.

---- 

**Further reading**

- Evans, C., Hardin, J., & Stoebel, D.M. (2018). Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions. *Briefings in Bioinformatics*, 19(5), 776–792. <a href="https://doi.org/10.1093/bib/bbx008" target="_blank">doi.org/10.1093/bib/bbx008</a>
- Zhao, S., Ye, Z., & Stanton, R. (2020). Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. *RNA*, 26, 903–909. <a href="https://doi.org/10.1261/rna.074922.120" target="_blank">doi.org/10.1261/rna.074922.120</a>
- Zhao, Y., Li, M.C., Konaté, M.M., et al. (2021). TPM, FPKM, or Normalized Counts? A Comparative Study of Quantification Measures for the Analysis of RNA-seq Data from the NCI Patient-Derived Models Repository. *Journal of Translational Medicine*, 19, 269. <a href="https://doi.org/10.1186/s12967-021-02936-w" target="_blank">doi.org/10.1186/s12967-021-02936-w</a>
