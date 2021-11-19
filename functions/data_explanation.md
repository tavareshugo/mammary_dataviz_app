<details><summary>**Explanation of filters** (click to expand)</summary>

- **FDR** is the false-discovery rate (i.e. corrected p-value) threshold to use for filtering differentially expressed genes in the ZFP57 WT-vs-KO comparisons. A 5% threshold is commonly used to consider a gene differentially expressed. 
<span style="color:#003300"><b>FDR = 1 is equivalent to no filter.</b></span>

- **Fold change** threshold is the minimum ratio of expression to consider between the two genotypes. If a threshold of 2x is chosen, it will select genes where either WT/KO > 2 or KO/WT > 2. 
<span style="color:#003300"><b>FC = 1 is equivalent to no filter.</b></span>

- **Isolde status** is the significance determined by _Isolde_. "Undetermined" is when there was some ambiguity in the data (e.g. a gene with a bias of ~0.6 or with high variation across replicates). "Filtered out" are genes which did not pass the filtering criteria used by _Isolde_. 
<span style="color:#003300"><b>Selecting all status is equivalent to no filter.</b></span>

- **Allele bias threshold** is the threshold used for the difference between maternal and paternal alleles. If a threshold of 0.8 is chosen, it will pick genes biased in either direction (i.e. either `maternal - paternal > 0.8` or `maternal - paternal < -0.8`).
<span style="color:#003300"><b>Bias = 0 is equivalent to no filter.</b></span>
</details>

<br>

**NOTE:** Filters apply cumulatively. 
For example, the default will return genes with FDR < 0.05 _and_ fold-change > 2 _and_ biased expression _and_ allele bias > 0.7 _and_ any cell type _and_ any stage. 

----

**Output is an Excel file with 4 worksheets.** 
This includes data from the two datasets, with two tables each:

- Zfp57 KO/WT data
    - Normalised expression in each sample
    - Differential expression analysis between KO/WT
- Hybrid data
  - Normalised expression in each sample
  - Allele-specific expression analysis (Isolde)

