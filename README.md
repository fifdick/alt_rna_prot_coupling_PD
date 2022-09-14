# Altered transcriptome-proteome coupling indicates aberrant proteostasis in Parkinson’s disease 

This repository holds all code and data to reproduce results, figures and and supplementary data presented in "Altered transcriptome-proteome coupling indicates aberrant proteostasis in Parkinson’s disease".
The paper is currently available on [medrxiv](https://www.medrxiv.org/content/10.1101/2021.03.18.21253875v1) (not peer reviewed).


## Prerequisites

* [R >= 3.6.0](https://www.r-project.org/)
* All R packages used in the analysis are listed at the top of each Rmd file
* Many of these are bioconductor packages so it helps to have BiocManager installed:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(<"packagename">)
```

## How to..  

### Pathway enrichment analysis

* The paper discusses results of pathway enrichment analysis that were performed on ranked gene lists.
* The ranking is based on the Pearson correlation coefficient (gene-wise, across-samples and within groups (HA, YG, PD). 
* The ranking takes into account the difference in correlation between groups (see paper methods).
* The genelist ranking (based on the correlation coefficients given in `./results/rds/gene_cl.rds`) are generated in `./pea.Rmd`.
* The pathway enrichment analysis is performed with [fgsea](https://github.com/ctlab/fgsea) in `./pea.Rmd`. 
* Requirements for this analysis are the `.rds` files in `./results/rds/` and the GO_simplified and KEGG annotations for fgsea in `./referenceData/`. R packages needed are listed at the beginning of the `./pea.Rmd` file.  
* Permutation analysis of these results are also performed in `pea.Rmd`, requirements for this is the permuted (simulated) dataset in `rdsData/perms/gene_cl_perm_9_15_4_5000.rds` and the pathway annotations. Permutation analysis was performed on the analyses where we compared correlations between groups and on the analyses where we did pathway enrichment analysis on the correlation coefficients in each group separately.   
* Would not advice to knit this file. Hasnt been edited appropriately or tested. The execution of permutation analyses is very time consuming. The purpose is mainly to document code.

### Reproduce Paper Figures

* The main figures (excluding figure 1, which is a schematic figure) can be reproduced with `./analysis.Rmd`.
* The code loads `.rds` files need for the analysis form `./rdsData/`. These are generated in `./preprocess.Rmd` as described below.
* The figures are available in `./result_figs/Draft/` and `./result_figs/Supp/`.
* The `analysis.Rmd` also creates and saves the `gene_cl.rds` file to `./results/rds/` (needed in `./pea.Rmd`).
* Would not advice to knit this file. Hasnt been edited appropriately or tested. The purpose is mainly to document code.

### Rerun analysis from scratch

* The rawest data available are transcript-level counts in `./salmonOut` and protein intensities in `./rdsData/rawData.rds`.
* With these and the metadata in `./rdsData/info.rds` `./rdsData/info_rna.rds` and `./rdsData/info_rna_pd.rds` all results can be reproduced.
* The file `./preprocess.Rmd`, loads the required data and:
  * aggregates transcript-level counts from `./salmonOut` to gene-level using [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) and the annotation file `./referenceData/gencode.v32.annotation.gtf` (needs to be unzipped).
  * filters out zeros and applies batch correction to the proteomics data in `./rdsData/rawData.rds`
  * saves all generated `.rds` files to `./rdsData`.
* Utility functions (also the batch correction function) are loaded from `./functions.R`.
* Next steps would be to go through `./analysis.Rmd` and then `./pea.Rmd`
* Would not advice to knit `./analysis.Rmd`. Hasnt been edited appropriately or tested. The purpose is mainly to document code.

## Content

* `./preprocess.Rmd` Code to preprocess data (Transcript counts and raw protein intensities)
* `./analysis.Rmd` Code to reproduce paper figures 
* `./pea.Rmd`, `./pea.html` Code to run pathway enrichment analysis and knitted html to browse through results
* `./referenceData/` Holds reference data needed for the analysis. Remember to unzip compressed files (before running `./preprocess.Rmd`)
* `./salmonOut/` Output from [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) for each sample in the analysis. Details on parameters of the salmon command see methods section of the paper. Sample_ids (Salmon dir. names) correspond to RNAseq_id (see also supplementary File 1, or './rdsData/info.rds`.
* `./rdsData/` 
  * `info.rds` Holds all sample metadata. (The proteomics data uses reporter.intensity.id as sample identifier, this can be mapped to RNAseq_id for the RNASeq data).
  * `info_rna.rds` and `info_rna_pd.rds` Holds sample metadata for RNA data (including RIN) for the groups HA, YG and PD respectively. 
  * `cpmAll.rds` and `cpmAll_pd.rds` Holds gene-level counts in counts per million (CPM), after filtering. These files are created in `./preprocess.Rmd`
  * `rawData.rds` Holds protein intensities before filtering and before normalization. This file is used in `./preprocess.Rmd` to generate the batch corrected matrix which is then used in the analysis.
* `./result_figs` All paper figures (some might have been edited with graphics software afterwards and thus dont resemble exactly the same figure as in the paper)
* `./results` Holds the rds files for the pathway enrichment analysis (i.e. correlation coefficients for the groups (YG, HA and PD).
* Each `.Rmd` file loads functions from `./functions.R`. 

## Version info 

### R sessionInfo() output  

```
> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8          LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8           LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8       LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8          LC_NAME=en_GB.UTF-8
 [9] LC_ADDRESS=en_GB.UTF-8        LC_TELEPHONE=en_GB.UTF-8
[11] LC_MEASUREMENT=en_GB.UTF-8    LC_IDENTIFICATION=en_GB.UTF-8

attached base packages:
 [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets
 [9] methods   base

other attached packages:
 [1] fgsea_1.21.2                ggbiplot_0.55
 [3] plyr_1.8.6                  infer_1.0.0
 [5] pbapply_1.5-0               ggplotify_0.0.5
 [7] ggrepel_0.8.2               pheatmap_1.0.12
 [9] stringr_1.4.0               knitr_1.29
[11] readr_1.3.1                 preprocessCore_1.48.0
[13] tidyr_1.1.1                 gridExtra_2.3
[15] readxl_1.3.1                ggfortify_0.4.10
[17] DESeq2_1.26.0               SummarizedExperiment_1.16.1
[19] DelayedArray_0.12.3         BiocParallel_1.20.1
[21] matrixStats_0.56.0          ggpubr_0.4.0
[23] tibble_3.0.3                mixOmics_6.19.4
[25] lattice_0.20-40             MASS_7.3-51.5
[27] boot_1.3-24                 cowplot_1.1.0
[29] ensembldb_2.10.2            AnnotationFilter_1.10.0
[31] GenomicFeatures_1.38.2      AnnotationDbi_1.48.0
[33] Biobase_2.46.0              GenomicRanges_1.38.0
[35] GenomeInfoDb_1.22.1         IRanges_2.20.2
[37] S4Vectors_0.24.4            BiocGenerics_0.32.0
[39] scales_1.1.1                RColorBrewer_1.1-2
[41] dplyr_1.0.2                 colorspace_1.4-1
[43] xlsx_0.6.4.2                kableExtra_1.1.0
[45] colormap_0.1.4              dendextend_1.14.0
[47] igraph_1.2.5                ermineR_1.0.1.9000
[49] patchwork_1.0.1             ComplexHeatmap_2.2.0
[51] viridis_0.6.0               viridisLite_0.4.0
[53] ggplot2_3.3.2               reshape2_1.4.4
[55] magrittr_1.5                nvimcom_0.9-102

loaded via a namespace (and not attached):
  [1] tidyselect_1.1.0         RSQLite_2.2.0            htmlwidgets_1.5.1
  [4] munsell_0.5.0            withr_2.2.0              rstudioapi_0.11
  [7] ggsignif_0.6.0           rJava_0.9-13             GenomeInfoDbData_1.2.2
 [10] bit64_4.0.2              vctrs_0.3.2              generics_0.0.2
 [13] xfun_0.16                BiocFileCache_1.10.2     R6_2.4.1
 [16] clue_0.3-57              locfit_1.5-9.4           bitops_1.0-6
 [19] gridGraphics_0.5-0       assertthat_0.2.1         nnet_7.3-13
 [22] gtable_0.3.0             rlang_0.4.10             genefilter_1.68.0
 [25] GlobalOptions_0.1.2      splines_3.6.3            rtracklayer_1.46.0
 [28] rstatix_0.6.0            lazyeval_0.2.2           gemmaAPI_2.0.3.9000
 [31] broom_0.7.0              checkmate_2.0.0          BiocManager_1.30.10
 [34] abind_1.4-5              backports_1.1.8          Hmisc_4.4-1
 [37] tools_3.6.3              ellipsis_0.3.1           Rcpp_1.0.5
 [40] base64enc_0.1-3          progress_1.2.2           zlibbioc_1.32.0
 [43] purrr_0.3.4              RCurl_1.98-1.2           prettyunits_1.1.1
 [46] rpart_4.1-15             openssl_1.4.2            GetoptLong_1.0.2
 [49] haven_2.3.1              cluster_2.1.4            data.table_1.13.0
 [52] RSpectra_0.16-0          openxlsx_4.1.5           circlize_0.4.10
 [55] ProtGenerics_1.18.0      hms_0.5.3                xlsxjars_0.6.1
 [58] evaluate_0.14            xtable_1.8-4             XML_3.99-0
 [61] rio_0.5.16               jpeg_0.1-8.1             shape_1.4.4
 [64] compiler_3.6.3           biomaRt_2.42.1           ellipse_0.4.2
 [67] V8_3.2.0                 crayon_1.3.4             unixtools_0.1-1
 [70] htmltools_0.5.1.1        corpcor_1.6.9            Formula_1.2-3
 [73] geneplotter_1.64.0       DBI_1.1.0                dbplyr_1.4.4
 [76] rappdirs_0.3.1           Matrix_1.2-18            car_3.0-9
 [79] forcats_0.5.0            pkgconfig_2.0.3          rvcheck_0.1.8
 [82] GenomicAlignments_1.22.1 foreign_0.8-75           xml2_1.3.2
 [85] rARPACK_0.11-0           annotate_1.64.0          webshot_0.5.2
 [88] XVector_0.26.0           rvest_0.3.6              digest_0.6.25
 [91] Biostrings_2.54.0        fastmatch_1.1-0          rmarkdown_2.3
 [94] cellranger_1.1.0         htmlTable_2.0.1          curl_4.3
 [97] Rsamtools_2.2.3          rjson_0.2.20             lifecycle_0.2.0
[100] jsonlite_1.7.0           carData_3.0-4            askpass_1.1
[103] pillar_1.4.6             httr_1.4.2               survival_3.1-8
[106] glue_1.4.1               zip_2.1.0                png_0.1-7
[109] bit_4.0.4                stringi_1.4.6            blob_1.2.1
[112] latticeExtra_0.6-29      memoise_1.1.0
```
 
