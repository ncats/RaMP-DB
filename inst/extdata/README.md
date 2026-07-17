# extdata File Descriptions

This directory contains example data files to aid in learning how to use RaMP

## metLinkR
These files were used for benchmarking the metLinkR package in the [original publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC12053952/). See the [github repository](https://github.com/ncats/MetLinkR) for more information, including a vignette.

The "manifest" file in this directory (HarmInputFiles.csv) is properly formatted for MetLinkR to analyze the input metabolomics data sets (the other 5 csvs in the directory). MetLinkR recognizes common names, HMDB IDs, KEGG IDs, LIPID MAPS IDs, PubChem IDs, and ChEBI IDs for human metabolites.

## CCLE data

The Cancer Cell Line Encyclopedia (CCLE) is an initiative from the Broad Institute to perform large-scale characterize of ~1000 cancer cell lines. Our example data is a list of analyte IDs which were found to be significantly different between cancer types in a subset of the data. 

Herein, we took a subset of data and focused on primary tumors of: 

- Acute Myeloid Leukaemia within the HAEMATOPOIETIC AND LYMPHOID TISSUE (n = 30)
- Glioma within the CENTRAL NERVOUS SYSTEM (n = 45)
- FIBROBLAST (n = 10)
- OESOPHAGUS (n = 25)

within 3 data files: 
- RNAseq: CCLE_RNAseq_rsem_genes_tpm_20180929.csv
- Proteomics: CCLE_RPPA_20181003.csv
- Metablomics: CCLE_metabolomics_20190502.csv

We identified analytes differentiating cancer types by: 
 1) An ANOVA was performed across all 4 cancer types and p-values were FDR adjusted
 2) Tukey HSD (Honestly Significant Difference) was performed for each analyte
     - A post hoc statistical test used after a significant ANOVA to determine which specific group means differ from each other
 3) Analytes which have a significant FDR adjusted ANOVA p-value 
    AND 
    significant p-adj between >4 comparisons in the Tukey's test were selected
