[![Build Status](https://api.travis-ci.com/ncats/RaMP-DB.svg?branch=sqlite)](https://travis-ci.com/github/ncats/RaMP-DB)

# New!  RaMP-DB 3.0!

RaMP-DB 3.0 is now released and includes an updated backend database with expanded annotations for >250,000 metabolites, ~34,000 genes/proteins, >122,000 pathways, and > 66,000 reactions.  Annotations include biological pathways, chemical classes and structures (for metabolites only), ontologies (metabolites only), and enzyme-metabolite relationships based on chemical reactions. Annotations are drawn from HMDB, KEGG (through HMDB), Lipid-MAPS, WikiPathways, Reactome, CheBI, PFOCR, RefMet, and Rhea reaction database. 

This R package allows users to interface with this up-do-date and comprehensive resource.  Functionalities include 1) Queries: simple and batch queries for pathways, ontologies, chemical annotations, and reaction-level gene-metabolite relationships; 2) Enrichment Analysis: pathway and chemical enrichment analyses; 3) Exploration of Results: interactive visualizations of metabolite/gene annotations.

The code used to build the backend RaMP-DB database is freely available at [https://github.com/ncats/RaMP-Backend](https://github.com/ncats/RaMP-Backend).

Please [click here to view our latest manuscript](https://pubmed.ncbi.nlm.nih.gov/36373969/).

# Web Interface
Our new revamped web interface can be found at [https://rampdb.nih.gov/](https://rampdb.nih.gov/).  The code for the front-end is publicly available at [https://github.com/ncats/ncats-frontend-library/tree/main/apps/ramp-client](https://github.com/ncats/ncats-frontend-library/tree/main/apps/ramp-client).

# APIs
API access is now available [here](https://rampdb.nih.gov/api).

# Why RaMP-DB (Relational Database of Metabolomic Pathways)

The purpose of RaMP-DB is to provide a publicly available database that integrates metabolite and gene/protein biological, 
chemical, reaction, and ontology information from multiple sources. The database structure and data is available as an SQLite database file and it is directly downloaded when using the RaMP-DB package.
Please see the Installation Instructions for further information.
Please note that this project is in continuous development and we appreciate any feedback (through email below or create a GitHub issue)!  Thanks to those that have provided some thus far.

## Contact Info:
For any questions or feedback, please send us a note at [NCATSRaMP@mail.nih.gov](NCATSRaMP@mail.nih.gov). 

If you find a bug, please submit an issue through this GitHub repo. 

## Vignette
Detailed instructions for installing RaMP locally are below.  We've also put together a vignette to get you started on the analyses.  Click here for [vignette](https://ncats.github.io/RaMP-DB/Updated_RaMP_Vignette.html).

## Citation
If you use RaMP-DB, please cite the following work:

Braisted J, Patt A, Tindall C, Sheils T, Neyra J, Spencer K, Eicher T, Mathé EA. RaMP-DB 2.0: a renovated knowledgebase for deriving biological and chemical insight from metabolites, proteins, and genes. Bioinformatics. 2023 Jan 1;39(1):btac726. doi: 10.1093/bioinformatics/btac726. PMID: 36373969; PMCID: PMC9825745.
To access, [click here](https://pubmed.ncbi.nlm.nih.gov/36373969/)

Zhang, B., et al., RaMP: A Comprehensive Relational Database of Metabolomics Pathways for Pathway Enrichment Analysis of Genes and Metabolites. Metabolites, 2018. 8(1). PMID: 29470400; PMCID: PMC5876005; DOI: 10.3390/metabo8010016
To access, [click here](https://www.mdpi.com/2218-1989/8/1/16)

## Installation Instructions
In order to use this R package locally, you will need to install the R code under this repository.

*Special Note:*
There is incompatibility (reported here: https://stat.ethz.ch/pipermail/bioc-devel/2023-October/020003.html) between 
the version of BiocFileCache installed using BiocManager (2.8.0) and the actual latest version (2.10.1).  The 
latter is needed to be compatible with other dependencies in RaMP-DB.  To install the latest version, you will 
need to download the source file from Bioconductor (https://bioconductor.org/packages/release/bioc/html/BiocFileCache.html), 
then install using the install.packages() function.  For a Mac, this looks like this:

```
install.packages('https://www.bioconductor.org/packages/3.18/bioc/src/contrib/Archive/BiocFileCache/BiocFileCache_2.10.1.tar.gz', repos=NULL, method='libcurl', dependencies=TRUE)
```

### Install and load the RaMP package 
You can install this package directly from GitHub using the install_github() function available through the [devtools package](https://cran.r-project.org/web/packages/devtools/index.html). In the R Console, type the following:
```R
# Locally install RaMP
install.packages("devtools")
library(devtools)
install_github("ncats/RAMP-DB")

# Load the package
library(RaMP)

# initializes the RaMP database object, downloading and caching the latest SQLite database
# if no version already exists in local cache.
rampDB <- RaMP()

# note that you can use the following method to check database versions hosted in your
# computer's local cache and databases that are available to download in our remote repository.
RaMP::listAvailableRaMPDbVersions()

# using that list of available RaMP DB versions, one can specify the database version to use
# if the selected version is not available on your computer, but is in our remote repository at GitHub,
# the SQLite DB file will be automatically downloaded into local file cache.
# RaMP is using the BiocFileCache package to manage a local file cache.
rampDB <- RaMP(version = "2.5.4")

```

### Important Notes

When gene or metabolite ids are input for queries, IDs should be prepended with their database of origin, e.g. kegg:C02712, hmdb:HMDB04824, or CAS:2566-39-4. The list of metabolite or gene/protien IDs may be of mixed source. Remember to include the colon in the prefix. The id prefixes that are currently included in RaMP are: 

| Analyte Type | ID Prefix Types |
|--------------|-----------------|
| Metabolites | hmdb, pubchem, chebi, chemspider, kegg, CAS, LIPIDMAPS, swisslipids, lipidbank, wikidata, plantfa, kegg_glycan |
| Genes/Proteins | ensembl, entrez, gene_symbol, uniprot, hmdb, ncbiprotein, EN, wikidata, chebi

The following RaMP functions can be used to list all represented id prefix types.
```
rampDB <- RaMP()
RaMP::getPrefixesFromAnalytes(db = rampDB, analyteType = 'metabolite')
RaMP::getPrefixesFromAnalytes(db = rampDB, analyteType = 'gene')
```

## Current Authors and Testers
* **Ewy Mathé** - ewy.mathe@nih.gov
* **Andrew Patt** - andy.patt@nih.gov
* **Khyati Mehta** - khyati.meta@nih.gov
* **Adam Tisch** - adam.tisch@nih.gov
* **Keith Kelleher** - keith.kelleher@nih.gov
* **Tim Sheils** - tim.sheils@nih.gov


## Previous Authors/Testers
* **John Braisted** - 
* **Tara Eicher** - 
* **Kyle Spencer** - 
* **Cole Tindall** - 
* **Bofei Zhang** - [Bofei5675](https://github.com/Bofei5675)
* **Shunchao Wang** - 
* **Rohith Vanam** - 
* **Jorge Neyra** - [Jorgeso](https://github.com/jorgeso)
