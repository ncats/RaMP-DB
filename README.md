[![Build Status](https://api.travis-ci.com/ncats/RaMP-DB.svg?branch=master&status=passed)](https://travis-ci.com/github/ncats/RaMP-DB)

# New!  RaMP app is accessible via a server (no installation needed!).
Coming Soon!

# RaMP - Relational Database of Metabolomic Pathways

The purpose of RaMP is to provide a publicly available database that integrates metabolite and gene biological pathways from multiple sources. Currently, we have integrated information from HMDB, KEGG, Reactome, and WikiPathways. The relational structure of RaMP enables complex and batch queries.  To facilitate its usage, we have created this R shiny web application that includes a user-friendly R Shiny web application.  Please note that this project is in continuous development and we certainly appreciate your feedback (through [our](https://github.com/ncats/RaMP-DB) GitHub site). More details can be found in <a href="http://www.mdpi.com/2218-1989/8/1/16" target="_blank">our manuscript</a>.

Also note that we are working on a server version of RaMP so that users do not have to install anything on their local machines.  Stay tuned!

## Basic Features:
The R packages and associated app performs some complex queries (e.g. retrieve all metabolites and/or genes that belong to a user input pathway or list of pathways).  It also performs pathway enrichment analysis given a list of metabolites and/or genes. Run the app to get further details.

Last date of dump file update: 03/02/2018

## Vignette
Detailed instructions for installing RaMP locally are below.  We've also put together a vignette to get you started on the analyses.  Click here for [vignette](https://ncats.github.io/RaMP-DB/RaMP_Vignette.html).


## Citation
If you use RaMP, please cite the following work:

Zhang, B., et al., RaMP: A Comprehensive Relational Database of Metabolomics Pathways for Pathway Enrichment Analysis of Genes and Metabolites. Metabolites, 2018. 8(1).

PMID: 29470400; PMCID: PMC5876005; DOI: 10.3390/metabo8010016

To access, [click here](https://www.mdpi.com/2218-1989/8/1/16)

## Installation Instructions
In order to use the web application, you will need the following:
* The R code under this repo
* The mysql dump file that contains the RaMP database (in the folder inst/extdata/)

If you would like to know how to build RaMP database from scratch, please check another GitHub site at [RaMP-BackEnd](https://github.com/ncats/RaMP-BackEnd)

### MySQL set-up
**Warning:** RaMP will not readily work with the new MySQL version (8.X.X), but is fully functional with MySQL version 5.7. We will be working toward updating RaMP to work with the newest MySQL version.

RaMP requires that MySQL and the RaMP database be set up on the machine that you will be running the app from.
To download MySQL, you can go to the [MySQL Downloads page](https://www.mysql.com/downloads/)

When installing, you will be prompted to create a password for the user "root", or it will create one automatically for you.  **Importantly, remember your MySQL password!**  You will need to get into mysql and to pass it as an argument to the RaMP R shiny web application.

If you want to reset your password , you can go to [MySQL References 5.7 - How to reset root password ] (https://dev.mysql.com/doc/refman/5.7/en/resetting-permissions.html)

Please note that you will need administrator privileges for this step..

### Creating the database locally
Once your MySQL environment is in place, creating the RaMP database locally is trivial.
First, launch MySQL and create the database:
```
> mysql -u root -p
mysql> create database ramp;
mysql> exit;
```

Here, we are naming the database RaMP but you can use any name you'd like.  It is worth noting though that the R package assumes that the name of the database is 'ramp' by default.  So if you change the name, remember to pass that name as arguments in the R package functions.

Second, populate the named database with the mysql dump file (which you can get from  inst/extdata/ramp180302.sql):
```
> mysql -u root -p ramp < ramp180302.sql
```

You're done!

Your "ramp" database should contain the following 8 tables:
1. analyte
2. analyehasontology
3. analytehaspathway
4. analytesynonym
5. catalyzed
6. ontology
7. pathway
8. source

If you want to explore this in MySQL, you can try:
```
mysql -u root -p
use ramp;
show tables;
select * from analytesynonym where synonym = "glucose";
```

### Install and load the RaMP package 
You can install this package directly from GitHub using the install_github() function available through the devtools package. In the R Console, type the following:
```R
# Locally install RaMP
install.packages("devtools")
library(devtools)
install_github("ncats/RAMP-DB")

# Load the package
library(RaMP)

# Set up your connection to the RaMP2.0 database:
pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
```

Note that prior to using RaMP functions, users much establish required parameters to appropriately connect to your local database (if you are not using the web app).  This step is simplified by a single function call (last line in the above code snippet).

If the username is different then root, then specify the username in the "username" parameter.  Similarly, if the name of the database is different than "ramp2", then specify the "dbname" parameter.

### Important Notes

If you reinstall the latest version of the RaMP package, be sure to also install the latest version of the mysql RaMP dump file.  

Also, when gene or metabolite ids are input for queries, IDs should be prepended with their database of origin, e.g. kegg:C02712, hmdb:HMDB04824, or CAS:2566-39-4. The list of metabolite or gene/protien IDs may be of mixed source. Remember to include the colon in the prefix. The id prefixes that are currently included in RaMP are: 

| Analyte Type | ID Prefix Types |
|--------------|-----------------|
| Metabolites | CAS, chebi, chemspider, enzymeNomenclature, hmdb, kegg, LIPIDMAPS, pubchem |
| Genes/Proteins |ensembl, entrez, uniprot |


## Current Authors
* **John Braisted** - john.braisted@nih.gov
* **Ewy Mathé** - ewy.mathe@nih.gov
* **Jorge Neyra** -jorge.neyra@nih.gov
* **Andrew Patt** - andy.patt@nih.gov
* **Tim Sheils** - tim.sheils@nih.gov
* **Kyle Spencer** - kyle.spencer@nih.gov
* **Cole Tindall ** - cole.tindall@nih.gov

## Previous Authors
* **Bofei Zhang** - [Bofei5675](https://github.com/Bofei5675)
* **Shunchao Wang** - shunchao.wang@osumc.edu
* **Rohith Vanam** - rohith.vanam@osumc.edu
* ** Jorge Neyra** - jorge.neyra@nih.gov
