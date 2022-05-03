[![Build Status](https://api.travis-ci.com/ncats/RaMP-DB.svg?branch=custom_bg_fix2)](https://travis-ci.com/github/ncats/RaMP-DB)

# New!  RaMP 2.0!

RaMP 2.0 is now released and includes an updated backend database with expanded annotations for >150,000 metabolites and ~14,000 genes/proteins.  Annotations include biological pathways, chemical classes and structures (for metabolites only), ontologies (metabolites only), and enzyme-metabolite relationships based on chemical reactions. Annotations are drawn from HMDB, KEGG (through HMDB), Lipid-MAPS, WikiPathways, Reactome, and CheBI. 

This R package includes functions that allow users to interface with this up-do-date and comprehensive resource.  Functionalities include 1) simple and batch queries for pathways, ontologies, chemical annotations, and reaction-level gene-metabolite relationships; 2) pathway and chemical enrichment analyses.

The code used to build the backend RaMP database is freely available at https://github.com/ncats/RaMP-Backend.

Please [click here to view our preprint](https://www.biorxiv.org/content/10.1101/2022.01.19.476987v1).

# Web Interface
Our new revamped web interface can be found at https://rampdb.nih.gov/.  The code is publicly available at https://github.com/ncats/RaMP-Client/.

# APIs
API access is now available [here](https://ramp-api-alpha.ncats.io/__docs__/).

# Why RaMP (Relational Database of Metabolomic Pathways)

The purpose of RaMP is to provide a publicly available database that integrates metabolite and gene/protein biological, chemical and other from multiple sources. The database is relational and it can be directly downloaded as a MySQL dump (under the folder inst/extdata) for integration into any tools.  Please note that this project is in continuous development and we appreciated any feedback. 

## Contact Info:
For any questions or feedback, please send us a note at [NCATSRaMP@mail.nih.gov](NCATSRaMP@mail.nih.gov).  

If you find a bug, please submit an issue through this GitHub repo.   

## Basic Features:
The R packages and associated app perform  the following queries:

	1. Retrieve analytes (genes, proteins, metabolites) given pathway(s) as input.
	2. Retrieve pathway annotations given analytes as input.
	3. Retrieve chemical annotations/structures given metabolites as input.
	4. Retrieve analytes involved in the same reaction (e.g. enzymes catalyzing reactions involving input metabolites)
	5. Retrieve ontologies (e.g. biospecimen location, disease, etc.) given input meteabolites.

The following analyses are also supported:

	1. Multi-omic pathway enrichment analysis
	2. Chemical enrichment analyses

Last date of dump file update: 12/13/2021

## Vignette
Detailed instructions for installing RaMP locally are below.  We've also put together a vignette to get you started on the analyses.  Click here for [vignette](https://ncats.github.io/RaMP-DB/RaMP_Vignette.html).


## Citation
If you use RaMP, please cite the following work:

Zhang, B., et al., RaMP: A Comprehensive Relational Database of Metabolomics Pathways for Pathway Enrichment Analysis of Genes and Metabolites. Metabolites, 2018. 8(1).

PMID: 29470400; PMCID: PMC5876005; DOI: 10.3390/metabo8010016

To access, [click here](https://www.mdpi.com/2218-1989/8/1/16)

## Installation Instructions
In order to use this R package locally, you will need the following:
* The R code under this repo
* Our latest MySQL database dump file. **Download [here](https://figshare.com/ndownloader/files/34941486).**

If you would like to know how to build RaMP database from scratch, please check another GitHub site at [RaMP-BackEnd](https://github.com/ncats/RaMP-BackEnd)

### MySQL set-up
RaMP requires that MySQL and the RaMP database be set up on the machine that you will be running the R package from.
To download MySQL, you can go to the [MySQL Downloads page](https://www.mysql.com/downloads/)

When installing, you will be prompted to create a password for the user "root", or it will create one automatically for you.  **Importantly, remember your MySQL password!**  You will need to get into mysql and to pass it as an argument to the RaMP R shiny web application.

If you want to reset your password , you can go to [MySQL References 5.7 - How to reset root password ] (https://dev.mysql.com/doc/refman/5.7/en/resetting-permissions.html)

Please note that you will need administrator privileges for this step..

If you are using a Mac, we recommend using brew to install MySQL.  Here's a good tutorial: https://www.novicedev.com/blog/how-install-mysql-macos-homebrew.

### Creating the database locally
Once your MySQL environment is in place, creating the RaMP database locally is trivial.
First, launch MySQL and create the database:
```
> mysql -u root -p
mysql> create database ramp;
mysql> exit;
```

Here, we are naming the database "ramp" but you can use any name you'd like.  It is worth noting though that the R package assumes that the name of the database is "ramp" by default.  So if you change the name, remember to pass that name as arguments in the R package functions.

Second, download and unzip the latest RaMP database from the inst/extdata folder.

Third, populate the named database with the mysql dump file that was previously downloaded. ([MySQL dump file.](https://figshare.com/ndownloader/files/34941486))
Unzip the sql.gz file and run the command below to load the unzipped databased sql file. 

```
> mysql -u root -p ramp < /path/to/mysqldump/rampXXXX.sql  
```

You're done!

Your "ramp" database should contain the following 8 tables:
1. analyte
2. analyehasontology
3. analytehaspathway
4. analytesynonym
5. binary_source_matrix_view
5. catalyzed
6. chem_props
7. db_version
8. metabolite_class
9. ontology
10. pathway
11. source
12. version_info

If you want to explore this in MySQL, you can try:
```
mysql -u root
use ramp;
show tables;
select * from source limit 4; 
select * from source where commonName = "creatine riboside";
select distinct(HMDBOntologyType) from ontology;
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

# Set up your connection to the RaMP2.0 database:
pkg.globals <- setConnectionToRaMP(dbname="ramp",username="root",conpass="",host = "localhost")
```

Note that prior to using RaMP functions, users much establish required parameters to appropriately connect to your local database (if you are not using the web app).  This step is simplified by a single function call (last line in the above code snippet).

If the username is different then root, then specify the username in the "username" parameter.  Similarly, if the name of the database is different than "ramp2", then specify the "dbname" parameter.

### Important Notes

If you reinstall the latest version of the RaMP package, be sure to also install the latest version of the MySQL RaMP dump file.  

Also, when gene or metabolite ids are input for queries, IDs should be prepended with their database of origin, e.g. kegg:C02712, hmdb:HMDB04824, or CAS:2566-39-4. The list of metabolite or gene/protien IDs may be of mixed source. Remember to include the colon in the prefix. The id prefixes that are currently included in RaMP are: 

| Analyte Type | ID Prefix Types |
|--------------|-----------------|
| Metabolites | hmdb, pubchem, chebi, chemspider, kegg, CAS, LIPIDMAPS, swisslipids, lipidbank, wikidata, plantfa, kegg_glycan |
| Genes/Proteins | ensembl, entrez, gene_symbol, uniprot, hmdb, ncbiprotein, EN, wikidata, chebi |

To query the ID types supports in MySQL:

```
select distinct(IDtype) from source where geneOrCompound ="compound";
mysql> select distinct(IDtype) from source where geneOrCompound ="gene";
```

## Current Authors and Testers
* **John Braisted** - john.braisted@nih.gov
* **Tara Eicher** - tara.eicher@nih.gov
* **Ewy Mathé** - ewy.mathe@nih.gov
* **Andrew Patt** - andy.patt@nih.gov
* **Tim Sheils** - tim.sheils@nih.gov
* **Kyle Spencer** - kyle.spencer@nih.gov
* **Cole Tindall** - cole.tindall@nih.gov

## Previous Authors
* **Bofei Zhang** - [Bofei5675](https://github.com/Bofei5675)
* **Shunchao Wang** - 
* **Rohith Vanam** - 
* **Jorge Neyra** - [Jorgeso](https://github.com/jorgeso)
