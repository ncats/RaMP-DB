# RaMP - Relational Database of Metabolomic Pathways

The purpose of RaMP is to provide a publicly available database that integrates metabolite and gene biological pathways from multiple sources. So far, we have integrated information from HMDB, KEGG Reactome, and WikiPathways. The relational structure of RaMP enables complex and batch queries.  To facilitate its usage, we have created an R package that includes a user-friendly R Shiny web application.  Please note that this project is in continuous development and we certainly appreciate your feedback.  

Also note that we are working on a server version of RaMP so that users do not have to install anything on their local machines.  Stay tuned!

## Basic Features:
The app performs some complex queries (e.g. retrieve allmetabolites and/or genes that belong to a user input pathway or list of pathways).  It also performs pathway enrichment analysis given a list of metabolites and/or genes.

Here are some screenshots to give you an idea.

Front page:

<img src="img/Picture1.png" alt = "FrontPage" width="400"/>

Search WorkFlow:<br/>

<img src="img/Picture2.png" alt = "FrontPage" width="400"/>
<img src="img/Picture3.png" alt = "FrontPage" width="400"/>
<img src="img/Picture4.png" alt = "FrontPage" width="400"/>

## Installation Instructions
In order to use the web application, you will need the following:
* The R code under this repo
* The mysql dump file that contains the RaMP database (in the folder inst/extdata/)

### MySQL set-up
RaMP requires that MySQL and the RaMP database be set up on the machine that you will be running the app from.
To download MySQL, you can go to the [MySQL Downloads page](https://www.mysql.com/downloads/)

When installing, you will be prompted to create a password for the user "root", or it will create one automatically for you.  **Importantly, remember your MySQL password!**  You will need to get into mysql and to pass it as an argument to the RaMP R shiny web application.

If you want to reset your password , you can go to [MySQL References 5.7 - How to reset root password ] (https://dev.mysql.com/doc/refman/5.7/en/resetting-permissions.html)

Please note that you will need administrator priviledges for this step..

### Creating the database locally
Once your MySQL environment is in place, creating the RaMP database locally is trivial.
First, launch MySQL and create the database:
```
> mysql -u root -p
mysql> create database ramp;
mysql> exit;
```

Here, we're calling the database ramp but you can use any name you'd like.  Realize though that R package assumes that the name of the database is RaMP by default.  So if you change the name, remember to pass that name as arguments in the R package functions.

Second, populate the named database with the mysql dump file (which you can get from  inst/extdata/ramp.sql):
```
> mysql -u root -p ramp < ramp.sql
```

You're done!

Note that your "ramp" database should contain the following 8 tables:
1. analyte
1. analyehasontology
1. analytehaspathway
1. analytesynonym
1. catalyzed
1. ontology
1. pathway
1. source

If you want to explore this in MySQL, you can try:
```
mysql -u root -p
use ramp;
view tables;
select * from analytesynonym where synonym = "glucose";
```


### Install and load the RaMP package 
You can install this package directly from GitHub. In the R Console, type the following:
```R
install_github("Bofei5675/RaMP2")
# Load the package
library(RaMP)
```

Now, you're set to use the web application locally.  Just type:
```R
RaMP::runRaMPapp(conpass="mysql_password")
```

If the username is different then root, then specify the username in the "username" parater.  Similarly, if the name of the database is different than "ramp", then specify the "dbname" parameter.

## Reporting Issues
If you encounter any problems, or find installation problems or bugs, please start an issue on the Issues tab or email Ewy.Mathe@osumc.edu directly. Thanks!

## Authors
* **Bofei Zhang** - [Bofei5675](https://github.com/Bofei5675)
* **Ewy Mathé** - [Mathelab](https://github.com/MatheLab)


