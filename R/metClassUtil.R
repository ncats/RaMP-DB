#' Query to retrieve database ID prefixes for analyte types
#' @param analyteType value to indicate the desired analyte type. Value one of 'gene' or 'metabolite'
#' @return Returns list of database ID prefixes for selected 'gene' or 'metabolite'
#' @export  getPrefixesFromAnalytes
getPrefixesFromAnalytes<-function(analyteType="gene") {
  con <- connectToRaMP()
  if (analyteType=="gene"){
    query1 <- "select distinct(IDtype) from source where geneOrCompound ='gene';"
    df1<- RMariaDB::dbGetQuery(con,query1)
    df1 <- data.frame(analyteType="Genes/Proteins", idTypes=paste(df1$IDtype,collapse=", "))
    }
  else if (analyteType=="metabolite"){
    query2 <- "select distinct(IDtype) from source where geneOrCompound ='compound';"
    df1 <- RMariaDB::dbGetQuery(con,query2)
    df1 <- data.frame(analyteType="Metabolites", idTypes=paste(df1$IDtype,collapse=", "))
    }
  else{
    print("analyteType must be 'gene' or 'metabolite'")
  }
  RMariaDB::dbDisconnect(con)
  return(df1)
}


#' Returns class data sources for metabolites
#' @return Returns list of data sources for metabolites

getMetabClassDataSources<-function(){
  con <- connectToRaMP()
  query1<-"select distinct(source) from metabolite_class order by source asc"
  results<- RMariaDB::dbGetQuery(con,query1)
  return(results)
}


#' Returns class categories for metabolites
#' @return Returns metabolite class types
#' @export getMetabClassTypes
getMetabClassTypes<-function(){
  con <- connectToRaMP()
  query1<-"select distinct(class_level_name) from metabolite_class order by class_level_name asc"
  results<- RMariaDB::dbGetQuery(con,query1)
  return(results)
}

#' Returns chemical classes for classType
#' @param classType one of the metab class types as returned by getMetabClassTypes() function; if null
#' will return a list of available classes for each class type
#' @return Returns metabolite classes for classTypes
#' @export getMetabChemClass
getMetabChemClass<-function(classType= 'ClassyFire_super_class') {
  if (!is.null(classType)) {
    con <- connectToRaMP()
    query1<- paste0("select class_level_name, class_name from metabolite_class where class_level_name = '",classType,"' group by class_level_name, class_name")
    res<-RMariaDB::dbGetQuery(con,query1)
    res<-split(res$class_name,res$class_level_name)
  }

  else if (is.null(classType)){
    con <- connectToRaMP()
    query1<- "select class_level_name, class_name from metabolite_class group by class_level_name, class_name"
    res<-RMariaDB::dbGetQuery(con,query1)
    res<-split(res$class_name,res$class_level_name)
  }

  if (length(res)==0){
      message(paste0(classType, ' Not Defined. Available Class Types:'))
      res<-getMetabClassTypes()

  }
  return(res)
  }


  #' Returns list of available ontologies from database
  #' @return Returns ontologies listed in database
  #' @export getOntologies

getOntologies<-function(){
    con <- connectToRaMP()
    query1<- "select distinct(commonName) from ontology"
    OntoRes<-RMariaDB::dbGetQuery(con,query1)
    return(OntoRes)
}


