#' Query to retrieve database ID prefixes for analyte types
#' @param analyteType value to indicate the desired analyte type. Value one of 'gene' or 'metabolite'
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @return Returns list of database ID prefixes for selected 'gene' or 'metabolite'
#' @examples
#' \dontrun{
#'   metabprefixes <- getPrefixesFromAnalytes( "metabolite", db = rampDB )
#' }
#'
#' @export  getPrefixesFromAnalytes
getPrefixesFromAnalytes<-function(analyteType="gene", db = RaMP()) {
  if (analyteType=="gene"){
    df1<- db@api$getGeneIDTypes()
    df1 <- data.frame(analyteType="Genes/Proteins", idTypes=paste(df1$IDtype,collapse=", "))
    }
  else if (analyteType=="metabolite"){
    df1 <- db@api$getMetaboliteIDTypes()
    df1 <- data.frame(analyteType="Metabolites", idTypes=paste(df1$IDtype,collapse=", "))
    }
  else{
    print("analyteType must be 'gene' or 'metabolite'")
  }
  return(df1)
}


#' Returns class data sources for metabolites
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @return Returns list of data sources for metabolites
getMetabClassDataSources<-function(db = RaMP()){
  return(db@api$getMetaboliteClassSources())
}


#' Returns class categories for metabolites
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @return Returns metabolite class types
#' @export getMetabClassTypes
getMetabClassTypes<-function(db = RaMP()){
  return(db@api$getMetaboliteClassTypes())
}

#' Returns chemical classes for classType
#' @param classType one of the metab class types as returned by getMetabClassTypes() function; if null
#' will return a list of available classes for each class type
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @return Returns metabolite classes for classTypes
#' @export getMetabChemClass
getMetabChemClass <- function( classType= 'ClassyFire_super_class', db = RaMP() ) {
  if (!is.null(classType)) {
    res <- db@api$getMetaboliteClassesForType(classType=classType)
    res <- split(res$class_name,res$class_level_name)
  }

  else if (is.null(classType)){
    res <- db@api$getAllMetaboliteClasses()
    res <- split(res$class_name,res$class_level_name)
  }

  if (length(res)==0){
    message(paste0(classType, ' Not Defined. Available Class Types:'))
    res<-getMetabClassTypes()

  }

  return(res)
}


#' Returns list of available ontologies from database
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @return Returns ontologies listed in database
#' @export getOntologies

getOntologies <- function(db = RaMP()) {
  return (db@api$getOntologies())
}


