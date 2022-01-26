#' Function that query database to find ontology information based on
#' the given list of analytes
#' @param analytes a vector of analytes or a analytes delimited by new line character
#' @param NameOrIds specify the type of given data
#' @return dataframe that contains searched ontology from given analytes
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname = "ramp2", username = "root",
#' 	conpass = "", host = "localhost")
#' getOntoFromMeta("hmdb:HMDB0071437")
#' }
#' @export
getOntoFromMeta <- function(analytes, NameOrIds = "ids") {
  if (!(NameOrIds %in% c("ids", "name"))) {
    stop("Specifiy the type of given data to 'ids' or 'name'")
  }

  now <- proc.time()
  if (is.character(analytes)) {
    if (grepl("\n", analytes)[1]) {
      list_metabolite <- strsplit(analytes, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else if (grepl(",", analytes)[1]) {
      list_metabolite <- strsplit(analytes, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- analytes
    }
  } else if (is.data.frame(analytes)) {
    list_metabolite <- unlist(analytes)
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite, shQuote)
  list_metabolite <- paste(list_metabolite, collapse = ",")
  con <- connectToRaMP()
  if (NameOrIds == "ids") {
    sql <- paste0("select * from source where sourceId in (", list_metabolite, ");")
  } else if (NameOrIds == "name") {
    sql <- paste0("select * from source where rampId in (select * from (select rampId from analytesynonym where Synonym in (", list_metabolite, ")) as subquery);")
    cat(file = stderr(), "query sql in Package call with -- ", sql, "\n")
  }
  df <- RMariaDB::dbGetQuery(con, sql)
  # print(colnames(df))
  RMariaDB::dbDisconnect(con)
  if (nrow(df) == 0) {
    message("This source id
            does not exist in the source table")
    return(NULL)
  }
  rampid <- unique(df$rampId)
  rampid <- sapply(rampid, shQuote)
  rampid <- paste(rampid, collapse = ",")
  con <- connectToRaMP()

  sql <- paste0(
    "select * from analytehasontology where rampCompoundId in (",
    rampid, ");"
  )
  df2 <- RMariaDB::dbGetQuery(con, sql)
  if (nrow(df2) == 0) {
    message("No searching result because these metabolites are not linked to ontology")
    return(NULL)
  }
  RMariaDB::dbDisconnect(con)
  rampontoid <- unique(df2$rampOntologyId)
  rampontoid <- sapply(rampontoid, shQuote)
  rampontoid <- paste(rampontoid, collapse = ",")
  sql <- paste0(
    "select * from ontology where rampOntologyId in (",
    rampontoid, ");"
  )
  con <- connectToRaMP()
  df3 <- RMariaDB::dbGetQuery(con, sql)
  # print(colnames(df3))
  RMariaDB::dbDisconnect(con)
  mdf <- unique(merge(df3, df2, all.x = T))
  mdf <- unique(merge(mdf, df,
    all.x = T, by.x = "rampCompoundId",
    by.y = "rampId"
  ))
  colnames(mdf)[colnames(mdf) == "commonName.x"] <- "Ontology"
  colnames(mdf)[colnames(mdf) == "commonName.y"] <- "Metabolites"

  mdf <- mdf[c("Metabolites", "sourceId", "IDtype", "Ontology", "HMDBOntologyType")]

  # need to make unique list (JB, 12/1/21)
  mdf <- unique(mdf)

  # colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
  return(mdf)
}
#' function that query database to find analytes in given ontologies
#' @param ontology a vector of ontology or ontologies delimited by new line character
#'
#' @return dataframe that  contains searched analytes from given ontology
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname = "ramp2", username = "root",
#' 	conpass = "", host = "localhost")
#' getMetaFromOnto("Adiposome")
#' }
#' @importFrom rlang .data
#' @export
getMetaFromOnto <- function(ontology) {
  now <- proc.time()
  if(is.character(ontology)){
    if(grepl("\n",ontology)[1]){
      list_ontology <- strsplit(ontology,"\n")
      list_ontology <- unlist(list_ontology)
    } else if(grepl(",",ontology)[1]){
      list_ontology <- strsplit(ontology,",")
      list_ontology <- unlist(list_ontology)
    } else {
      list_ontology <- ontology
    }
  } else if(is.data.frame(ontology)){
    list_ontology <- unlist(ontology)
  }
  list_ontology <- unique(list_ontology)
  #list_ontology <- sapply(list_ontology,shQuote)

  #JCB Mod 12/20/21 - drop collapse at this point, until we determine matched ontologies
  # Comment out
  #list_ontology <- paste(list_ontology,collapse = ",")

  # Find the ontology in ramp2
  start = Sys.time()

  allontos <- getOntologies()
  matched_ontos <- unlist(lapply(list_ontology,
	function(x) grep(x, allontos$commonName)))
  # Create df which will be used later to create the output
  df <- allontos[matched_ontos, ]

  end <- Sys.time()
  print(paste0("Hey it took this long to get my onto id: ", end-start) ))
  start = Sys.time()

  # con <- connectToRaMP()
  # sql <- paste0('select * from ontology where commonName in (',
  #              list_ontology,');')

  # df <- RMariaDB::dbGetQuery(con,sql)
  # RMariaDB::dbDisconnect(con)
  # print(colnames(df))


  # Only proceed if df has anything returned
  if (length(matched_ontos) > 0) {
    # rampontoid <- paste(sapply(unique(df$rampOntologyId),shQuote),
    #              collapse = ',')

    rampontoid <- paste(sapply(unique(allontos$rampOntologyId[matched_ontos]), shQuote),
      collapse = ","
    )

    # Find RaMP IDs associated with ontology
    con <- connectToRaMP()
    sql <- paste0(
      "select * from analytehasontology where rampOntologyId in (",
      rampontoid, ");"
    )

    df2 <- RMariaDB::dbGetQuery(con, sql)
    RMariaDB::dbDisconnect(con)

    end = Sys.time()
    print(paste0("Hey it took this long to get * from analytehasontology: ", end-start) ))


    # connect ramp IDs with source table to get source IDs and names
    # print(colnames(df2))
    rampid <- paste(sapply(unique(df2$rampCompoundId), shQuote), collapse = ",")
    con <- connectToRaMP()
    sql <- paste0("select * from source where rampId in (", rampid, ");")
    df3 <- RMariaDB::dbGetQuery(con, sql)
    df3 <- unique(df3)
    # print(colnames(df3))
    RMariaDB::dbDisconnect(con)
    if (nrow(df3) == 0) {
      message("No searching result since this ramp id
	            does not exist in source table")
      return(NULL)
    }

    ## Merge ontology with nice metabolite names
    # print('Merging 1...')
    mdf <- unique(merge(df3, df2, by.x = "rampId", by.y = "rampCompoundId", all.x = T))
    # print('Merging 2...')
    mdf <- unique(merge(mdf, df,
      by.x = "rampOntologyId",
      by.y = "rampOntologyId",
      all.x = T
    ))
    # mdf <- mdf[c('sourceId','commonName.y','biofluidORcellular')]
    # colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
    colnames(mdf)[colnames(mdf) == "commonName.x"] <- "Metabolites"
    colnames(mdf)[colnames(mdf) == "commonName.y"] <- "Ontology"
      ## Put RaMP ids back in and concatenate by that. Use HMDB name > Reactome name > Wiki name
      mdf_final <- lapply(unique(mdf$rampId), function(x){
          temp_df<-mdf %>%
              dplyr::filter(.data$rampId == x)
          if (nrow(temp_df)==1){
              return(temp_df)
          }else{
              if ("hmdb" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "hmdb"))
              }else if ("LIPIDMAPS" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "LIPIDMAPS"))
              }else if ("kegg" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "kegg"))
              }else if ("wikidata" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "wikidata"))
              }else if ("pubchem" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "pubchem"))
              }else if ("chebi" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "chebi"))
              }else if ("CAS" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "CAS"))
              }else if ("chemspider" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "chemspider"))
              }else if ("lipidbank" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "lipidbank"))
              }else if ("swisslipids" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "swisslipids"))
              }else if ("plantfa" %in% temp_df$dataSource){
                  return(temp_df %>%
                         dplyr::filter(.data$IDtype == "plantfa"))
              }
          }
      })
      mdf_final <- do.call(rbind.data.frame, mdf_final)
      mdf_final <- mdf_final[c("Metabolites", "sourceId", "Ontology", "HMDBOntologyType")]
      mdf_final <- unique(mdf_final)
    # Collapse rows for each metabolite:
    # Metabolites<-c()
    # print("NO ERROR YET")
    # print(dim(mdf))
    # out=data.table::setDT(mdf)[,lapply(data.table::.SD, function(x) paste0(unique(x),collapse=", ")),
                                        # 	by=Metabolites]
    return(mdf_final)
  } else {
    warning("Ontology not found in RaMP.  Run the getOntologies() function to see available ontologies.")
    return(NA)
  }
}



#' function that query database to find analytes in given ontologies
#' @param ontology a vector of ontology or ontologies delimited by new line character
#'
#' @return dataframe that  contains searched analytes from given ontology
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname = "ramp2", username = "root",
#' 	conpass = "", host = "localhost")
#' getMetaFromOnto("Adiposome")
#' }
#' @importFrom rlang .data
#' @export
getMetaFromOntoV2 <- function(ontology) {
  now <- proc.time()
  if(is.character(ontology)){
    if(grepl("\n",ontology)[1]){
      list_ontology <- strsplit(ontology,"\n")
      list_ontology <- unlist(list_ontology)
    } else if(grepl(",",ontology)[1]){
      list_ontology <- strsplit(ontology,",")
      list_ontology <- unlist(list_ontology)
    } else {
      list_ontology <- ontology
    }
  } else if(is.data.frame(ontology)){
    list_ontology <- unlist(ontology)
  }
  list_ontology <- unique(list_ontology)
  #list_ontology <- sapply(list_ontology,shQuote)

  #JCB Mod 12/20/21 - drop collapse at this point, until we determine matched ontologies
  # Comment out
  #list_ontology <- paste(list_ontology,collapse = ",")

  # Find the ontology in ramp2
  start = Sys.time()

  allontos <- getOntologies()
  matched_ontos <- unlist(lapply(list_ontology,
                                 function(x) grep(x, allontos$commonName)))
  # Create df which will be used later to create the output
  df <- allontos[matched_ontos, ]

  end <- Sys.time()
  print(paste0("Hey it took this long to get my onto id: ", end-start) )
start = Sys.time()

# con <- connectToRaMP()
# sql <- paste0('select * from ontology where commonName in (',
#              list_ontology,');')

# df <- RMariaDB::dbGetQuery(con,sql)
# RMariaDB::dbDisconnect(con)
# print(colnames(df))


# Only proceed if df has anything returned
if (length(matched_ontos) > 0) {
  # rampontoid <- paste(sapply(unique(df$rampOntologyId),shQuote),
  #              collapse = ',')

  rampontoid <- paste(sapply(unique(allontos$rampOntologyId[matched_ontos]), shQuote),
                      collapse = ","
  )

  # Find RaMP IDs associated with ontology
  con <- connectToRaMP()
  # sql <- paste0(
  #   "select * from analytehasontology where rampOntologyId in (",
  #   rampontoid, ");"
  # )


  sql <- paste0("select rampId,
  group_concat(distinct s.sourceId order by s.sourceId separator ', ') as source_ids,
  group_concat(distinct s.commonName order by s.commonName separator ', ') as common_names, o.commonName
  from source s, analytehasontology ao, ontology o where ao.rampOntologyId in (
  select distinct rampOntologyId from ontology where commonName in (",list_ontology,")
  )and o.rampOntologyId = ao.rampOntologyId
  and s.rampId = ao.rampCompoundId
  group by o.commonName, s.rampId")


  df2 <- RMariaDB::dbGetQuery(con, sql)
  RMariaDB::dbDisconnect(con)

  end = Sys.time()
  print(paste0("Hey it took this long to get * from analytehasontology: ", end-start) )


# connect ramp IDs with source table to get source IDs and names
# print(colnames(df2))
rampid <- paste(sapply(unique(df2$rampCompoundId), shQuote), collapse = ",")
con <- connectToRaMP()
sql <- paste0("select * from source where rampId in (", rampid, ");")
df3 <- RMariaDB::dbGetQuery(con, sql)
df3 <- unique(df3)
# print(colnames(df3))
RMariaDB::dbDisconnect(con)
if (nrow(df3) == 0) {
  message("No searching result since this ramp id
	            does not exist in source table")
  return(NULL)
}

## Merge ontology with nice metabolite names
# print('Merging 1...')
mdf <- unique(merge(df3, df2, by.x = "rampId", by.y = "rampCompoundId", all.x = T))
# print('Merging 2...')
mdf <- unique(merge(mdf, df,
                    by.x = "rampOntologyId",
                    by.y = "rampOntologyId",
                    all.x = T
))
# mdf <- mdf[c('sourceId','commonName.y','biofluidORcellular')]
# colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
colnames(mdf)[colnames(mdf) == "commonName.x"] <- "Metabolites"
colnames(mdf)[colnames(mdf) == "commonName.y"] <- "Ontology"
## Put RaMP ids back in and concatenate by that. Use HMDB name > Reactome name > Wiki name
mdf_final <- lapply(unique(mdf$rampId), function(x){
  temp_df<-mdf %>%
    dplyr::filter(.data$rampId == x)
  if (nrow(temp_df)==1){
    return(temp_df)
  }else{
    if ("hmdb" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "hmdb"))
    }else if ("LIPIDMAPS" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "LIPIDMAPS"))
    }else if ("kegg" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "kegg"))
    }else if ("wikidata" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "wikidata"))
    }else if ("pubchem" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "pubchem"))
    }else if ("chebi" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "chebi"))
    }else if ("CAS" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "CAS"))
    }else if ("chemspider" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "chemspider"))
    }else if ("lipidbank" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "lipidbank"))
    }else if ("swisslipids" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "swisslipids"))
    }else if ("plantfa" %in% temp_df$dataSource){
      return(temp_df %>%
               dplyr::filter(.data$IDtype == "plantfa"))
    }
  }
})
mdf_final <- do.call(rbind.data.frame, mdf_final)
mdf_final <- mdf_final[c("Metabolites", "sourceId", "Ontology", "HMDBOntologyType")]
mdf_final <- unique(mdf_final)
# Collapse rows for each metabolite:
# Metabolites<-c()
# print("NO ERROR YET")
# print(dim(mdf))
# out=data.table::setDT(mdf)[,lapply(data.table::.SD, function(x) paste0(unique(x),collapse=", ")),
# 	by=Metabolites]
return(mdf_final)
} else {
  warning("Ontology not found in RaMP.  Run the getOntologies() function to see available ontologies.")
  return(NA)
}



