#' Retrieves analytes that involved in same reaction as input metabolite
#'
#' @param analytes a vector of analytes that need to be searched
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @param NameOrIds whether input is "names" or "ids" (default is "ids")
#' @return a dataframe containing query results. If the input is a metabolite, the function will output
#' gene transcript common names and source IDs that are known to catalyze
#' reactions in the same pathway as that metabolite. Conversely, if the input
#' is a gene, the function will return the common name and source id of metabolites
#' known to be catalyzed directly or indirectly by that gene.
#'
#' @examples
#' \dontrun{
#' rampFastCata(analytes="creatine",conpass="mypassword",NameOrIds="names")
#' }
#' @export
rampFastCata <- function(analytes=NULL,conpass=NULL,
                         dbname="ramp",username="root",
                         host = "localhost",
                         NameOrIds="ids") {
  
  if(is.null(analytes))
    stop("Please provide input analytes")
  if (!(NameOrIds %in% c('names','ids')))
    stop('Please specify search by "names" or "ids"')
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }
  
  now <- proc.time()
  if(is.character(analytes)){
    if(grepl("\n",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- analytes
    }
  } else if(is.data.frame(analytes)){
    list_metabolite <- unlist(analytes)
  } else {stop("The input 'analytes' is not a recognized format. Please check input.")}
  
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  
  #  print(list_metabolite)
  
  # Retrieve RaMP analyte ids 
  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                        password = conpass,
                        dbname = dbname,
                        host = host)
  if (NameOrIds == 'names'){
    #    query1 <- paste0("select Synonym as analyte1,rampId,geneOrCompound as type1 from analytesynonym where Synonym in (",list_metabolite,");")
    query1 <- paste0("select rampId,geneOrCompound as type1,Synonym as InputAnalyte from analytesynonym where Synonym in (",list_metabolite,");")
  } else if (NameOrIds == 'ids'){
    #    query1 <- paste0('select rampId,geneOrCompound as type1,commonName as InputMetabolite from analytesynonym where rampId in (select rampId from source where sourceId in (',list_metabolite,'));')
    query1 <- paste0('select rampId,geneOrCompound as type1,commonName as InputMetabolite from source where sourceId in (',list_metabolite,');')
  }
  
  # Retrieves Name, RaMPID and type (gene or compound) for input
  df1<- DBI::dbGetQuery(con,query1)
  DBI::dbDisconnect(con)
  
  #print(df1$rampId)
  df_c <- df_g <- NULL
  mdf_c <- mdf_g <- NULL
  # Process metabolite ids
  mdf_cfin2 <- mdf_gfin2 <- c()
  if(length(grep("RAMP_C",df1$rampId))!=0){
    df_c <- df1[grep("RAMP_C",df1$rampId),]
    print("Get Compound ...")
    c_id <- unique(df_c$rampId)
    if(length(c_id) == 0){
      message("Input metabolites do not have catalyzation information")
      mdf_cfin2 <- NULL #return(NULL)
    } else {
      print(length(c_id))
      c_id <- sapply(c_id,shQuote)
      c_id <- paste(c_id,collapse = ",")
      
      # Retrieve rampid of genes that are in same reaction
      query_c <- paste0("select rampCompoundId as rampId,rampGeneId as rampId2 from catalyzed where rampCompoundId in (",c_id,");")
      print("Geting gene Id from Compound Id ...")
      con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                            password = conpass,
                            dbname = dbname,
                            host = host)
      df_c2 <- DBI::dbGetQuery(con,query_c)
      DBI::dbDisconnect(con) 
      if(nrow(df_c2) == 0){
        message("No genes found in same reaction as input metabolite")
        mdf_cfin2 <- NULL 
      } else {
        print("Getting names from gene Id ...")
        analyte2_list <- unique(df_c2$rampId2)
        analyte2_list <- sapply(analyte2_list,shQuote)
        analyte2_list <- paste(analyte2_list,collapse = ",")
        # Get names for metabolite ids
        query2 <- paste0("select * from source 
             		where rampId in (",analyte2_list,");")
        con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                              password = conpass,
                              dbname = dbname,
                              host = host)
        df_c3 <- DBI::dbGetQuery(con,query2)
        DBI::dbDisconnect(con)
        
        if(nrow(df_c3) == 0){
          message("Cannot retrieve names for those metabolites")
          mdf_cfin2 <- NULL
        } else {
          # Now merge it all:
          mdc_c <- merge(df_c,df_c2)
          colnames(df_c3)[which(colnames(df_c3)=="rampId")]="rampId2"
          mdf_cfin <- merge(mdc_c,df_c3)
          #print(colnames(mdf_cfin))
          if (NameOrIds == 'names'){
            mdf_cfin <- mdf_cfin[,c("InputAnalyte","sourceId","IDtype","commonName")]
          } else if (NameOrIds == 'ids'){
            mdf_cfin <- mdf_cfin[,c("InputMetabolite","sourceId","IDtype","commonName")]
          }
          colnames(mdf_cfin) <- c("Input_Metabolite","Gene_sourceId","Gene_IDtype",
                                  "Gene_CommonName")
          
          # Collapse source ids:
          mdf_cfin$temp <- paste(mdf_cfin[,"Input_Metabolite"],mdf_cfin[,"Gene_CommonName"])
          tempout <- data.frame(Input_Metabolite=NA,Gene_CommonName=NA,Gene_sourceIds=NA)
          mdf_cfin2=c()
          for (i in unique(mdf_cfin$temp)) {
            temp <- mdf_cfin[which(mdf_cfin$temp==i),]
            tempout$Input_Metabolite=temp[1,"Input_Metabolite"]
            tempout$Gene_sourceIds <- 
              paste(paste(temp$Gene_IDtype,temp$Gene_sourceId,sep=": "),
                    collapse="; ")
            tempout$Gene_CommonName=temp[1,"Gene_CommonName"]
            mdf_cfin2 <- rbind(mdf_cfin2,tempout)
          }
          colnames(mdf_cfin2) <- c("Input_Analyte","Input_CatalyzedBy_CommonName",
                                   "Input_CatalyzedBy_SourceIds")
        } # end else couldn't retrieve names for metabolites
      } # end else couldn't find metabolite ids
    } # no catalyzation information
  } # end if compound ids found
  
  # Do analagous for genes 
  if(length(grep("RAMP_G",df1$rampId))!=0){
    print("Also find gene inside")
    df_g <- df1[grep("RAMP_G",df1$rampId),]
    print("Get gene ...")
    g_id <- df_g$rampId
    g_id <- sapply(unique(g_id),shQuote)
    g_id <- paste(g_id,collapse = ",")
    
    if(length(g_id) == 0){
      message("No IDs found for input genes")
      mdf_gfin2 <- NULL #return(NULL)
    } else {
      # Get rampID for genes and catalyzed metabolites
      query_g <- paste0("select * from catalyzed where rampGeneId in (",g_id,");")
      
      con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                            password = conpass,
                            dbname = dbname,
                            host = host)
      df_g2 <- DBI::dbGetQuery(con,query_g)
      DBI::dbDisconnect(con)
      if(nrow(df_g2) == 0){
        message("Could not find metabolites in same reaction as input genes")
        mdf_gfin2=c()
      } else {
        analyte2_list <- df_g2$rampCompoundId
        analyte2_list <- sapply(analyte2_list,shQuote)
        analyte2_list <- paste(analyte2_list,collapse = ",")
        
        # Get names for metabolite IDs
        query2 <- paste0("select * from source where rampId in (",analyte2_list,");")
        con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                              password = conpass,
                              dbname = dbname)
        df_g3 <-DBI::dbGetQuery(con,query2)
        DBI::dbDisconnect(con)
        if(nrow(df_g3) == 0){
          message("Cannot retrieve names for those genes")
          mdf_gfin2 <- NULL
        } else {
          # Now merge it all:
          mdc_g <- merge(df_g,df_g2)
          colnames(df_g3)[which(colnames(df_g3)=="rampId")]="rampId2"
          mdf_gfin <- merge(mdc_g,df_g3)
          colnames(mdf_gfin)[colnames(mdf_gfin) == 'InputMetabolite'] = 'InputAnalyte'
          print(colnames(mdf_gfin))
          mdf_gfin <- mdf_gfin[,c("InputAnalyte","sourceId","IDtype","commonName")]
          colnames(mdf_gfin) <- c("Input_Gene","Gene_sourceId","Gene_IDtype",
                                  "Gene_CommonName")
          mdf_gfin <- mdf_gfin[!duplicated(mdf_gfin),]
          
          # Collapse source ids:
          mdf_gfin$temp <- paste(mdf_gfin[,"Input_Gene"],mdf_gfin[,"Gene_CommonName"])
          tempout <- data.frame(Input_Analyte=NA,Input_CatalyzedBy_CommonName=NA, 
                                Input_CatalyzedBy_SourceIds=NA)
          mdf_gfin2=c()
          for (i in unique(mdf_gfin$temp)) {
            temp <- mdf_gfin[which(mdf_gfin$temp==i),]
            tempout$Input_Analyte=temp[1,"Input_Gene"]
            tempout$Input_CatalyzedBy_CommonName=temp[1,"Gene_CommonName"]
            tempout$Input_CatalyzedBy_SourceIds <-
              paste(paste(temp$Gene_IDtype,temp$Gene_sourceId,sep=": "),collapse="; ")
            mdf_gfin2 <- rbind(mdf_gfin2,tempout)
          }
        } # else couldn't retrieve names for those genes
      } # couldn't find metabolites catalyzed by input genes
    } # end else couldn't find ids for input gene names
  } # end gene
  mdf <- rbind(mdf_cfin2,mdf_gfin2)
  print("Done ...")
  print("timing ...")
  print(proc.time() - now)
  #  if(!is.null(mdf)) {
  #  	colnames(mdf) <- c("Input_Analyte","Input_CatalyzedBy_CommonName",
  #		"Input_CatalyzedBy_SourceIds")
  #   }
  return(mdf)
}

#' Generate dataframe from given files for shiny app input list of metabolites
#'
#' @param infile a file object given from files
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#'
#' @return a dataframe either from multiple csv file
rampFileOfAnalytes_tab4 <- function(infile,conpass=NULL,
	dbname="ramp",username="root",host = "localhost"){
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }

  name <- infile[[1,'name']]
  summary <- data.frame(pathway  = character(0),id = character(0),
                        source = character(0),metabolite = character(0))
  rampOut <- list()
  for (i in 1:length(infile[,1])){
    rampOut <- readLines(infile[[i,'datapath']])
    summary <- rampFastCata(analytes=rampOut,conpass=conpass,
                            host = host)
  }
  return(summary)
}

#' Plots a network based on gene-metabolite relationships
#' @importFrom magrittr %>%
#'
#' @param catalyzedf a data.frame output by rampFastCata() that contains analytes that are in the same reaction
#' @return  An interactive HTML plot that allows the user to pan/zoom into regions of interest. User genes/metabolites are highlighted in blue, whereas analytes found by the function are orange.
#' @examples
#' \dontrun{
#' catalyzedf <- rampFastCata(analytes="creatine",conpass="mypassword",NamesOrIds="names")
#' plotCataNetwork(catalyzedf)
#' }
#' @export
plotCataNetwork <- function(catalyzedf = NULL) {

        if(is.null(catalyzedf) ||
        (length(intersect(c("Input_Analyte","Input_CatalyzedBy_CommonName",
                "Input_CatalyzedBy_SourceIds"),colnames(catalyzedf)))!=3)) {
                stop("Please make sure that the input is the resulting data.frame returned by the rampFastCata() function")
        }

        toplot = catalyzedf[,c("Input_Analyte","Input_CatalyzedBy_CommonName")]
        colnames(toplot)<-c("from","to")

        #colopts <- brewer.pal(12,"Set3")
        mycol=rep("black")
        mysize<-rep(8,nrow(toplot))
        mynames<-rep(NA,nrow(toplot))

        myedges <- cbind(toplot,mycol,mysize,mynames)

        # Now set nodes
        mynodes=c(unique(toplot$from),unique(toplot$to))
        mycol=c(rep("blue",length(unique(toplot$from))),
                rep("orange",length(unique(toplot$to))))
        mysize <- rep(8,length(mynodes))
        mynames <- mynodes

        mynodes=data.frame(color=mycol,size=mysize,id=mynames,label=mynames)

        # Now plot
        visNetwork::visNetwork(mynodes, myedges, width = "100%",height="1000px") %>%
		visNetwork::visInteraction(dragNodes = FALSE,
                 dragView = TRUE,hideEdgesOnDrag=TRUE,hideNodesOnDrag=TRUE,
                 navigationButtons=TRUE,zoomView = TRUE) %>%
  		visNetwork::visLayout(randomSeed = 123) %>%
		visNetwork::visPhysics(
		  barnesHut = list(
		    gravitationalConstant = -100,
		    centralGravity = 0,
		    springConstant = 0
		  ),
		  stabilization = TRUE)
        #return(NULL) #return(list(nodes=mynodes,edges=myedges))
}
