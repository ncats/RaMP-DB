#' The function write a data.frame to a csv files.
#' 
#' @param df a data frame returned by functions that queires database.
#' @param f.name a string that represents output file name.
#' @param write.from a string that specifies what type of tables you want to write
#' @param perc_analyte_overlap parameter for pathway enrichment analysis
#' @param min_pathway_tocluster parameter for pathway enrichment analysis
#' @param perc_pathway_overlap parameter for pathway enrichment analysis
#' @export
write_to_csv <- function(df,f.name,write.from,
                         perc_analyte_overlap,min_pathway_tocluster,
                         perc_pathway_overlap
                         ){
  if(!(write.from %in% c('tab3_s1','tab3_s2','tab3_s2_summary_search',
                         'tab3_s2_summary_fisherResults')))
    stop('Wrong arguments for write.from')
  if(write.from == 'tab3_s1'){
    df = df[,c("pathwayName","pathwaysourceId","pathwaysource")]
    write.csv(df,file = f.name,row.names = F)
  } else if(write.from == "tab3_s2"){
    df= df[,c("pathwayName","pathwaysourceId",
              "pathwaysource","commonName")]
    write.csv(df,file = f.name,row.names = F)
  } else if(write.from == 'tab3_s2_summary_search'){
    df <- as.data.frame(table(df$commonName))
    colnames(df) <- c("Query","Num_Pathways")
    write.csv(df,file = f.name,row.names = F)
  } else if(write.from == 'tab3_s2_summary_fisherResults'){
    # if data.frame is from fisherTestResultsSignificant() in shiny app
    # so the data.frame is from output of FilterFishersResults()
    data <- df
    cluster_output<-RaMP::findCluster(fishers_df=data,perc_analyte_overlap=as.numeric(perc_analyte_overlap),
                           min_pathway_tocluster=as.numeric(min_pathway_tocluster),
                           perc_pathway_overlap=as.numeric(perc_pathway_overlap))
    cluster_list<-out$cluster_list
    if(length(unique(cluster_list))>1){
      print(paste0(length(cluster_list)," clusters found"))
    }else{
      print("Clustering failed")
    }
    # cluster_list is all clusters that found by findCluster()$cluster_list
    # cluster_output is all clusters from findCluster()
    
    out <- df
    rampout <- out$fishresults
    if(!is.null(rampOut)) {
      if(out$analyte_type=="both"){
        rampOut<-rampOut[,c("pathwayName","Pval.Metab","Num_In_Path.Metab","Total_In_Path.Metab",
                            "Pval.Gene", "Num_In_Path.Gene","Total_In_Path.Gene", "Pval_combined",
                            "Pval_combined_FDR","Pval_combined_Holm","pathwaysourceId","pathwaysource",
                            "cluster_assignment","rampids")]
        
        colnames(rampOut)<-c("Pathway Name", "Raw Fisher's P Value (Metabolites)","User Metabolites in Pathway",
                             "Total Metabolites in Pathway","Raw Fisher's P Value (Genes)","User Genes in Pathway",
                             "Total Genes in Pathway","Raw Fisher's P Value (Combined)","FDR Adjusted P Value (Combined)",
                             "Holm Adjusted P Value (Combined)","Source ID","Source DB","In Cluster","rampids")
        rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value (Combined)"]),]
      }else{
        results_fisher<-rampOut[,c("pathwayName","Pval","Pval_FDR","Pval_Holm","pathwaysourceId","pathwaysource",
                                   "Num_In_Path","Total_In_Path","cluster_assignment","rampids")]
        colnames(rampOut)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
                             "Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway",
                             "In Cluster","rampids")
        rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value"]),]
      }
      write.csv(rampOut,f.name,row.names = FALSE)
    }else{
      write.csv(c("No significant results"),f.name,row.names = FALSE)
    }
  }
}