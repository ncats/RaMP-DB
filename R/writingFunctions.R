#' Function that writes output from getPathwayFromAnalyte() to a CSV file
#'
#' @param mypathways data frame returned by function getPathwayFromAnalyte()
#' @param outputfile name of output file
#' @export
writePathwaysToCSV <- function(mypathways=NULL,outputfile=NULL) {
	if(is.null(mypathways) || is.null(outputfile)) {
		stop("Be sure to specify the output of the function getPathwayFromAnalyte() and an output file")
	}
	if(!all(c("pathwayName","pathwaysourceId",
		"pathwaysource","commonName") %in% colnames(mypathways))) {
		stop("Make sure that your input data is the output of the function getPathwayFromAnalyte()")
	}
	out= mypathways[,c("pathwayName","pathwaysourceId",
       	      "pathwaysource","commonName")]
    	utils::write.csv(out,file = outputfile,row.names = F)
}


#' Function that writes Fishers Test results, after clustering to a CSV file
#' 
#' @param fishResults a data frame returned by function runCombinedFisherTest()
#' @param outputfile name of output file
#' @param rampid whether or not to include rampId (default is FALSE)
#' @export
write_FishersResults <- function(fishResults=NULL,outputfile=NULL,rampid=FALSE){
        if(is.null(fishResults)) {
                stop("Be sure to specify the output of the function findCluster()")
        }
	clusters <- fishResults$cluster_list
	if(is.null(clusters)) {
	
		out <- fishResults$fishresults
		mycols <- setdiff(colnames(out),c("pathwayRampId","pathwayName"))
		mycols <- c("pathwayName",mycols)
		utils::write.csv(out[,mycols],file=outputfile,row.names=F)
	}
	else {
		cluster_list<-fishResults$cluster_list
		out <- fishResults
		rampOut=out$fishresults
		if(!is.null(rampOut)) {
			if(out$analyte_type=="both"){
				rampOut<-rampOut[,c("pathwayName","Pval.Metab","Num_In_Path.Metab","Total_In_Path.Metab",
                            "Pval.Gene", "Num_In_Path.Gene","Total_In_Path.Gene", "Pval_combined",
                            "Pval_combined_FDR","Pval_combined_Holm","pathwaysourceId","pathwaysource",
                            "cluster_assignment","rampids")]
        
        		    colnames(rampOut)<-c("Pathway Name", "Raw Fisher's P Value (Metabolites)","User Metabolites in Pathway",
                             "Total Metabolites in Pathway","Raw Fisher's P Value (Genes)","User Genes in Pathway",
                             "Total Genes in Pathway","Raw Fisher's P Value (Combined)","FDR Adjusted P Value (Combined)",
                             "Holm Adjusted P Value (Combined)","Source ID","Source DB","Cluster","rampids")
        		    rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value (Combined)"]),]
      			} else{
        		  results_fisher<-rampOut[,c("pathwayName","Pval","Pval_FDR","Pval_Holm","pathwaysourceId","pathwaysource",
                                   "Num_In_Path","Total_In_Path","cluster_assignment","rampids")]
        colnames(rampOut)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
                             "Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway",
                             "Cluster","rampids")
        		  rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value"]),]
      			}
#	      utils::write.csv(rampOut,outputfile,row.names = FALSE)
    		}else{
      			rampOut <- "No significant results"		
#	utils::write.csv(c("No significant results"),outputfile,row.names = FALSE)
    		}
	}

	if (!rampid) {
		rampOut=rampOut[,-ncol(rampOut)]
	}
	if(is.null(outputfile)) {
		return(rampOut[order(rampOut[,"Cluster"]),])
	} else {
		utils::write.csv(rampOut[order(rampOut[,"Cluster"]),],outputfile,row.names = FALSE)
	}

}
