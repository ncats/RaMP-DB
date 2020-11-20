library(methods)

#' @export FisherTest
FisherTest <- setRefClass("FisherTest",
    fields = list (
        db_connection = "list",
        fishresults = "list",
        analyte_type = "character"
    ),
    methods = list(
        initialize = function() {
            if (is.null(db_connection)) {
                stop("Please provide a db_connection")
            }
        },
        #' Do fisher test for only one pathway from search result
        #' clicked on highchart
        #' @param pathwaydf a data frame resulting from getPathwayFromAnalyte
        #' @param total_metabolites number of metabolites analyzed in the experiment (e.g. background) (default is 1000; set to 'NULL' to retrieve total number of metabolites that map to any pathway in RaMP). Assumption that analyte_type is "metabolite")
        #' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
        #' @param analyte_type "metabolites" or "genes" (default is "metabolites")
        #' @param conpass password for database access (string)
        #' @param dbname name of the mysql database (default is "ramp")
        #' @param username username for database access (default is "root")
        #' @param host host name for database access (default is "localhost")
        #' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway
        run_fisher_test = function(
            pathwaydf,
            total_metabolites=NULL,
            total_genes=20000,
            analyte_type="metabolites"
        ) {
            now <- proc.time()
            print("Fisher Testing ......")

            if(analyte_type=="metabolites") {
                total_analytes=total_metabolites
            } else if (analyte_type=="genes") {
                total_analytes=total_genes
            } else {
                stop("Please define the analyte_type variable as 'metabolites' or 'genes'")
            }

            print(paste("early in fisher... total analytes =", toString(total_analytes),sep=" "))
            print(paste("early in fisher... total mets=", toString(total_metabolites),sep=" "))

            contingencyTb <- matrix(0,nrow = 2,ncol = 2)

            colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
            rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")

            # Get pathway ids that contain the user analytes
            pid <- unique(pathwaydf$pathwayRampId);

            analyte_pathway = AnalytePathwayRepository(db_connection)

            # Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
            allids <- analyte_pathway$get_analytes_pathways()

            allids <- allids[!duplicated(allids),]

            if((analyte_type == "metabolites") && (is.null(total_metabolites))) {

                # JCB replaced these lines. Reducing to a source, extracting compound indices, then applying to the full set of rampIds
                # caused an error in the tally of analytes
                #
                # wiki_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="wiki"),"rampId"])]))
                # react_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="reactome"),"rampId"])]))
                # kegg_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="kegg"),"rampId"])]))

                # first extract source-specific ids, then select for compound ids from the source-specific ids
                sourceIds <- allids[which(allids$pathwaySource=="wiki"),"rampId"]
                wiki_totanalytes <- length(unique(sourceIds[grep("RAMP_C",sourceIds)]))

                sourceIds <- allids[which(allids$pathwaySource=="reactome"),"rampId"]
                react_totanalytes <- length(unique(sourceIds[grep("RAMP_C",sourceIds)]))

                sourceIds <- allids[which(allids$pathwaySource=="kegg"),"rampId"]
                kegg_totanalytes <- length(unique(sourceIds[grep("RAMP_C",sourceIds)]))

            } else {
                print("either not metabolites or total_metabolites is not null")
            }

            if(analyte_type=="genes") {

                # for now we're using a fixed population size for genes
                # this can be enhanced to take a list of all measured genes
                # or use a subset of genes having pathway annotations within each source
                wiki_totanalytes <- react_totanalytes <- kegg_totanalytes <- total_genes
            }

            input_RampIds <- analyte_pathway$get_analytes_pathways(pathway_ramp_ids=pid)

            if(is.null(input_RampIds)) {

                stop("Data doesn't exist")

            } else {

            # data frames for metabolites with pathawayRampID, Freq based  on Source(kegg, reactome, wiki)

                input_RampId_C <- input_RampIds[grep("RAMP_C", input_RampIds$rampId), ]
                unique_input_RampId_C <- unique(input_RampId_C[,c("rampId", "pathwayRampId")])
                unique_pathwayRampId_source <- unique(input_RampId_C[,c("pathwayRampId", "pathwaySource")])

                freq_unique_input_RampId_C <- as.data.frame(table(unique_input_RampId_C[,"pathwayRampId"]))

                names(freq_unique_input_RampId_C)[1] = 'pathwayRampId'
                merge_Pathwayfreq_source <- merge(freq_unique_input_RampId_C, unique_pathwayRampId_source, by="pathwayRampId")

                # subset metabolite data based on source -  kegg, reactome, wiki

                input_kegg_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "kegg")
                input_reactome_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "reactome")
                input_wiki_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "wiki")

            # data frames for Genes with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)

                input_RampId_G <- input_RampIds[grep("RAMP_G", input_RampIds$rampId), ]
                unique_input_RampId_G <- unique(input_RampId_G[,c("rampId", "pathwayRampId")])
                unique_pathwayG_source <- unique(input_RampId_G[,c("pathwayRampId", "pathwaySource")])

                freq_unique_input_RampId_G <- as.data.frame(table(unique_input_RampId_G[,"pathwayRampId"]))

                names(freq_unique_input_RampId_G)[1] = 'pathwayRampId'
                merge_PathwayG_source <- merge(freq_unique_input_RampId_G, unique_pathwayG_source, by="pathwayRampId")

            # subset gene data based on source -  kegg, reactome, wiki

                input_kegg_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "kegg")

                input_reactome_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "reactome")

                input_wiki_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "wiki")
            }

            # Loop through each pathway, build the contingency table, and calculate Fisher's Exact
            # test p-value
            pval=totinpath=userinpath=pidused=c()
            for (i in pid) {

                if(analyte_type=="metabolites") {

                if ((!is.na(input_kegg_metab$pathwayRampId[1])) && i %in% input_kegg_metab$pathwayRampId) {
                    tot_in_pathway <- input_kegg_metab[which(input_kegg_metab[,"pathwayRampId"]==i),"Freq"]
                    total_analytes <- kegg_totanalytes
                } else if ((!is.na(input_wiki_metab$pathwayRampId[1])) && i %in% input_wiki_metab$pathwayRampId) {
                    tot_in_pathway <- input_wiki_metab[which(input_wiki_metab[,"pathwayRampId"]==i),"Freq"]
                    total_analytes <- wiki_totanalytes
                } else if ((!is.na(input_reactome_metab$pathwayRampId[1])) && i %in% input_reactome_metab$pathwayRampId) {
                    tot_in_pathway <- input_reactome_metab[which(input_reactome_metab[,"pathwayRampId"]==i),"Freq"]
                    total_analytes <- react_totanalytes
                } else {
                    tot_in_pathway = 0
                }

                } else {

                if ((!is.na(input_kegg_gene$pathwayRampId[1])) && i %in% input_kegg_gene$pathwayRampId) {
                    tot_in_pathway <- input_kegg_gene[which(input_kegg_gene[,"pathwayRampId"]==i),"Freq"]
                    total_analytes <- kegg_totanalytes
                } else if ((!is.na(input_wiki_gene$pathwayRampId[1])) && i %in% input_wiki_gene$pathwayRampId) {
                    tot_in_pathway <- input_wiki_gene[which(input_wiki_gene[,"pathwayRampId"]==i),"Freq"]
                    total_analytes <- wiki_totanalytes
                } else if ((!is.na(input_reactome_gene$pathwayRampId[1])) &&i %in% input_reactome_gene$pathwayRampId) {
                    tot_in_pathway <- input_reactome_gene[which(input_reactome_gene[,"pathwayRampId"]==i),"Freq"]
                    total_analytes <- react_totanalytes
                } else {
                    tot_in_pathway = 0
                }

                }
                tot_out_pathway <- total_analytes - tot_in_pathway
                # fill the rest of the table out
                user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
                user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
                contingencyTb[1,1] <- tot_in_pathway - user_in_pathway
                contingencyTb[1,2] <- tot_out_pathway - user_out_pathway
                contingencyTb[2,1] <- user_in_pathway
                contingencyTb[2,2] <- user_out_pathway

                # Put the test into a try catch in case there's an issue, we'll have some details on the contingency matrix
                tryCatch({
                    result <- stats::fisher.test(contingencyTb)
                }, error = function(e) {
                    print(toString(e))
                    print(i)
                    print(contingencyTb)
                })

                pval <- c(pval,result$p.value )
                userinpath<-c(userinpath,user_in_pathway)
                totinpath<-c(totinpath,tot_in_pathway)
                pidused <- c(pidused,i)
            } # end for loop


            # Now run fisher's tests for all other pids
            allpids <- analyte_pathway$get_pathway_ramp_ids(not_in_pathway_sources=c("hmdb"))
            pidstorun <- setdiff(allpids[,1],pid)
            pidstorunlist <- sapply(pidstorun,shQuote)
            pidstorunlist <- paste(pidstorunlist,collapse = ",")

            #print(paste0(length(pidstorun),"pathways"))

            # We're collecting p-values for all pathways, now those with no analyte support at all - JCB:?

            # calculating p-values for all other pathways
            kegg_metab <- kegg_metab
            kegg_gene <- kegg_gene
            wiki_metab <- wiki_metab
            wiki_gene <- wiki_gene
            reactome_metab <- reactome_metab
            reactome_gene <- reactome_gene
            hmdb_metab <- hmdb_metab
            hmdb_gene <- hmdb_gene
            count=1;
            pval2=userinpath2=totinpath2=c()

            for (i in pidstorun) {
            if(( count %% 100) ==0) {print(paste0("Processed ",count))}
            count=count+1
            user_in_pathway=0
            if(analyte_type=="metabolites") {

                if (i %in% kegg_metab$pathwayRampId) {
                tot_in_pathway <- kegg_metab[which(kegg_metab[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- kegg_totanalytes

                } else if (i %in% wiki_metab$pathwayRampId) {
                tot_in_pathway <- wiki_metab[which(wiki_metab[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- wiki_totanalytes

                } else if (i %in% reactome_metab$pathwayRampId) {
                tot_in_pathway <- reactome_metab[which(reactome_metab[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- react_totanalytes

                } else if (i %in% hmdb_metab$pathwayRampId) {
                tot_in_pathway <- hmdb_metab[which(hmdb_metab[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- NULL
                } else {
                tot_in_pathway=0
                total_analytes <- NULL
                }

            } else {

                if (i %in% kegg_gene$pathwayRampId) {
                tot_in_pathway <- kegg_gene[which(kegg_gene[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- kegg_totanalytes
                } else if (i %in% wiki_gene$pathwayRampId) {
                tot_in_pathway <- wiki_gene[which(wiki_gene[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- wiki_totanalytes
                } else if (i %in% reactome_gene$pathwayRampId) {
                tot_in_pathway <- reactome_gene[which(reactome_gene[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- react_totanalytes
                } else if (i %in% hmdb_gene$pathwayRampId) {
                tot_in_pathway <- hmdb_gene[which(hmdb_gene[,"pathwayRampId"]==i),"Freq"]
                total_analytes <- NULL
                } else {
                tot_in_pathway=0
                total_analytes <- NULL
                }
            }

            # Check that the pathway being considered has your analyte type, if not, move on

                if(is.null(total_analytes)) {next;}
                tot_out_pathway <- total_analytes - tot_in_pathway
                # fill the rest of the table out

                # JCB: Another issue 10/7/2020
                # This line was used for user_out_pathway
                # This section of code is for all pathways that have no analyte support.
                user_out_pathway <- length(unique(pathwaydf$rampId))

                # This line was commented out in production *but* now we have total_analytes set properly
                # not sure why this line was changed.
                # user_out_pathway <- total_analytes - user_in_pathway

                contingencyTb[1,1] <- tot_in_pathway - user_in_pathway
                contingencyTb[1,2] <- tot_out_pathway - user_out_pathway
                contingencyTb[2,1] <- user_in_pathway
                contingencyTb[2,2] <- user_out_pathway

                # Added try catch
                tryCatch({
                result <- stats::fisher.test(contingencyTb)
                }, error = function(e) {
                print(toString(e))
                print(i)
                print(contingencyTb)
                })

                pval2 <- c(pval2,result$p.value )
                userinpath2<-c(userinpath2,user_in_pathway)
                totinpath2<-c(totinpath2,tot_in_pathway)
                # pidused <- c(pidused,i)
            } # end for loop

            # only keep pathways that have > 8 or < 100 compounds
            keepers <- intersect(which(c(totinpath,totinpath2)>=8),
                                which(c(totinpath,totinpath2)<100))
            #hist(totinpath,breaks=1000)
            print(paste0("Keeping ",length(keepers)," pathways"))
            #fdr <- stats::p.adjust(c(pval,pval2)[keepers],method="fdr")
            #holm <- stats::p.adjust(c(pval,pval2)[keepers],method="holm")
            print(paste0("Calculated p-values for ",length(c(pval,pval2))," pathways"))

            # format output (retrieve pathway name for each unique source id first
            out <- data.frame(pathwayRampId=c(pidused,pidstorun)[keepers],
                                Pval=c(pval,pval2)[keepers],   #FDR.Adjusted.Pval=fdr,
                                # Holm.Adjusted.Pval=holm,
                                Num_In_Path=c(userinpath,userinpath2)[keepers],
                                Total_In_Path=c(totinpath,totinpath2)[keepers])
            print(dim(out))
            #out2 <- merge(out,pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
            #	"pathwaysource","pathwayRampId")],
            #	by="pathwayRampId",all.x=TRUE)
            #finout <- out[,c("pathwayName", "Pval", #"FDR.Adjusted.Pval",
            # "Holm.Adjusted.Pval",
            #	"pathwaysourceId",
            #	"pathwaysource","Num_In_Path","Total_In_Path","pathwayRampId")]
            #  finout=finout[!duplicated(finout),]
            out = out[!duplicated(out),]
            print(colnames(fishers_results))

            # foruser is the output needed, based on what user input
            return(out)
        },
        #' Perform fuzzy multiple linkage partitioning clustering on pathways identified by
        #' Fisher's test
        #'
        #' @param fishers_df The data frame generated by runFisherTest
        #' @param perc_analyte_overlap Minimum overlap for pathways to be considered similar
        #' (Default = 0.5)
        #' @param min_pathway_tocluster Minimum number of 'similar' pathways required to start
        #' a cluster (medoid) (Default = 3)
        #' @param perc_pathway_overlap Minimum overlap for clusters to merge (Default = 0.5)
        #'
        #' @return list:[[1]] Pathway enrichment result dataframe with cluster assignment column added
        #' [[2]] analyte type
        #' [[3]] cluster assignment in the list form
        #'@examples
        #'\dontrun{
        #' pathwaydf<-getPathwayFromAnalyte(c("MDM2","TP53","glutamate","creatinine"),
        #'                 NameOrIds="names", conpass=conpass)
        #' fisher.results <- runCombinedFisherTest(pathwaydf=pathwaydf,conpass=conpass)
        #' filtered.fisher.results <- FilterFishersResults(fisher.results,p_holmadj_cutoff=0.05)
        #' filteredclust.fisher.results <- findCluster(filtered.fisher.results)
        #'}
        #' @export
        run_combined_fisherTest = function(
            pathwaydf,
            total_metabolites=NULL,
            total_genes=20000,
            min_analyte=2
        ) {
            G <- M <- 0

            # Grab pathways that contain metabolites to run Fisher on metabolites
            # This will return all pathways that have at 8-120 metabolites/genes in them
            fishmetab <- pathwaydf[grep("RAMP_C_",pathwaydf$rampId),]
            if(nrow(fishmetab) == 0) {outmetab=NULL} else{
                M=1
                print("Running Fisher's tests on metabolites")
                outmetab <- run_fisher_test(pathwaydf=fishmetab,analyte_type="metabolites",
                                        total_metabolites=total_metabolites,total_genes=total_genes,
                                        conpass=conpass,dbname=dbname,
                                        username=username,host=host)
            }

            # Grab pathways that contain genes to run Fisher on genes
            fishgene <- pathwaydf[grep("RAMP_G_",pathwaydf$rampId),]
            if(nrow(fishgene) == 0) {outgene=NULL} else{
                G=1
                print("Running Fisher's tests on genes")
                outgene <- run_fisher_test(pathwaydf=fishgene,analyte_type="genes",
                                        total_metabolites=total_metabolites,total_genes=total_genes,
                                        conpass=conpass,dbname=dbname,
                                        username=username,host=host)
            }

            if(is.null(outgene) & !is.null(outmetab)) {
                out <- outmetab
                fdr <- stats::p.adjust(out$Pval,method="fdr")
                out<-cbind(out,fdr);colnames(out)[ncol(out)]="Pval_FDR"
                holm <- stats::p.adjust(out$Pval,method="holm")
                out<-cbind(out,holm);colnames(out)[ncol(out)]="Pval_Holm"
                keepers <- which(out$Num_In_Path>=min_analyte)
                out2 <- merge(out[keepers,],
                            pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
                                        "pathwaysource")],by="pathwayRampId")
            } else if (!is.null(outgene) & is.null(outmetab)) {
                out <- outgene
                fdr <- stats::p.adjust(out$Pval,method="fdr")
                out<-cbind(out,fdr);colnames(out)[ncol(out)]="Pval_FDR"
                holm <- stats::p.adjust(out$Pval,method="holm")
                out<-cbind(out,holm);colnames(out)[ncol(out)]="Pval_Holm"
                keepers <- which(out$Num_In_Path>=min_analyte)
                out2 <- merge(out[keepers,],
                            pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
                                        "pathwaysource")],by="pathwayRampId")
            } else {
                # merge the results if both genes and metabolites were run
                G = M = 1
                allfish <- merge(outmetab,outgene,
                                by="pathwayRampId",all.x=T,all.y=T)
                colnames(allfish)[which(colnames(allfish)=="Pval.x")]="Pval.Metab"
                colnames(allfish)[which(colnames(allfish)=="Pval.y")]="Pval.Gene"
                colnames(allfish)[which(colnames(allfish)=="Total_In_Path.x")]="Total_In_Path.Metab"
                colnames(allfish)[which(colnames(allfish)=="Total_In_Path.y")]="Total_In_Path.Gene"
                colnames(allfish)[which(colnames(allfish)=="Num_In_Path.x")]="Num_In_Path.Metab"
                colnames(allfish)[which(colnames(allfish)=="Num_In_Path.y")]="Num_In_Path.Gene"

                # Calculate combined p-values for pathways that have both genes and metabolites
                gm <- intersect(which(!is.na(allfish$Pval.Metab)),which(!is.na(allfish$Pval.Gene)))
                combpval <- stats::pchisq(-2 * (log(allfish$Pval.Metab[gm])+log(allfish$Pval.Gene[gm])),
                                        df=2,lower.tail=FALSE)

                g <- which(is.na(allfish$Pval.Metab))
                gpval <- allfish$Pval.Gene[g]
                m <- which(is.na(allfish$Pval.Gene))
                mpval <- allfish$Pval.Metab[m]

                out <- rbind(allfish[gm,],allfish[g,],allfish[m,])
                out <- cbind(out,c(combpval,gpval,mpval))
                colnames(out)[ncol(out)]="Pval_combined"
                fdr <- stats::p.adjust(out$Pval_combined,method="fdr")
                out <- cbind(out,fdr)
                colnames(out)[ncol(out)]="Pval_combined_FDR"
                holm <- stats::p.adjust(out$Pval_combined,method="holm")
                out <- cbind(out,holm)
                colnames(out)[ncol(out)]="Pval_combined_Holm"

                keepers <- intersect(c(which(out$Num_In_Path.Metab>=min_analyte),
                                    which(is.na(out$Num_In_Path.Metab))),
                                    c(which(out$Num_In_Path.Gene>=min_analyte),
                                    which(is.na(out$Num_In_Path.Gene)))
                )


                # Now that p-values are calculated, only return pathways that are in the list
                # of pathways that contain user genes and metabolites
                out2 <- merge(out[keepers,],
                            pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
                                        "pathwaysource")],by="pathwayRampId")
            } # end merging when genes and metabolites were run
            fishresults <- out2[!duplicated(out2),]

            analyte_type=c()
            if(G==1 && M==1) {
                analyte_type="both"
            } else if (G==1 && M==0) {
                analyte_type="genes"
            } else if (G==0 && M==1) {
                analyte_type="metabolites"
            }

            return(list(fishresults=fishresults,analyte_type=analyte_type))
        },
        #' @param fishers_df The data frame generated by runFisherTest
        #' @param p_holmadj_cutoff return pathways where Holm adjusted pvalues are < p_holmadj_cutoff
        #' @param p_fdradj_cutoff return pathways where FDR adjusted pvalues are < p_fdradj_cutoff
        #' @return list:[[1]]Dataframe with pathway enrichment results, only significant pathways
        #' [[2]]analyte type
        #'@examples
        #'\dontrun{
        #' pathwaydf<-getPathwayFromAnalyte(c("MDM2","TP53","glutamate","creatinine"),
        #'                 NameOrIds="names", conpass=conpass)
        #' fisher.results <- runCombinedFisherTest(pathwaydf=pathwaydf,conpass=conpass)
        #' filtered.fisher.results <- FilterFishersResults(fisher.results,p_holmadj_cutoff=0.05)
        #'}
        #' @export
        filter_fishers_results = function(
            p_holmadj_cutoff=NULL,
            p_fdradj_cutoff=NULL
        ){

        if(length(grep("Pval_combined",colnames(fishresults)))==0) {
            if(!is.null(p_holmadj_cutoff)) {
            return(list(fishresults=fishresults[which(fishresults[,"Pval_Holm"] <=
                                                        p_holmadj_cutoff),],analyte_type=analyte_type))
            } else if (!is.null(p_fdradj_cutoff)) {
            return(list(fishresults=fishresults[which(fishresults[,"Pval_FDR"] <=
                                                        p_fdradj_cutoff),],analyte_type=analyte_type))
            } else {
            stop("Please set a cutoff for Holm Adjusted pvalues
                    (p_holmadj_cutoff paramter) or FDR Adjusted pvalues
                    (p_fdradj_cutoff)")
            }
        }  else { # ORA was performed on both genes and metabolites:
            if(!is.null(p_holmadj_cutoff)) {
            return(list(fishresults=fishresults[which(fishresults[,"Pval_combined_Holm"] <=
                                                        p_holmadj_cutoff),],analyte_type=analyte_type))
            } else if (!is.null(p_fdradj_cutoff)) {
            return(list(fishresults=fishresults[which(fishresults[,"Pval_combined_FDR"] <=
                                                        p_fdradj_cutoff),],analyte_type=analyte_type))
            } else {
            stop("Please set a cutoff for Holm Adjusted pvalues
                                (p_holmadj_cutoff paramter) or FDR Adjusted pvalues
                                (p_fdradj_cutoff)")
            }
        }
        }
    )
)
