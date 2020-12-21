
#' Pathway enrichment analysis
#' @import stats
#' @import sqldf
#' @param analyte vector of analyte names or Id's
#' @param identifier_type names or ids
#' @param p_holmadj_cutoff p_holmadj_cutoff
#' @param p_fdradj_cutoff p_fdradj_cutoff
#' @param perc_analyte_overlap perc_analyte_overlap
#' @param perc_pathway_overlap perc_pathway_overlap
#' @param min_pathway_tocluster min_pathway_tocluster
#' @return list with pathway enrichment analysis results
pathway_enrichment_analysis <- function(
    analytes="",
    identifier_type="names",
    p_holmadj_cutoff=0.05,
    p_fdradj_cutoff=NULL,
    perc_analyte_overlap=0.2,
    perc_pathway_overlap=0.2,
    min_pathway_tocluster=2,
    host=NULL,
    dbname=NULL,
    username=NULL,
    conpass=NULL
) {

    if (is.null(host)) {
        host <- .GlobalEnv$db_host
    }
    if (is.null(dbname)) {
        dbname <- .GlobalEnv$db_dbname
    }
    if (is.null(username)) {
        username <- .GlobalEnv$db_username
    }
    if (is.null(conpass)) {
        conpass <- .GlobalEnv$db_password
    }
    

    pathways_df <- getPathwayFromAnalyte(
        analytes = analytes,
        conpass=conpass,
        host=host,
        dbname=dbname,
        username=username,
        NameOrIds = identifier_type
    )
    fishers_results_df <- runCombinedFisherTest(
        pathwaydf = pathways_df,
        conpass=conpass,
        host=host,
        dbname=dbname,
        username=username
    )
    filtered_results <- FilterFishersResults(
        fishers_df=fishers_results_df,
        p_holmadj_cutoff = p_holmadj_cutoff,
        p_fdradj_cutoff = p_fdradj_cutoff
    )
    clustering_results <- findCluster(
        filtered_results,
        perc_analyte_overlap=perc_analyte_overlap,
        min_pathway_tocluster=min_pathway_tocluster,
        perc_pathway_overlap=perc_pathway_overlap
    )
    fishresults <- clustering_results$fishresults
    ids_no_cluster <- fishresults[
        fishresults$cluster_assignment != 'Did not cluster', 'pathwayRampId'
    ]
    pathway_matrix <- clustering_results$pathway_matrix[ids_no_cluster,ids_no_cluster]

    cluster_coordinates <- c()

    if (!is.null(pathway_matrix)) {
        distance_matrix <- dist(1 - pathway_matrix)

        fit <- cmdscale(distance_matrix, eig=TRUE, k=2)

        cluster_coordinates <- data.frame(fit$points)
        cluster_coordinates <- cbind(
            pathwayRampId = rownames(cluster_coordinates),
            cluster_coordinates
        )
        rownames(cluster_coordinates) <- NULL

        names(cluster_coordinates)[2] <- "x"
        names(cluster_coordinates)[3] <- "y"

        options(sqldf.driver = "SQLite")
        cluster_coordinates <- sqldf(
            "select
                cluster_coordinates.pathwayRampId,
                cluster_coordinates.x,
                cluster_coordinates.y,
                fishresults.cluster_assignment
            from cluster_coordinates
            left join fishresults on
                cluster_coordinates.pathwayRampId = fishresults.pathwayRampId"
        )
    }

    analyte_ids <- sapply(analytes,shQuote)
    analyte_ids <- paste(analyte_ids,collapse = ",")

    where_clause = "where s.sourceId in"

    if (identifier_type == "names") {
        where_clause = "where s.commonName in"
    }

    query <- paste0(
        "select s.sourceId, commonName, GROUP_CONCAT(p.sourceId) as pathways ",
        "from source as s ",
        "left join analyte as a on s.rampId = a.rampId ",
        "left join analytehaspathway as ap on a.rampId = ap.rampId ",
        "left join pathway as p on ap.pathwayRampId = p.pathwayRampId ",
        where_clause, " (", analyte_ids, ") ",
        "group by s.sourceId, s.commonName"
    )
    
    con <- RaMP::connectToRaMP(dbname=dbname,username=username,conpass=conpass,host = host)
    cids <- DBI::dbGetQuery(con,query)
    DBI::dbDisconnect(con)

    response <- list(
        fishresults = clustering_results$fishresults,
        cluster_coordinates = cluster_coordinates,
        analytes = cids
    )

    return(response)
}
