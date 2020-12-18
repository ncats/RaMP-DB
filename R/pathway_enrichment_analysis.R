
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
    analyte="",
    identifier_type="names",
    p_holmadj_cutoff=0.05,
    p_fdradj_cutoff=NULL,
    perc_analyte_overlap=0.2,
    perc_pathway_overlap=0.2,
    min_pathway_tocluster=2
) {
    host <- "ramp-db.ncats.io"
    dbname <- "ramp"
    username <- "ramp_query_user"
    conpass <- "ramp_query_user"

    if (typeof(min_pathway_tocluster) == "character") {
        min_pathway_tocluster <- strtoi(min_pathway_tocluster, base = 0L)
    }
    db_connection <- DBConnection()
    on.exit(db_connection$disconnect())
    analytes <- c(analyte)
    pathways <- getPathwayFromAnalyte(
        analytes = analytes,
        con = db_connection$connection,
        NameOrIds = identifier_type
    )
    fisher_test <- FisherTest(db_connection=db_connection)
    fisher_test$run_combined_fisherTest(pathways)
    filtered_results <- fisher_test$filter_fishers_results(
        p_holmadj_cutoff = p_holmadj_cutoff,
        p_fdradj_cutoff = p_fdradj_cutoff
    )
    fisher_results_clustering <- FisherResultsClustering(fishers_results=filtered_results)
    fisher_results_clustering$find_cluster(
        perc_analyte_overlap=perc_analyte_overlap,
        min_pathway_tocluster=min_pathway_tocluster,
        perc_pathway_overlap=perc_pathway_overlap
    )
    fishresults <- fisher_results_clustering$fishers_results_df
    ids_no_cluster <- fishresults[
        fishresults$cluster_assignment != "Did not cluster", "pathwayRampId"
    ]
    pathway_matrix <- fisher_results_clustering$pathway_matrix[
        ids_no_cluster, ids_no_cluster
    ]

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

    analyte_repository = AnalyteRepository(db_connection=db_connection)

    analyte_names <- NULL
    analyte_ids <- NULL

    if (identifier_type == "names") {
        analyte_names <- analytes
    } else {
        analyte_ids <- analytes
    }

    cids <- analyte_repository$get_analytes(analyte_ids=analyte_ids, analyte_names=analyte_names)

    response <- list(
        fishresults = fisher_results_clustering$fishers_results_df,
        clusterCoordinates = cluster_coordinates,
        analytes = cids
    )

    return(response)
}
