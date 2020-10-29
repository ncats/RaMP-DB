library(plumber)
library(sqldf)

host <- "ramp-db.ncats.io"
dbname <- "ramp"
username <- "ramp_query_user"
conpass <- "ramp_query_user"

#* @filter cors
cors <- function(res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    plumber::forward()
}

#* Return analyte source intersects
#* @serializer unboxedJSON
#* @get /api/analyte_intersects
function() {

    intersects <- list(
        compounds=list(
            list(
                sets=list("KEGG"),
                size=9
            ),
            list(
                sets=list("REACTOME"),
                size=226
            ),
            list(
                sets=list("WP"),
                size=631
            ),
            list(
                sets=list("HMDB"),
                size=106778
            ),
            list(
                sets=list("KEGG", "REACTOME"),
                size=20
            ),
            list(
                sets=list("KEGG", "WP"),
                size=125
            ),
            list(
                sets=list("KEGG", "HMDB"),
                size=3791
            ),
            list(
                sets=list("REACTOME", "WP"),
                size=817
            ),
            list(
                sets=list("REACTOME", "HMDB"),
                size=675
            ),
            list(
                sets=list("WP", "HMDB"),
                size=591
            ),
            list(
                sets=list("KEGG", "REACTOME", "WP"),
                size=142
            ),
            list(
                sets=list("REACTOME", "WP", "HMDB"),
                size=973
            ),
            list(
                sets=list("KEGG", "REACTOME", "HMDB"),
                size=940
            ),
            list(
                sets=list("KEGG", "WP", "HMDB"),
                size=19
            ),
            list(
                sets=list("KEGG", "REACTOME", "WP", "HMDB"),
                size=33
            )
        ),
        genes=list(
            list(
                sets=list("KEGG"),
                size=936
            ),
            list(
                sets=list("REACTOME"),
                size=278
            ),
            list(
                sets=list("WP"),
                size=975
            ),
            list(
                sets=list("HMDB"),
                size=1207
            ),
            list(
                sets=list("KEGG", "REACTOME"),
                size=171
            ),
            list(
                sets=list("KEGG", "WP"),
                size=1026
            ),
            list(
                sets=list("KEGG", "HMDB"),
                size=1781
            ),
            list(
                sets=list("REACTOME", "WP"),
                size=3825
            ),
            list(
                sets=list("REACTOME", "HMDB"),
                size=19
            ),
            list(
                sets=list("WP", "HMDB"),
                size=14
            ),
            list(
                sets=list("KEGG", "REACTOME", "WP"),
                size=4441
            ),
            list(
                sets=list("REACTOME", "WP", "HMDB"),
                size=1088
            ),
            list(
                sets=list("KEGG", "REACTOME", "HMDB"),
                size=28
            ),
            list(
                sets=list("KEGG", "WP", "HMDB"),
                size=21
            ),
            list(
                sets=list("KEGG", "REACTOME", "WP", "HMDB"),
                size=1544
            )
        )
    )
}

#* Return analytes from pathway
#* @param pathwaySourceId
#* @get /api/analytes
function(pathwaySourceId="") {
    con <- DBI::dbConnect(RMySQL::MySQL(),
                          user = username,
                          dbname=dbname,
                          password = conpass,
                          host = host)
    on.exit(DBI::dbDisconnect(con))
    source_ids <- sapply(pathwaySourceId,shQuote)
    source_ids <- paste(source_ids,collapse = ",")

    query <- paste0(
        "select p.sourceId, p.pathwayName, GROUP_CONCAT(s.sourceId) as analytes ",
        "from pathway as p ",
        "left join analytehaspathway as ap on p.pathwayRampId = ap.pathwayRampId ",
        "join source as s on ap.rampId = s.rampId ",
        "where p.sourceId in (", source_ids, ") ",
        "group by p.sourceId, p.pathwayName"
    )

    cids <- DBI::dbGetQuery(con,query)
    return(cids)
}


#* Return pathway enrichment analysis results
#* @param analyte
#* @param identifier_type names or ids
#* @param p_holmadj_cutoff
#* @param p_fdradj_cutoff
#* @perc_analyte_overlap
#* @perc_pathway_overlap
#* @min_pathway_tocluster
#* @get /api/pathway_enrichment_analysis
function(
    analyte="",
    identifier_type="names",
    p_holmadj_cutoff=0.05,
    p_fdradj_cutoff=NULL,
    perc_analyte_overlap=0.2,
    perc_pathway_overlap=0.2,
    min_pathway_tocluster=2
) {

    if (typeof(min_pathway_tocluster) == "character") {
        min_pathway_tocluster = strtoi(min_pathway_tocluster, base = 0L)
    }

    con <- DBI::dbConnect(RMySQL::MySQL(),
                        user = username,
                        dbname=dbname,
                        password = conpass,
                        host = host)
    on.exit(DBI::dbDisconnect(con))
    pathways <- getPathwayFromAnalyte(
        analytes = c(analyte),
        con = con,
        NameOrIds = identifier_type
    )
    fishers_results <- runCombinedFisherTest(
        pathways,
        con = con
    )
    filtered_results <- FilterFishersResults(
        fishers_results,
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
    ids_no_cluster <- fishresults[fishresults$cluster_assignment != 'Did not cluster', 'pathwayRampId']
    pathway_matrix <- clustering_results$pathway_matrix[ids_no_cluster,ids_no_cluster]
    distance_matrix <- dist(1 - pathway_matrix)

    fit <- cmdscale(distance_matrix,eig=TRUE, k=2)

    cluster_coordinates <- data.frame(fit$points)
    cluster_coordinates <- cbind(pathwayRampId = rownames(cluster_coordinates), cluster_coordinates)
    rownames(cluster_coordinates) <- NULL

    names(cluster_coordinates)[2] <- "x"
    names(cluster_coordinates)[3] <- "y"

    cluster_coordinates <- sqldf("
                    select
                        cluster_coordinates.pathwayRampId,
                        cluster_coordinates.x,
                        cluster_coordinates.y,
                        fishresults.cluster_assignment
                    from cluster_coordinates
                    left join fishresults on cluster_coordinates.pathwayRampId = fishresults.pathwayRampId
                   ")

    analyte_ids <- sapply(analyte,shQuote)
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

    cids <- DBI::dbGetQuery(con,query)

    response <- list(fishresults=clustering_results$fishresults, clusterCoordinates=cluster_coordinates, analytes=cids)

    return(response)
}

#' Serve the default HTML file
#' @get /
#' @get /pathway-enrichment-analysis
#' @get /about
function(res) {
    plumber::include_html("../../client/index.html", res)
}

#' @assets ../../client /client
list()

#' @assets ../../client/assets /assets
list()
