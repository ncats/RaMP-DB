library(plumber)
library(sqldf)
library(config)

host <- "ramp-db.ncats.io"
dbname <- "ramp"
username <- "ramp_query_user"
conpass <- "ramp_query_user"

#* @filter cors
cors <- function(req, res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200
    return(list())
  } else {
    plumber::forward()
  }
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

#* Return analytes from source database
#* @param identifier
#* @serializer unboxedJSON
#* @get /api/source/analytes
function(identifier) {
    identifiers <- c(identifier)
    identifiers <- sapply(identifiers,shQuote)
    identifiers <- paste(identifiers, collapse=",")
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    con <- DBI::dbConnect(RMySQL::MySQL(),
                          user = username,
                          dbname = dbname,
                          password = conpass,
                          host = host)
    query <- paste0(
        "select s.sourceId, s.IDtype, s.geneOrCompound, s.commonName, min(ansyn.Synonym) as synonym ",
        "from source as s ",
        "left join analytesynonym as ansyn on s.rampId = ansyn.rampId and ansyn.Synonym in (", identifiers, ") ",
        "where s.sourceId in (", identifiers, ") ",
        "or s.commonName in (", identifiers, ") ",
        "or s.rampId in (",
            "select analytesynonym.rampId ",
            "from analytesynonym ",
            "where analytesynonym.Synonym in (", identifiers, ")",
        ") ",
        "group by s.sourceId, s.IDtype, s.geneOrCompound, s.commonName"
    )
    analytes <- DBI::dbGetQuery(con,query)
    DBI::dbDisconnect(con)
    return(analytes)
}

#* Return analytes from pathway
#* @param pathwaySourceId
#* @get /api/analytes_test
function(pathwaySourceId="") {
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    con <- DBI::dbConnect(RMySQL::MySQL(),
                          user = username,
                          dbname=dbname,
                          password = conpass,
                          host = host)
    on.exit(DBI::dbDisconnect(con))
    source_ids <- sapply(pathwaySourceId, shQuote)
    source_ids <- paste(source_ids, collapse = ",")

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

#* Return analytes from pathway
#* @param analyte
#* @get /api/analytes
function(analyte="") {
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password

    analytes <- c(analyte)

    analytes_df <- getOntoFromMeta(
        analytes = analytes,
        conpass = conpass,
        host = host,
        dbname = dbname,
        username = username
    )

    return(analytes_df)
}

#' Return pathway enrichment analysis results
#' @param analyte
#' @param identifier_type names or ids
#' @param p_holmadj_cutoff
#' @param p_fdradj_cutoff
#' @param perc_analyte_overlap
#' @param perc_pathway_overlap
#' @param min_pathway_tocluster
#' @get /api/pathway_enrichment_analysis
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
        min_pathway_tocluster <- strtoi(min_pathway_tocluster, base = 0L)
    }

    analytes <- c(analyte)

    analysisResults <- pathway_enrichment_analysis(
        analytes,
        identifier_type,
        p_holmadj_cutoff,
        p_fdradj_cutoff,
        perc_analyte_overlap,
        perc_pathway_overlap,
        min_pathway_tocluster
    )

    return(list(
        fishresults = analysisResults$fishresults,
        clusterCoordinates = analysisResults$cluster_coordinates,
        analytes = analysisResults$cids
    ))
}

#' Return pathways from given list of analytes
#' @param analyte
#' @get /api/pathways
function(analyte="") {
    analytes <- c(analyte)
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    pathways_df_ids <- tryCatch(
        {
            pathways_df <- getPathwayFromAnalyte(
                analytes = analytes,
                conpass=conpass,
                host=host,
                dbname=dbname,
                username=username,
                NameOrIds = 'ids'
            )
        },
        error=function(cond) {
            print(cond)
            return(data.frame(stringsAsFactors=FALSE))
        }
    )
    pathways_df_names <- tryCatch(
        {
            pathways_df <- getPathwayFromAnalyte(
                analytes = analytes,
                conpass=conpass,
                host=host,
                dbname=dbname,
                username=username,
                NameOrIds = 'names'
            )
        },
        error=function(cond) {
            print(cond)
            return(data.frame(stringsAsFactors=FALSE))
        }
    )
    pathways_df <- rbind(pathways_df_ids, pathways_df_names)
    return(unique(pathways_df))
}

#' Return combined Fisher's test results from given list of analytes query results
#' @parser json
#' @post /api/combined-fisher-test
function(req) {
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    pathways_df <- as.data.frame(req$body)
    fishers_results_df <- runCombinedFisherTest(
        pathwaydf = pathways_df,
        conpass=conpass,
        host=host,
        dbname=dbname,
        username=username
    )
    return(fishers_results_df)
}

#' Return filtered Fisher's test results from given list of Fisher's test results
#' @param p_holmadj_cutoff
#' @param p_fdradj_cutoff
#' @parser json
#' @post /api/filter-fisher-test-results
function(req, p_holmadj_cutoff=0.05, p_fdradj_cutoff=NULL) {
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    fishers_results <- req$body
    fishers_results$fishresults <- as.data.frame(fishers_results$fishresults)
    filtered_results <- FilterFishersResults(
        fishers_df=fishers_results,
        p_holmadj_cutoff = p_holmadj_cutoff,
        p_fdradj_cutoff = p_fdradj_cutoff
    )
    return(filtered_results)
}

#' Return filtered Fisher's test results from given list of Fisher's test results
#' @param perc_analyte_overlap
#' @param perc_pathway_overlap
#' @param min_pathway_tocluster
#' @parser json
#' @post /api/cluster-fisher-test-results
function(req, analyte_source_id, perc_analyte_overlap=0.2, perc_pathway_overlap=0.2, min_pathway_tocluster=2) {
    analytes <- c(analyte_source_id)
    if (typeof(min_pathway_tocluster) == "character") {
        min_pathway_tocluster <- strtoi(min_pathway_tocluster, base = 0L)
    }
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    fishers_results <- req$body
    fishers_results$fishresults <- as.data.frame(fishers_results$fishresults)
    clustering_results <- findCluster(
        fishers_results,
        perc_analyte_overlap=perc_analyte_overlap,
        min_pathway_tocluster=min_pathway_tocluster,
        perc_pathway_overlap=perc_pathway_overlap
    )

    return(clustering_results)
}

#' Return filtered Fisher's test results from given list of Fisher's test results
#' @param analyte_source_id
#' @param perc_analyte_overlap
#' @param perc_pathway_overlap
#' @param min_pathway_tocluster
#' @parser json
#' @post /api/cluster-fisher-test-results-extended
function(req, analyte_source_id, perc_analyte_overlap=0.2, perc_pathway_overlap=0.2, min_pathway_tocluster=2) {
    analytes <- c(analyte_source_id)
    if (typeof(min_pathway_tocluster) == "character") {
        min_pathway_tocluster <- strtoi(min_pathway_tocluster, base = 0L)
    }
    config <- config::get()
    host <- config$db_host
    dbname <- config$db_dbname
    username <- config$db_username
    conpass <- config$db_password
    fishers_results$fishresults <- as.data.frame(fishers_results$fishresults)
    clustering_results <- findCluster(
        fishers_results,
        perc_analyte_overlap=perc_analyte_overlap,
        min_pathway_tocluster=min_pathway_tocluster,
        perc_pathway_overlap=perc_pathway_overlap
    )

    # return(clustering_results)
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
                fishresults.cluster_assignment,
                fishresults.pathwayName
            from cluster_coordinates
            left join fishresults on
                cluster_coordinates.pathwayRampId = fishresults.pathwayRampId"
        )
    }

    analyte_ids <- sapply(analytes,shQuote)
    analyte_ids <- paste(analyte_ids,collapse = ",")

    query <- paste0(
        "select s.sourceId, commonName, GROUP_CONCAT(p.sourceId) as pathways ",
        "from source as s ",
        "left join analyte as a on s.rampId = a.rampId ",
        "left join analytehaspathway as ap on a.rampId = ap.rampId ",
        "left join pathway as p on ap.pathwayRampId = p.pathwayRampId ",
        "where s.sourceId in (", analyte_ids, ") ",
        "group by s.sourceId, s.commonName"
    )

    con <- RaMP::connectToRaMP(dbname=dbname,username=username,conpass=conpass,host = host)
    cids <- DBI::dbGetQuery(con,query)
    DBI::dbDisconnect(con)

    response <- list(
        fishresults = clustering_results$fishresults,
        clusterCoordinates = cluster_coordinates,
        analytes = cids
    )

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
