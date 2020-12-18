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


#' Return pathway enrichment analysis results
#' @param analyte
#' @param identifier_type names or ids
#' @param p_holmadj_cutoff
#' @param p_fdradj_cutoff
#' @param perc_analyte_overlap
#' @param perc_pathway_overlap
#' @param min_pathway_tocluster
#' @get /api/pathway_enrichment_analysis
pathway_enrichment_analysis

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
