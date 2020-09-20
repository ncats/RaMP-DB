library(plumber)

host <- "ramp-db.ncats.io"
dbname <- "ramp"
username <- "ramp_query_user"
conpass <- "ramp_query_user"

#* @filter cors
cors <- function(res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    plumber::forward()
}

#* Return analytes from pathway
#* @param pathway
#* @get /api/analyte
function(pathway="") {
    ramp_out <- getAnalyteFromPathway(
        pathway,
        conpass = conpass,
        host = host,
        dbname = dbname,
        username = username
    )
    ramp_out
}


#* Return pathway enrichment analysis results
#* @param analyte
#* @param identifier_type names or ids
#* @param p_holmadj_cutoff
#* @param p_fdradj_cutoff
#* @get /api/pathway_enrichment_analysis
function(
    analyte="",
    identifier_type="names",
    p_holmadj_cutoff=0.05,
    p_fdradj_cutoff=NULL
) {

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
    filtered_results
}

#' Serve the default HTML file
#' @get /
#' @get /pathway-enrichment-analysis
function(res) {
    plumber::include_html("../../client/index.html", res)
}

#' @assets ../../client /client
list()

#' @assets ../../client/assets /assets
list()
