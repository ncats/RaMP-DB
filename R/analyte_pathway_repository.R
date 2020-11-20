library(methods)

#' @export AnalytePathwayRepository
AnalytePathwayRepository <- setRefClass("AnalytePathwayRepository",
    fields = list (
        db_connection = "list"
    ),
    methods = list(
        initialize = function() {
            if (is.null(db_connection)) {
                stop("Please provide a db_connection")
            }
        },
        get_analytes_pathways = function(
            analyte_ids, pathway_ramp_ids
        ) {
            base_query <- "select * from analytehaspathway"

            predicates <- c()

            if (!is.null(analyte_ids)) {
                analyte_ids_query <- sapply(analyte_ids, shQuote)
                analyte_ids_query <- paste(analyte_ids_query, collapse = ",")
                analyte_ids_query <- paste0("rampId in (", analyte_ids_query, ")")
                predicates <- c(predicates, analyte_ids_query)
            }

            if (!is.null(pathway_ramp_ids)) {
                pathway_ids_query <- sapply(pathway_ramp_ids, shQuote)
                pathway_ids_query <- paste(pathway_ids_query, collapse = ",")
                pathway_ids_query <- paste0("pathwayRampId in (", pathway_ids_query, ")")
                predicates <- c(predicates, pathway_ids_query)
            }

            predicates_query <- paste(predicates, collapse = " and ")

            query <- paste(base_query, predicates_query, sep =  " where ")

            analyte_pathways <- db_connection$run_query(query)

            return(analyte_pathways)
        },
        get_pathway_ramp_ids = function(
            distinct = TRUE,
            in_pathway_sources = NULL,
            not_in_pathway_sources = NULL
        ) {

            distinct_query_1 = ""
            distinct_query_2 = ""

            if (distinct == TRUE) {
                distinct_query_1 <- "distinct("
                distinct_query_2 <- ")"
            }

            query <- paste0(
                "select ",
                distinct_query_1,
                "pathwayRampId",
                distinct_query_2,
                "from analytehaspathway"
            )

            if (!is.null(in_pathway_sources)) {
                in_pathway_sources_query <- sapply(in_pathway_sources, shQuote)
                in_pathway_sources_query <- paste(in_pathway_sources_query, collapse = ",")
                query <- paste0(query, "where pathwaySource in ", in_pathway_sources_query)
            } else if (!is.null(not_in_pathway_sources)) {
                not_in_pathway_sources_query <- sapply(not_in_pathway_sources, shQuote)
                not_in_pathway_sources_query <- paste(in_pathway_sources_query, collapse = ",")
                query <- paste0(query, "where pathwaySource not in ", not_in_pathway_sources_query)
            }

            pathway_ramp_ids = db_connection$run_query(query)

            return(pathway_ramp_ids)
        }
    )
)
