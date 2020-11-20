library(methods)

#' @import RMySQL
#' @export DBConnection
DBConnection <- setRefClass("DBConnection",
    fields = list(
        con = "MySQLConnection"
    ),
    methods = list(
        initialize = function() {
            con <<- DBI::dbConnect(
                RMySQL::MySQL(),
                user = .GlobalEnv$db_username,
                dbname = .GlobalEnv$db_dbname,
                password = .GlobalEnv$db_password,
                host = .GlobalEnv$db_host
            )
        },
        disconnect = function() {
            DBI::dbDisconnect(con)
        },
        run_query = function(query) {
            return(DBI::dbGetQuery(con, query))
        }
    )
)
