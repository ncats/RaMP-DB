
.onLoad <- function(libname, pkgname) {
}

.onAttach <- function(libname, pkgname) {

    message("Loading Ramp version: ", packageDescription("RaMP")$Version)
    db_local <- listRaMPVersions(local = TRUE)
    if (length(db_local)) {
        message("Locally cached RaMP Database versions:")
        message(paste0(" - ", db_local, "\n"))
    } else {
        message("No cached RaMP databases present. Use 'RaMP()' to download ",
                "the most recent release. See '?RaMP' for more information.")
    }
}



