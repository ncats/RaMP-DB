#' Produce RaMP input from metabolite metadata data.frame or csv
#'
#'@description
#' Converts data.frame or csv formatted metabolite metadata into RaMP data input format. The input should have ID sources (e.g. hmdb, kegg, entrez) as column names and the corresponding rows filled with IDs from that source.
#'
#' @details
#' RaMP input format is a vector of entries that look like prefix:ID where the prefix is a code for the ID system the ID belongs to (e.g. hmdb, kegg, entrez).
#' This function expects input with prefixes as column names, and IDs belonging to that prefix as rows in the corresponding column. Data can be supplied as a csv file or data.frame object.
#' NA values, columns named with an usupported prefix, and columns named with a supported prefix but no data will be ignored by the function.
#'
#'  A complete list of currently supported prefixes for metabolites can be found by running the getPrefixesFromAnalytes function, as shown in examples.
#'
#' @param data_frame a data.frame object where the column names are prefixes that will be prepended to each identifier in that column. Please specify either this argument or csv_path.
#' @param csv_path a string containing the file path to a csv that will be read and converted to RaMP data input format. Column names should be prefixes that will be prepended to each identifier in that column.
#' @param ... additional arguments that will be passed to readr::read_csv()
#'
#' @returns A vector object of type character() including each identifier in the input prepended to a supported prefix in the format "prefix:identifier". This result is designed to be input for several RaMP functions.
#'
#' @examples
#'  \dontrun{
#' rampDB <- RaMP()
#'
#' # Example use with data.frame
#' df <- data.frame(ensembl = c('ENSG00000135679', 'ENSG00000141510'),
#'                              hmdb = c('HMDB0000064', NA), fake_ID = c('123', NA))
#'
#' RaMPInput <- createRaMPInput(data_frame = df)
#'
#' getPathwayFromAnalyte(RaMPInput)
#'
#' # Example use with csv
#' RaMPInput <- createRaMPInput(csv_path = 'createRaMPInput_test_data.csv')
#'
#' getPathwayFromAnalyte(RaMPInput)
#'
#' # Retrieve supported ID source lists for metabolites and genes
#' RaMP::getPrefixesFromAnalytes(db = rampDB, analyteType = 'metabolite')
#'
#' RaMP::getPrefixesFromAnalytes(db = rampDB, analyteType = 'gene')
#' }
#'
#' @export
createRaMPInput <- function ( db = RaMP(), data_frame = NULL, csv_path = NULL, ... ) {

  #Produce dataframe of user data from csv path if not provided initially. User must specify exactly one of these.
  if (!is.null(csv_path) & is.null(data_frame)) {
    data_frame <- readr::read_csv(csv_path, ...)
  }

  else if (!is.null(csv_path) & !is.null(data_frame)) {
    stop('You have entered values for both the csv_path and data_frame parameters. Please only specify one of these.')
  }

  else if (is.null(csv_path) & is.null(data_frame)) {
    stop('You have not entered any input data to convert. Please specify either the csv_path or data_frame parameter.')
  }

  #Get acceptable prefixes from RaMP and convert them to character vector.
  prefix_c <- strsplit(getPrefixesFromAnalytes(db = db, analyteType = 'metabolite')$idTypes[1], ',')[[1]]

  prefix_c <- append(prefix_c, strsplit(getPrefixesFromAnalytes(db = db, analyteType = 'gene')$idTypes[1], ',')[[1]])

  prefix_c <- gsub(' ', '', prefix_c)

  prefix_c <- unique(prefix_c)


  #Parse data_frame into RaMP prefix:ID input
  final_id_c <- c()

  for (i in colnames(data_frame)) {

    if (i %in% prefix_c) {
      ids <- data_frame[,i]
      ids <- ids[!is.na(ids)]

      if (length(ids) > 0) {
        final_id_c <- append(final_id_c, paste0(i, ':', ids))
      }
    }

    else {
      warning('Column "', i, '" not processed because identifier name not recognized. Valid names are: ', paste(prefix_c, collapse = ', '))
    }
  }

  return(final_id_c)
}
