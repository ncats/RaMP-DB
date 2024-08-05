#' Produce RaMP input from metabolite metadata data.frame object or .csv/.xlsx files
#'
#' @description
#' Converts data.frame,  .csv, or .xlsx formatted metabolite metadata into RaMP data input format. The input should have ID sources (e.g. hmdb, kegg, entrez) as column names and the corresponding rows filled with IDs from that source.
#'
#' @details
#' RaMP input format is a character() of entries that look like prefix:ID where the prefix is a code for the ID system the ID belongs to (e.g. hmdb, kegg, entrez). This function expects input with prefixes as column names, with IDs belonging to that prefix as rows in the corresponding column. Data can be supplied as a csv file, excel file, or data.frame object. NA values, 'NA' values, columns named with an unsupported prefix, and columns named with a supported prefix but no data will be ignored by the function. Two IDs in the same string separated by a semicolon are separated into two entries.
#'
#' A complete and current list of currently supported prefixes for metabolites can be found by running the getPrefixesFromAnalytes() function, as shown in examples. These will include the following genes/proteins: ensembl, gene_symbol, uniprot, entrez, hmdb, wikidata, EN, ncbiprotein, brenda, chebi. They will also include the following metabolites: hmdb, chebi, chemspider, kegg, pubchem, CAS, wikidata, LIPIDMAPS, lipidbank, swisslipids, plantfa, kegg_glycan, rhea-comp, polymer
#'
#' The COMETS data format, which uses HMDB_ID and PUBCHEM as column names, is also supported.
#'
#' @param data_frame a data.frame object where the column names are prefixes that will be prepended to each identifier in that column. Please specify either this argument or file_path. Supported prefixes can be found with getPrefixesFromAnalytes()
#' @param file_path a string containing the file path to a csv that will be read and converted to RaMP data input format. Column names should be prefixes that will be prepended to each identifier in that column.
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @param ... additional arguments that will be passed to readr::read_csv() or readxl::read_excel()
#'
#' @returns A character() object including each identifier in the input prepended to a supported prefix in the format "prefix:identifier". This result is designed to be input for several RaMP functions.
#'
#' @examples
#'\dontrun{
#' # Retrieve supported ID source lists for metabolites and genes
#' RaMP::getPrefixesFromAnalytes(analyteType = 'metabolite')
#'
#' RaMP::getPrefixesFromAnalytes(analyteType = 'gene')
#'
#' # Example use with demo data
#' df <- data.frame(ensembl = c('ENSG00000135679', 'ENSG00000141510'), hmdb = c('HMDB0000064', NA), fake_ID = c('123', NA))
#'
#' dir <- system.file("extdata", package="RaMP", mustWork=TRUE)
#' data_frame <- file.path(dir,"ExampleInputFileName.csv")
#'
#' RaMPInput <- createRaMPInput(data_frame = df)
#'
#' getPathwayFromAnalyte(RaMPInput)
#'
#' # Example use with csv
#' RaMPInput <- createRaMPInput(file_path = 'createRaMPInput_test_data.csv')
#'
#' getPathwayFromAnalyte(RaMPInput)
#'}
#' @export
createRaMPInput <- function ( data_frame = NULL, file_path = NULL, db = RaMP(), ... ) {

  #Produce dataframe of user data from csv path if not provided initially. User must specify exactly one of these.
  if (!is.null(file_path) & is.null(data_frame)) {
    if (substr(file_path, (nchar(file_path) - 3), nchar(file_path)) == '.csv') {
      data_frame <- readr::read_csv(file_path, ...)
    }
    else if (substr(file_path, (nchar(file_path) - 4), nchar(file_path)) == '.xlsx') {
      data_frame <- readxl::read_excel(file_path, ...)
    }

  }

  else if (!is.null(file_path) & !is.null(data_frame)) {
    stop('You have entered values for both the file_path and data_frame parameters. Please only specify one of these.')
  }

  else if (is.null(file_path) & is.null(data_frame)) {
    stop('You have not entered any input data to convert. Please specify either the file_path or data_frame parameter.')
  }

  data_frame <- tibble::tibble(data_frame)

  #Get acceptable prefixes from RaMP and convert them to character vector.
  prefix_c <- strsplit(getPrefixesFromAnalytes(db = db, analyteType = 'metabolite')$idTypes[1], ',')[[1]]

  prefix_c <- append(prefix_c, strsplit(getPrefixesFromAnalytes(db = db, analyteType = 'gene')$idTypes[1], ',')[[1]])

  prefix_c <- gsub(' ', '', prefix_c)

  prefix_c <- unique(prefix_c)

  #Map from COMETS column names to RaMP prefixes
  COMETS_convert_df <- data.frame(COMETS = c('HMDB_ID', 'PUBCHEM'), RaMP = c('hmdb', 'pubchem'))

  for (i in colnames(data_frame)) {

    if (i %in% COMETS_convert_df$COMETS) {
      colnames_c <- colnames(data_frame)
      colnames_c[which(colnames_c == i)] <- COMETS_convert_df$RaMP[[which(COMETS_convert_df$COMETS == i)]]
      colnames(data_frame) <- colnames_c
    }
  }


  #Parse data_frame into RaMP prefix:ID input
  final_id_c <- c()

  for (i in colnames(data_frame)) {

    if (i %in% prefix_c) {
      ids <- dplyr::pull(data_frame[,i])
      ids <- as.character(ids)
      ids <- ids[ids != 'NA']
      ids <- ids[!is.na(ids)]
      ids <- unlist(strsplit(ids, split = ';'))

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
