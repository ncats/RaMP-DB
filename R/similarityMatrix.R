getSimMatrixFromSparseData <- function(sparseDataFrame, compressed = TRUE) {
  simMatrix <- SimilarityMatrix$new(names = sparseDataFrame$pathwayRampId)
  simMatrix$initializeFromSparseDataFrame(sparseDataFrame = sparseDataFrame, compressed = compressed)
  return (simMatrix$getMatrix())
}

getSimMatrixFromFullMatrix <- function(fullMatrix, pathwayRampIds) {
  simMatrix <- SimilarityMatrix$new(names = pathwayRampIds)
  simMatrix$initializeFromFullDataFrame(fullMatrix = fullMatrix)
  return (simMatrix$getMatrix())
}

#' @importFrom R6 R6Class
SimilarityMatrix <- R6::R6Class(
  "SimilarityMatrix",
  public = list(
    initialize = function(names) {
      private$names <- names
      private$simMatrix <- private$createEmptyMatrix(names = names)
    },

    initializeFromFullDataFrame = function(fullMatrix) {
      private$simMatrix <- fullMatrix[private$names, private$names]
    },

    initializeFromSparseDataFrame = function(sparseDataFrame, compressed = TRUE) {
      dataFrame <- private$processBlobColumn(unprocessedDF = sparseDataFrame, compressed = compressed)
      private$convertSparseDataToMatrix(dataFrame)
    },

    getSize = function() {
      return (dim(private$simMatrix))
    },
    getNames = function() {
      return (private$names)
    },
    getMatrix = function() {
      return (private$simMatrix)
    }
  ),
  private = list(
    simMatrix = NULL,
    names = NULL,
    createEmptyMatrix = function(names) {
      matrixSize <- length(names)
      simMatrix <- data.frame(row.names = names, matrix(data = 0, ncol = matrixSize, nrow=matrixSize))
      diag(simMatrix) <- 1
      colnames(simMatrix) <- names
      return (simMatrix)
    },
    processBlobColumn = function(unprocessedDF, compressed) {
      processedDF <- unprocessedDF  %>%
        dplyr::mutate(similarity_pairs = sapply(blob, function(x) {
          if (is.null(x)) {
            return (NA)
          }
          if (compressed) {
            return (memDecompress(from=x, type = 'gzip', asChar = T))
          }
          return (x)
        })) %>%
        dplyr::mutate(index = sapply(similarity_pairs, function(x) {
          firstCommaIndex <- regexpr(',', x)
          as.integer(substr(x, 1, firstCommaIndex - 1))
        })) %>%
        dplyr::select('pathwayRampId', 'index', 'similarity_pairs')
      return (processedDF)
    },
    convertSparseDataToMatrix = function(sparseDataFrame) {
      for (mapping_row_index in 1:nrow(sparseDataFrame)) {
          sim_pairs <- sparseDataFrame[mapping_row_index, 'similarity_pairs'][[1]]
          pairs <- strsplit(sim_pairs, '[|]')[[1]]
          column_pathway_index <- 0
          for (pair_index in 1:length(pairs)) {
            pair_pieces <- strsplit(pairs[pair_index], ',')[[1]]
            column_index_diff <- as.integer(pair_pieces[1])
            column_pathway_index <- column_pathway_index + column_index_diff
            list_index <- which(sparseDataFrame$index == column_pathway_index)
            if (length(list_index) > 0 && list_index > mapping_row_index) {
              similarity <- as.double(pair_pieces[2]) / 1000
              private$simMatrix[mapping_row_index, list_index] <- similarity
              private$simMatrix[list_index, mapping_row_index] <- similarity
            }
          }
        }
    }
  )
)
