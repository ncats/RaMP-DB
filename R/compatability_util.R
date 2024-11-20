getFunctionMetadata <- function(folder) {
  files <- list.files(folder, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  functionMetadata <- list()

  currentFile <- basename(sys.frame(1)$ofile)

  for (file in files) {
    if (grepl(currentFile, file)) {
      next # don't load recursively
    }

    tryCatch({
      blocks <- roxygen2::parse_file(file)
      for (block in blocks) {
        tags <- sapply(block$tags, function(x) x$tag)
        if ("export" %in% tags) {
          funcName <- block$object$alias
          functionMetadata[[funcName]] <- formals(block$object$value)
        }
      }
    }, error = function(e) {
      message("Error parsing file: ", file)
      message("Error message: ", e$message)
    })

  }
  return(functionMetadata)
}

getParameterInfo <- function(metaDataList, func) {
  sapply(names(metaDataList[[func]]), function(paramName) {
    paste0(paramName, "-", typeof(metaDataList[[func]][[paramName]]))
  })
}

findVersionChanges <- function(sourceFolder, referenceFolder) {
  sourceMetadata <- getFunctionMetadata(sourceFolder)
  referenceMetadata <- getFunctionMetadata(referenceFolder)

  sourceFunctionNames <- names(sourceMetadata)
  referenceFunctionNames <- names(referenceMetadata)

  newFunctionNames <- setdiff(sourceFunctionNames, referenceFunctionNames)
  removedFunctionNames <- setdiff(referenceFunctionNames, sourceFunctionNames)
  commonFunctions <- intersect(sourceFunctionNames, referenceFunctionNames)

  newFunctions <- list()
  removedFunctions <- list()
  functionsWithDifferentParams <- list()

  for (func in newFunctionNames) {
    params <- getParameterInfo(sourceMetadata, func)
    newFunctions[[func]] <- list(params = params)
  }

  for (func in removedFunctionNames) {
    params <- getParameterInfo(referenceMetadata, func)
    removedFunctions[[func]] <- list(params = params)
  }

  for (func in commonFunctions) {
    sourceParams <- getParameterInfo(sourceMetadata, func)
    referenceParams <- getParameterInfo(referenceMetadata, func)
    if (!identical(sourceParams, referenceParams)) {
      functionsWithDifferentParams[[func]] <- list(
        sourceParams = sourceParams,
        referenceParams = referenceParams
      )
    }
  }

  return (list(newFunctions = newFunctions, removedFunctions = removedFunctions, changedFunctions = functionsWithDifferentParams))
}

printVersionChanges <- function(sourceFolder, referenceFolder) {
  results <- findVersionChanges(sourceFolder = sourceFolder, referenceFolder = referenceFolder)

  cat('Removed Functions', '\n')
  for (func in names(results$removedFunctions)) {
    cat(func, ':')
    cat(paste(results$removedFunctions[[func]][['params']], collapse = ','), '\n')
  }

  cat('New Functions', '\n')
  for (func in names(results$newFunctions)) {
    cat(func, ':')
    cat(paste(results$newFunctions[[func]][['params']], collapse = ','), '\n')
  }

  cat('Changed Functions', '\n')
  for (func in names(results$changedFunctions)) {
    cat(func, '\n\t')
    cat(paste(results$changedFunctions[[func]][['sourceParams']], collapse = ','), '\n\t')
    cat(paste(results$changedFunctions[[func]][['referenceParams']], collapse = ','), '\n')
  }

}
# printVersionChanges("./R", "./../RaMP-DB-main/R")
