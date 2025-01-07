# CREATE TEST DATA

getFullMatrixTestData <- function(cellValue = 0) {
  nameList <- c('one', 'two', 'three', 'four')
  fullSize <- length(nameList)
  fullMatrix <- data.frame(row.names = nameList, matrix(data = cellValue, ncol = fullSize, nrow = fullSize))
  colnames(fullMatrix) <- nameList
  diag(fullMatrix) <- 1
  return(fullMatrix)
}

getSparseArrayTestData <- function() {
  data <- c(
    c('1,-|1,200|1,300|1,400|1,500|1,600|1,700|1,800|1,900|1,1000'), # similarity equals the others row number
    c('2,-|2,400'),                                                  # skip three, similar to 4
    c('4,-|6,1000'),                                                 # skip six, similar to 10
    c('10,-|1,1000')                                                 # similar to something outside
  )
  test_dataFrame <- data.frame(data = matrix(0, ncol = 2, nrow = 4))
  colnames(test_dataFrame) <- c('pathwayRampId', 'blob')
  test_dataFrame[, 1] <- c('one', 'two', 'three', 'four')
  test_dataFrame[, 2] <- data
  return(test_dataFrame)
}

# TEST CASES START HERE

test_that("similarity matrix initialized", {
  testData <- getFullMatrixTestData()

  simMat <- SimilarityMatrix$new(names = colnames(testData))

  expect_equal(simMat$getSize(), c(4, 4))
  expect_equal(length(simMat$getNames()), 4)
  expect_equal(rownames(simMat$getMatrix()), colnames(testData))
})

test_that("full matrix filtered by requested names", {
  cellValue <- 0.5
  testData <- getFullMatrixTestData(cellValue = cellValue)
  nameSubset <- c('one', 'three')

  simMat <- SimilarityMatrix$new(names = nameSubset)

  expect_equal(simMat$getSize(), c(2, 2))

  simMat$initializeFromFullDataFrame(testData)
  theMatrix <- simMat$getMatrix()

  expect_equal(theMatrix[1, 1], 1)
  expect_equal(theMatrix[1, 2], cellValue)
  expect_equal(theMatrix[2, 1], cellValue)
  expect_equal(theMatrix[2, 2], 1)
  expect_equal(simMat$getNames(), nameSubset)
})

test_that("full matrix works with real data", {
  db254 <- RaMP(version = "2.5.4")
  pathwayNames <- c("RAMP_P_000000003", "RAMP_P_000000005")
  simMat <- SimilarityMatrix$new(names = pathwayNames)
  simMat$initializeFromFullDataFrame(db254@api$legacySimMatrices$analyte_result)
  theMatrix <- simMat$getMatrix()

  expect_equal(simMat$getSize(), c(2, 2))
  expect_equal(simMat$getNames(), pathwayNames)

  expect_equal(theMatrix[1, 1], 1)
  expect_equal(theMatrix[1, 2], .3125, tolerance = .0001) # these numbers might change if you change the database version
  expect_equal(theMatrix[2, 1], .3125, tolerance = .0001)
  expect_equal(theMatrix[2, 2], 1)
})

test_that("sparse matrix construction works", {
  testData <- getSparseArrayTestData()
  simMat <- SimilarityMatrix$new(names = testData$pathwayRampId)
  simMat$initializeFromSparseDataFrame(testData, compress = FALSE)
  expect_equal(simMat$getNames(), testData$pathwayRampId)
  theMatrix <- simMat$getMatrix()

  expect_equal(theMatrix[1, 1], 1)
  expect_equal(theMatrix[1, 2], .2)
  expect_equal(theMatrix[1, 3], .4)
  expect_equal(theMatrix[1, 4], 1)
  expect_equal(theMatrix[2, 2], 1)
  expect_equal(theMatrix[2, 3], .4)
  expect_equal(theMatrix[2, 4], 0)
  expect_equal(theMatrix[3, 3], 1)
  expect_equal(theMatrix[3, 4], 1)
  expect_equal(theMatrix[4, 4], 1)
})

test_that('sparse matrix works with real data', {
  db300 <- RaMP(branch = 'ramp3.0')
  pathwayNames <- c("RAMP_P_000000003", "RAMP_P_000000005")
  simMat <- SimilarityMatrix$new(names = pathwayNames)
  rawBlobs <- db300@api$getPathwayOverlapBlobs(pathwayNames, "both")
  simMat$initializeFromSparseDataFrame(rawBlobs)
  theMatrix <- simMat$getMatrix()

  expect_equal(simMat$getSize(), c(2, 2))
  expect_equal(simMat$getNames(), pathwayNames)

  expect_equal(theMatrix[1, 1], 1)
  expect_equal(theMatrix[1, 2], .312, tolerance = .001) # these numbers might change if you change the database version
  expect_equal(theMatrix[2, 1], .312, tolerance = .001)
  expect_equal(theMatrix[2, 2], 1)
})
