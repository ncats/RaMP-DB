test_that("enrichment in mixed data returns correctly formatted output", {

  analytes <- c(
    "ensembl:ENSG00000135679",
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148",
    "ensembl:ENSG00000141510"
  )
  result <- runCombinedFisherTest(db = rampDB, analytes = analytes, background="NULL", background_type="database")
  expect_equal(
    (length(result)>1),
    TRUE
  )
  expect_equal(
    class(result[[1]]),
    "data.frame"
  )
  expect_equal(
    (ncol(result[[1]])>5),
    TRUE
  )
})



test_that("enrichment in metabolite data returns correctly formatted output", {

  analytes <- c(
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148"
  )

  pathwaydfids_metabolites <- getPathwayFromAnalyte(db = rampDB, analytes = analytes)
  result <- runCombinedFisherTest(db = rampDB, analytes = analytes, background="NULL", background_type="database")
  expect_equal(
    (length(result)>1),
    TRUE
  )
  expect_equal(
    class(result[[1]]),
    "data.frame"
  )
  expect_equal(
    (ncol(result[[1]])>5),
    TRUE
  )
})


test_that("enrichment in gene data returns correctly formatted output", {

    analytes <- c(
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148"
  )

  pathwaydfids_metabolites <- getPathwayFromAnalyte(db = rampDB, analytes = analytes)
  result <- runCombinedFisherTest(db = rampDB, analytes = analytes, background="NULL", background_type="database")
  expect_equal(
    (length(result)>1),
    TRUE
  )

  expect_equal(
    (ncol(result[[1]])>5),
    TRUE
  )
})



