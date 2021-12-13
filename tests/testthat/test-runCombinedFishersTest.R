test_that("enrichment in mixed data returns correctly formatted output", {
  pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2",
    username = "root",
    conpass = "",
    host = "localhost",
  )
  pathwaydfids_mixed <- getPathwayFromAnalyte(c(
    "ensembl:ENSG00000135679",
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148",
    "ensembl:ENSG00000141510"
  ))
  result <- runCombinedFisherTest(pathwaydf = pathwaydfids_mixed)
  expect_equal(
    length(result),
    2
  )
  expect_equal(
    class(result[[1]]),
    "data.frame"
  )
  expect_equal(
    ncol(result[[1]]),
    14
  )
})



test_that("enrichment in metabolite data returns correctly formatted output", {
  pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2", username = "root", conpass = "",
    host = "localhost",
  )
  pathwaydfids_metabolites <- getPathwayFromAnalyte(c(
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148"
  ))
  result <- runCombinedFisherTest(pathwaydf = pathwaydfids_metabolites)
  expect_equal(
    length(result),
    2
  )
  expect_equal(
    class(result[[1]]),
    "data.frame"
  )
  expect_equal(
    ncol(result[[1]]),
    10
  )
})


test_that("enrichment in gene data returns correctly formatted output", {
  pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2", username = "root", conpass = "",
    host = "localhost",
  )
  pathwaydfids_genes <- getPathwayFromAnalyte(c(
    "ensembl:ENSG00000135679",
    "ensembl:ENSG00000141510"
  ))
  result <- runCombinedFisherTest(pathwaydf = pathwaydfids_genes)
  expect_equal(
    length(result),
    2
  )
  expect_equal(
      class(result[[1]]),
    "data.frame"
  )
  expect_equal(
      ncol(result[[1]]),
    10
  )
})
