test_that("enrichment in mixed data returns correctly formatted output", {

  analytes <- c(
    "ensembl:ENSG00000135679",
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148",
    "ensembl:ENSG00000141510"
  )
  result <- runCombinedFisherTest(db = rampDB, analytes = analytes)
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
  result <- runCombinedFisherTest(db = rampDB, analytes = analytes)
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

  genes <- c(
    "ensembl:ENSG00000135679",
    "ensembl:ENSG00000141510"
  )

  result <- runCombinedFisherTest(db = rampDB, analytes = genes)
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
