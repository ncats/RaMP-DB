test_that("enrichment in mixed data returns correctly formatted output", {
  pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2",
    username = "root",
    conpass = "mysql123",
    host = "localhost",
    socket = paste0(
      "/lscratch/",
      Sys.getenv("SLURM_JOB_ID"),
      "/mysql/mysql.sock"
    )
  )
  pathwaydfids_mixed <- getPathwayFromAnalyte(c(
    "ensembl:ENSG00000135679",
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148",
    "ensembl:ENSG00000141510"
  ))
  backgrounddf_mixed <- getPathwayFromAnalyte(c(
      "ensembl:ENSG00000135679",
      "hmdb:HMDB0000064",
      "hmdb:HMDB0000148",
      "ensembl:ENSG00000141510",
      "hmdb:HMDB0075426",
      "hmdb:HMDB0085890"
  ))  
  result <- runCombinedFisherTest(pathwaydf = pathwaydfids_mixed)
  result_custom_bg  <- runCombinedFisherTest(pathwaydf = pathwaydfids_mixed,
                                             backgrounddf = backgrounddf_mixed)
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
    expect_equal(
    length(result_custom_bg),
    2
  )
  expect_equal(
    class(result_custom_bg[[1]]),
    "data.frame"
  )
  expect_equal(
    ncol(result_custom_bg[[1]]),
    14
  )
})



test_that("enrichment in metabolite data returns correctly formatted output", {
  pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2", username = "root", conpass = "mysql123",
    host = "localhost",
    socket = paste0(
      "/lscratch/",
      Sys.getenv("SLURM_JOB_ID"),
      "/mysql/mysql.sock"
    )
  )
  pathwaydfids_metabolites <- getPathwayFromAnalyte(c(
    "hmdb:HMDB0000064",
    "hmdb:HMDB0000148"
  ))
  backgrounddf_metabolites <- getPathwayFromAnalyte(c(
      "hmdb:HMDB0000064",
      "hmdb:HMDB0000148",
      "hmdb:HMDB0075426",
      "hmdb:HMDB0085890"
  ))  
  result <- runCombinedFisherTest(pathwaydf = pathwaydfids_metabolites)
  result_custom_bg  <- runCombinedFisherTest(pathwaydf = pathwaydfids_metabolites,
                                             backgrounddf = backgrounddf_metabolites)
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
    expect_equal(
    length(result_custom_bg),
    2
  )
  expect_equal(
    class(result_custom_bg[[1]]),
    "data.frame"
  )
  expect_equal(
    ncol(result_custom_bg[[1]]),
    10
  )
})


test_that("enrichment in gene data returns correctly formatted output", {
  pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2", username = "root", conpass = "mysql123",
    host = "localhost",
    socket = paste0(
      "/lscratch/",
      Sys.getenv("SLURM_JOB_ID"),
      "/mysql/mysql.sock"
    )
  )
  pathwaydfids_genes <- getPathwayFromAnalyte(c(
    "ensembl:ENSG00000135679",
    "ensembl:ENSG00000141510"
  ))
  backgrounddf_genes <- getPathwayFromAnalyte(c(
      "ensembl:ENSG00000135679",
      "ensembl:ENSG00000141510"
  ))  
  result <- runCombinedFisherTest(pathwaydf = pathwaydfids_genes)
  result_custom_bg  <- runCombinedFisherTest(pathwaydf = pathwaydfids_genes,
                                             backgrounddf = backgrounddf_genes)
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
    expect_equal(
    length(result_custom_bg),
    2
  )
  expect_equal(
    class(result_custom_bg[[1]]),
    "data.frame"
  )
  expect_equal(
    ncol(result_custom_bg[[1]]),
    10
  )
})
