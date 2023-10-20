test_that("Table Returned Is Not Null,getMetaFromOnto", {

  ontologies <- c("Colon", "Liver", "Lung")

  Metabolites <- getMetaFromOnto(db = rampDB, ontology = ontologies)

  expect_true(NROW(Metabolites) != 0)
  expect_true(!is.null(Metabolites))
})
