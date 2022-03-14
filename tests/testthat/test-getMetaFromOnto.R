test_that("Table Returned Is Not Null,getMetaFromOnto", {
  library(properties)
  dbpass <- properties::read.properties('../../dbprops.txt')
  pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
  assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  ontologies <- c("Colon", "Liver", "Lung")

  Metabolites <- getMetaFromOnto(ontology = ontologies)

  expect_true(NROW(Metabolites) != 0)
  expect_true(!is.null(Metabolites))
})
