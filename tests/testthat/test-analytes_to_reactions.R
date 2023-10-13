test_that("analytes_to_reactions", {
  library(properties)
  dbpass <- properties::read.properties('../../dbprops.txt')
  pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
  assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  analytes = c("uniprot:P04406", "uniprot:P00338", "uniprot:P11413", "uniprot:Q99798", "uniprot:P08559",
               "chebi:15361", "chebi:14314", "chebi:24996")

  rxnResult <- RaMP::getReactionsForAnalytes(analytes=analytes, analyteType="both", humanProtein = T)

  expect_true(ncol(rxnResult[[1]]) == 17, label="met2rxn 17 columns")

  expect_true(nrow(rxnResult[[1]]) > 0, label = "have met2rxn results")

  expect_true(ncol(rxnResult[[2]]) == 15, label="prot2rxn 15 columns")

  expect_true(nrow(rxnResult[[2]]) > 0, label = "have prot2rxn results")

  expect_true(ncol(rxnResult[[3]]) == 12, label="mp_common_rxn 12 columns")

  expect_true(nrow(rxnResult[[3]]) > 0, label = "have mp_common_rxn results")
})


test_that("metabolites_to_reactions", {
  library(properties)
  dbpass <- properties::read.properties('../../dbprops.txt')
  pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
  assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  analytes = c("chebi:15361", "chebi:14314", "chebi:24996")

  rxnResult <- RaMP::getReactionsForAnalytes(analytes=analytes, analyteType="metabolites", humanProtein = T)

  expect_true(ncol(rxnResult[[1]]) == 17, label="met2rxn 17 columns")

  expect_true(nrow(rxnResult[[1]]) > 0, label = "have met2rxn results")

  expect_true(ncol(rxnResult[[2]]) == 0, label="empty prot2rxn 15 columns")

  expect_true(nrow(rxnResult[[2]]) == 0, label = "empty have prot2rxn results")

  expect_true(ncol(rxnResult[[3]]) == 0, label="empty mp_common_rxn")

  expect_true(nrow(rxnResult[[3]]) == 0, label = "empty mp_common_rxn")
})


test_that("proteins_to_reactions", {
  library(properties)
  dbpass <- properties::read.properties('../../dbprops.txt')
  pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
  assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  analytes = c("uniprot:P04406", "uniprot:P00338", "uniprot:P11413", "uniprot:Q99798", "uniprot:P08559")

  rxnResult <- RaMP::getReactionsForAnalytes(analytes=analytes, analyteType="both", humanProtein = T)

  expect_true(ncol(rxnResult[[1]]) == 0, label="empty met2rxn")

  expect_true(nrow(rxnResult[[1]]) == 0, label = "empty met2rxn")

  expect_true(ncol(rxnResult[[2]]) == 15, label="prot2rxn 15 columns")

  expect_true(nrow(rxnResult[[2]]) > 0, label = "have prot2rxn results")

  expect_true(ncol(rxnResult[[3]]) == 0, label="empty mp_common_rxn")

  expect_true(nrow(rxnResult[[3]]) == 0, label = "empty mp_common_rxn")
})


