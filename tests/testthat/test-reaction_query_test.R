for (rampDB in test_databases) {
  test_that("test_mixed_analyte_reaction_query", {

    analytes = c(
      'chebi:58115',
      'chebi:456215',
      'chebi:58245',
      'chebi:58450',
      'chebi:17596',
      'chebi:16335',
      'chebi:16750',
      'chebi:172878',
      'chebi:62286',
      'chebi:77897',
      'uniprot:P30566',
      'uniprot:P30520',
      'uniprot:P00568',
      'uniprot:P23109',
      'uniprot:P22102',
      'uniprot:P15531'
    )

    # rampDB has been already been instantiated in helper-setup..
    rxns = RaMP::getReactionsForAnalytes(db = rampDB, analytes = analytes, humanProtein = T)
    expect_true(nrow(rxns$met2rxn) > 100)
    expect_true(nrow(rxns$prot2rxn) > 15)
    expect_true(nrow(rxns$metProteinCommonReactions) > 3)

  })


  test_that("test_metabolite_reaction_query", {

    mets = c(
      'chebi:58115',
      'chebi:456215',
      'chebi:58245',
      'chebi:58450',
      'chebi:17596',
      'chebi:16335',
      'chebi:16750',
      'chebi:172878',
      'chebi:62286',
      'chebi:77897'
    )

    # rampDB has been already been instantiated in helper-setup..
    rxns = RaMP::getReactionsForAnalytes(db = rampDB, analytes = mets, humanProtein = T)
    expect_true(nrow(rxns$met2rxn) > 100)
    expect_true(is.null(nrow(rxns$protein2rxn)))

  })

  test_that("test_protein_reaction_query", {

    proteins = c(
      'uniprot:P30566',
      'uniprot:P30520',
      'uniprot:P00568',
      'uniprot:P23109',
      'uniprot:P22102',
      'uniprot:P15531'
    )

    # rampDB has been already been instantiated in helper-setup..
    rxns = RaMP::getReactionsForAnalytes(db = rampDB, analytes = proteins, humanProtein = T)
    expect_true(nrow(rxns$prot2rxn) > 15)

  })


  test_that("test_mixed_analyte_reaction_class_query", {
    analytes = c(
      'chebi:58115',
      'chebi:456215',
      'chebi:58245',
      'chebi:58450',
      'chebi:17596',
      'chebi:16335',
      'chebi:16750',
      'chebi:172878',
      'chebi:62286',
      'chebi:77897',
      'uniprot:P30566',
      'uniprot:P30520',
      'uniprot:P00568',
      'uniprot:P23109',
      'uniprot:P22102',
      'uniprot:P15531'
    )

    rxnClasses <- RaMP:::getReactionClassesForAnalytes(db = rampDB, analytes = analytes)
    expect_true(nrow(rxnClasses$class_ec_level_1) > 4)
    expect_true(nrow(rxnClasses$class_ec_level_2) > 15)
    expect_true(nrow(rxnClasses$class_ec_level_3) > 25)

  })


  test_that("test_metabolite_reaction_class_query", {
    mets = c(
      'chebi:58115',
      'chebi:456215',
      'chebi:58245',
      'chebi:58450',
      'chebi:17596',
      'chebi:16335',
      'chebi:16750',
      'chebi:172878',
      'chebi:62286',
      'chebi:77897'
    )

    rxnClasses <- RaMP:::getReactionClassesForAnalytes(db = rampDB, analytes = mets)
    expect_true(nrow(rxnClasses$class_ec_level_1) > 4)
    expect_true(nrow(rxnClasses$class_ec_level_2) > 15)
    expect_true(nrow(rxnClasses$class_ec_level_3) > 20)

  })

  test_that("test_protein_reaction_class_query", {
    proteins = c(
      'uniprot:P30566',
      'uniprot:P30520',
      'uniprot:P00568',
      'uniprot:P23109',
      'uniprot:P22102',
      'uniprot:P15531'
    )
    rxnClasses <- RaMP:::getReactionClassesForAnalytes(db = rampDB, analytes = proteins)
    expect_true(nrow(rxnClasses$class_ec_level_1) > 2)
    expect_true(nrow(rxnClasses$class_ec_level_2) > 3)
    expect_true(nrow(rxnClasses$class_ec_level_3) > 4)

  })
}