test_that("reaction enrichment data returns correctly formatted output", {

  analytes = c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450',
               'chebi:17596', 'chebi:16335', 'chebi:16750', 'chebi:172878',
               'chebi:62286', 'chebi:77897', 'uniprot:P30566','uniprot:P30520',
               'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')

  result <- runReactionClassTest(analytes = analytes, db = rampDB)
  expect_equal(
    (length(result)>1),
    TRUE
  )
  expect_true(
    NROW(result)!=0
  )

  expect_true(
    !is.null(result))
})

test_that("reaction class enrichment data is returned correctly when selecting for specific EC level, runReactionClassTest", {

  analytes = c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450',
               'chebi:17596', 'chebi:16335', 'chebi:16750', 'chebi:172878',
               'chebi:62286', 'chebi:77897', 'uniprot:P30566','uniprot:P30520',
               'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')

  reaction.classes<-getReactionClassesForAnalytes(db = rampDB, analytes = analytes)
  enrichedReactionSets <- runReactionClassTest(db = rampDB, analytes = analytes)
  EC_level1 <- enrichedReactionSets$EC_Level1Stats


  expect_true(
    NROW(EC_level1)!=0
  )

  expect_true(
    !is.null(EC_level1))
})


test_that("enrichment in chebi data returns correctly formatted output", {

  chebi = c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450',
               'chebi:17596', 'chebi:16335', 'chebi:16750', 'chebi:172878',
               'chebi:62286', 'chebi:77897')

  result <- runReactionClassTest(db = rampDB, analytes = chebi)
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


test_that("enrichment in uniprot data returns correctly formatted output", {

  uniprot = c('uniprot:P30566','uniprot:P30520',
               'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')

  result <- runReactionClassTest(db = rampDB, analytes = uniprot)
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
