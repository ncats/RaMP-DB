for (rampDB in test_databases) {
  test_that("chemical class enrichment data is returned correctly, ChemicalClassEnrichment", {

    metabolites.of.interest <- c('hmdb:HMDB0000056',
                                 'hmdb:HMDB0000439',
                                 'hmdb:HMDB0000479',
                                 'hmdb:HMDB0000532',
                                 'hmdb:HMDB0011211')

    chemical.classes <- getChemClass(db = rampDB, mets = metabolites.of.interest, background = "NULL", backgroundType = "database")

    enrichedClassSets <- runEnrichChemClass(db = rampDB, mets = metabolites.of.interest, background = "NULL", backgroundType = "database")

    expect_true(
      NROW(enrichedClassSets) != 0
    )

    expect_true(
      !is.null(enrichedClassSets))
  })


  test_that("chemical class enrichment data is returned correctly, ChemicalClassEnrichment, where background is Saliva and inferIdMapping is FALSE", {

    metabolites.of.interest <- c('hmdb:HMDB0000056',
                                 'hmdb:HMDB0000439',
                                 'hmdb:HMDB0000479',
                                 'hmdb:HMDB0000532',
                                 'hmdb:HMDB0011211')

    chemical.classes <- getChemClass(db = rampDB, mets = metabolites.of.interest, background = "NULL", backgroundType = "database")

    enrichedClassSets <- runEnrichChemClass(db = rampDB, mets = metabolites.of.interest, background = "Saliva", backgroundType = "biospecimen", inferIdMapping = F)

    expect_true(
      NROW(enrichedClassSets) != 0
    )

    expect_true(
      !is.null(enrichedClassSets))
  })


  test_that("chemical class enrichment data is returned correctly when selecting for specific chemical classes, ChemicalClassEnrichment", {

    metabolites.of.interest <- c('hmdb:HMDB0000056',
                                 'hmdb:HMDB0000439',
                                 'hmdb:HMDB0034365',
                                 'hmdb:HMDB0035227',
                                 'hmdb:HMDB0008057',
                                 'hmdb:HMDB0011211')
    chemical.classes <- getChemClass(db = rampDB, mets = metabolites.of.interest, background = "NULL", backgroundType = "database")
    enrichedClassSets <- runEnrichChemClass(db = rampDB, mets = metabolites.of.interest, background = "NULL", backgroundType = "database")
    classy_fire_classes <- enrichedClassSets$ClassyFire_class


    expect_true(
      NROW(classy_fire_classes) != 0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })


  test_that("chemical class enrichment data is returned correctly when selecting biospecimen background for Urine metabolites", {

    metabolites.of.interest <- c('hmdb:HMDB0000056',
                                 'hmdb:HMDB0000439',
                                 'hmdb:HMDB0034365',
                                 'hmdb:HMDB0035227',
                                 'hmdb:HMDB0008057',
                                 'hmdb:HMDB0011211')

    chemical.classes <- getChemClass(db = rampDB, mets = metabolites.of.interest, background = "Urine", backgroundType = "biospecimen")
    enrichedClassSets <- runEnrichChemClass(db = rampDB, mets = metabolites.of.interest, background = "Urine", backgroundType = "biospecimen")
    classy_fire_classes <- enrichedClassSets$ClassyFire_class


    expect_true(
      NROW(classy_fire_classes) != 0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })


  test_that("chemical class enrichment data is returned correctly when performing enrichment used inferIdMapping = T", {

    metabolites.of.interest <- c('hmdb:HMDB0000056',
                                 'hmdb:HMDB0000439',
                                 'hmdb:HMDB0034365',
                                 'hmdb:HMDB0035227',
                                 'hmdb:HMDB0008057',
                                 'hmdb:HMDB0011211')
    chemical.classes <- getChemClass(db = rampDB, mets = metabolites.of.interest, backgroundType = "database", inferIdMapping = T)
    enrichedClassSets <- runEnrichChemClass(db = rampDB, mets = metabolites.of.interest, backgroundType = "database", inferIdMapping = T)
    classy_fire_classes <- enrichedClassSets$ClassyFire_class


    expect_true(
      NROW(classy_fire_classes) != 0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })


  test_that("chemical class enrichment data is returned correctly when performing enrichment used inferIdMapping = T, and biospecimen background, Urine", {

    metabolites.of.interest <- c('hmdb:HMDB0000056',
                                 'hmdb:HMDB0000439',
                                 'hmdb:HMDB0034365',
                                 'hmdb:HMDB0035227',
                                 'hmdb:HMDB0008057',
                                 'hmdb:HMDB0011211')

    chemical.classes <- getChemClass(db = rampDB, mets = metabolites.of.interest, backgroundType = "database", inferIdMapping = T)
    enrichedClassSets <- runEnrichChemClass(db = rampDB, mets = metabolites.of.interest, background = 'Urine', backgroundType = "biospecimen", inferIdMapping = T)
    classy_fire_classes <- enrichedClassSets$ClassyFire_class

    expect_true(
      NROW(classy_fire_classes) != 0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })
}