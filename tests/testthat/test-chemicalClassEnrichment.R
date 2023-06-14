test_that("chemical class enrichment data is returned correctly, ChemicalClassEnrichment", {
  library(properties)
  dbpass <- properties::read.properties('../../dbprops.txt')
  pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
  assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  metabolites.of.interest<- c('hmdb:HMDB0000056',
                              'hmdb:HMDB0000439',
                              'hmdb:HMDB0000479',
                              'hmdb:HMDB0000532',
                              'hmdb:HMDB0011211')




  chemical.classes<-chemicalClassSurvey(mets = metabolites.of.interest, background="NULL", background_type="database")

  enrichedClassSets <- chemicalClassEnrichment(mets = metabolites.of.interest, background="NULL", background_type="database")

  expect_true(
    NROW(enrichedClassSets)!=0
  )

  expect_true(
    !is.null(enrichedClassSets))
  })


test_that("chemical class enrichment data is returned correctly, ChemicalClassEnrichment, where background is Saliva and inferIdMapping is FALSE", {
  library(properties)
  dbpass <- properties::read.properties('../../dbprops.txt')
  pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
  assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  metabolites.of.interest<- c('hmdb:HMDB0000056',
                              'hmdb:HMDB0000439',
                              'hmdb:HMDB0000479',
                              'hmdb:HMDB0000532',
                              'hmdb:HMDB0011211')

  chemical.classes<-chemicalClassSurvey(mets = metabolites.of.interest, background="NULL", background_type="database")

  enrichedClassSets <- chemicalClassEnrichment(mets = metabolites.of.interest, background="Saliva", background_type="biospecimen", inferIdMapping = F)

  expect_true(
    NROW(enrichedClassSets)!=0
  )

  expect_true(
    !is.null(enrichedClassSets))
})


  test_that("chemical class enrichment data is returned correctly when selecting for specific chemical classes, ChemicalClassEnrichment", {
    library(properties)
    dbpass <- properties::read.properties('../../dbprops.txt')
    pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
    assign("pkg.globals", pkg.globals, envir = .GlobalEnv)


  metabolites.of.interest<- c('hmdb:HMDB0000056',
                              'hmdb:HMDB0000439',
                              'hmdb:HMDB0034365',
                              'hmdb:HMDB0035227',
                              'hmdb:HMDB0008057',
                              'hmdb:HMDB0011211')
   chemical.classes<-chemicalClassSurvey(mets = metabolites.of.interest, background="NULL", background_type="database")
   enrichedClassSets <- chemicalClassEnrichment(mets = metabolites.of.interest, background="NULL", background_type="database")
  classy_fire_classes <- enrichedClassSets$ClassyFire_class


  expect_true(
   NROW(classy_fire_classes)!=0
  )

  expect_true(
   !is.null(classy_fire_classes))
  })




  test_that("chemical class enrichment data is returned correctly when selecting biospecimen background for Urine metabolites", {
    library(properties)
    dbpass <- properties::read.properties('../../dbprops.txt')
    pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
    assign("pkg.globals", pkg.globals, envir = .GlobalEnv)


    metabolites.of.interest<- c('hmdb:HMDB0000056',
                                'hmdb:HMDB0000439',
                                'hmdb:HMDB0034365',
                                'hmdb:HMDB0035227',
                                'hmdb:HMDB0008057',
                                'hmdb:HMDB0011211')
    chemical.classes<-chemicalClassSurvey(mets = metabolites.of.interest, background="Urine", background_type="biospecimen")
    enrichedClassSets <- chemicalClassEnrichment(mets = metabolites.of.interest, background="Urine", background_type="biospecimen")
    classy_fire_classes <- enrichedClassSets$ClassyFire_class


    expect_true(
      NROW(classy_fire_classes)!=0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })



  test_that("chemical class enrichment data is returned correctly when performing enrichment used inferIdMapping = T", {
    library(properties)
    dbpass <- properties::read.properties('../../dbprops.txt')
    pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
    assign("pkg.globals", pkg.globals, envir = .GlobalEnv)


    metabolites.of.interest<- c('hmdb:HMDB0000056',
                                'hmdb:HMDB0000439',
                                'hmdb:HMDB0034365',
                                'hmdb:HMDB0035227',
                                'hmdb:HMDB0008057',
                                'hmdb:HMDB0011211')
    chemical.classes<-chemicalClassSurvey(mets = metabolites.of.interest, background_type="database", inferIdMapping = T)
    enrichedClassSets <- chemicalClassEnrichment(mets = metabolites.of.interest, background_type="database", inferIdMapping = T)
    classy_fire_classes <- enrichedClassSets$ClassyFire_class


    expect_true(
      NROW(classy_fire_classes)!=0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })


  test_that("chemical class enrichment data is returned correctly when performing enrichment used inferIdMapping = T, and biospecimen background, Urine", {
    library(properties)
    dbpass <- properties::read.properties('../../dbprops.txt')
    pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
    assign("pkg.globals", pkg.globals, envir = .GlobalEnv)


    metabolites.of.interest<- c('hmdb:HMDB0000056',
                                'hmdb:HMDB0000439',
                                'hmdb:HMDB0034365',
                                'hmdb:HMDB0035227',
                                'hmdb:HMDB0008057',
                                'hmdb:HMDB0011211')
    chemical.classes<-chemicalClassSurvey(mets = metabolites.of.interest, background_type="database", inferIdMapping = T)
    enrichedClassSets <- chemicalClassEnrichment(mets = metabolites.of.interest, background = 'Urine', background_type="biospecimen", inferIdMapping = T)
    classy_fire_classes <- enrichedClassSets$ClassyFire_class

    expect_true(
      NROW(classy_fire_classes)!=0
    )

    expect_true(
      !is.null(classy_fire_classes))
  })
