  test_that("chemicalClassSurvey returns not NULL for metabolite classes output", {

    library(properties)
    dbpass <- properties::read.properties('../../dbprops.txt')
    pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
    assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  metabolites.of.interest<- c('hmdb:HMDB0000056','hmdb:HMDB0000439','hmdb:HMDB0000479','hmdb:HMDB0000532',
                              'hmdb:HMDB0001015','hmdb:HMDB0001138','hmdb:HMDB0029159','hmdb:HMDB0029412',
                              'hmdb:HMDB0034365','hmdb:HMDB0035227','hmdb:HMDB0007973','hmdb:HMDB0008057',
                              'hmdb:HMDB0011211')

  chemical.classes<- chemicalClassSurvey(mets = metabolites.of.interest, background="NULL", background_type="database")

  metabolite.classes <- as.data.frame(chemical.classes$met_classes)

  expect_true(
   NROW(metabolite.classes)!=0
  )

  expect_true(
    !is.null(metabolite.classes)
  )
  })

  test_that("chemicalClassSurvey returns not NULL for count_summary output", {

    library(properties)
    dbpass <- properties::read.properties('../../dbprops.txt')
    pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
    assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

  metabolites.of.interest<- c('hmdb:HMDB0000056','hmdb:HMDB0000439','hmdb:HMDB0000479','hmdb:HMDB0000532',
                              'hmdb:HMDB0001015','hmdb:HMDB0001138','hmdb:HMDB0029159','hmdb:HMDB0029412',
                              'hmdb:HMDB0034365','hmdb:HMDB0035227','hmdb:HMDB0007973','hmdb:HMDB0008057',
                              'hmdb:HMDB0011211')

 chemical.classes<- chemicalClassSurvey(mets = metabolites.of.interest, background="NULL", background_type="database")
 count_summary<-as.data.frame(chemical.classes$count_summary$ClassyFire_sub_class)

  expect_true(
    NROW(count_summary)!=0
  )

  expect_true(
    !is.null(count_summary))
  })

