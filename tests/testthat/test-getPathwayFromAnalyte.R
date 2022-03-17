  test_that("getPathwayFromAnalyte,Returns table for analytes not NULL",
    {
      library(properties)
      dbpass <- properties::read.properties('../../dbprops.txt')
      pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
      assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

       analytes <- c(
        'hmdb:HMDB0000056',
        'hmdb:HMDB0000439',
        'hmdb:HMDB0000479'
        )

      Pathways <- getPathwayFromAnalyte(analytes)

      expect_true(
        !is.null(Pathways))
      expect_true(
        NROW(Pathways) != 0)
    })
