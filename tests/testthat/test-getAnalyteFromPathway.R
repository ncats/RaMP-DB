test_that("Table returned shows correct output for single pathway ,getAnalyteFromPathway",
          {
            library(properties)
            dbpass <- properties::read.properties('../../dbprops.txt')
            pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
            assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

            my_analytes <- getAnalyteFromPathway('sphingolipid metabolism')

            expect_true(
              !is.null(my_analytes)
              )
            expect_true(
              NROW(my_analytes) != 0)
            })



test_that("Table returned shows correct output for multiple pathways ,getAnalyteFromPathway",
       {
         library(properties)
         dbpass <- properties::read.properties('../../dbprops.txt')
         pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
         assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

        my_analytes <-
          getAnalyteFromPathway(c(
            'sphingolipid metabolism',
            "De Novo Triacylglycerol Biosynthesis",
            "Creatine metabolism"
          ))

        expect_true(
          !is.null(my_analytes)
        )
        expect_true(
          NROW(my_analytes) != 0)
      })






















