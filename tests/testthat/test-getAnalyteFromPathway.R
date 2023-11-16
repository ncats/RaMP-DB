test_that("Table returned shows correct output for single pathway ,getAnalyteFromPathway",
          {

            my_analytes <- getAnalyteFromPathway(db = rampDB, pathway='sphingolipid metabolism')

            expect_true(
              !is.null(my_analytes)
            )
            expect_true(
              NROW(my_analytes) != 0)
          })



test_that("Table returned shows correct output for multiple pathways ,getAnalyteFromPathway",
          {

            my_analytes <-
              getAnalyteFromPathway(db = rampDB, pathway=c(
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





test_that("Fuzzy match test for TCA and Creatine",
          {

            my_analytes <-
              getAnalyteFromPathway(db = rampDB, pathway=c(
                "TCA",
                "Creatine"
              ), match="fuzzy")

            print(dim(my_analytes))
            print(unique(my_analytes$pathwayName))

            expect_true(
              !is.null(my_analytes)
            )
            expect_true(
              NROW(my_analytes) != 0)
          })

