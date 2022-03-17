test_that("Fisher test results does not equal filtered fisher test results, FilterFisherResults",
          {
            library(properties)
            dbpass <- properties::read.properties('../../dbprops.txt')
            pkg.globals <- setConnectionToRaMP(host=dbpass$hostname, dbname=dbpass$dbname, username=dbpass$username, conpass=dbpass$conpass)
            assign("pkg.globals", pkg.globals, envir = .GlobalEnv)

            analytes <-
                c(
                  "ensembl:ENSG00000135679",
                  "hmdb:HMDB0000064",
                  "hmdb:HMDB0000148",
                  "ensembl:ENSG00000141510"
                )

            fisher.results <-
              runCombinedFisherTest(analytes = analytes)

            filtered.fisher.results <-
              FilterFishersResults(fisher.results, pval_type='holm', pval_cutoff  = 0.05)

             fisher.results <-
              fisher.results$fishresults[,
                c("pathwayName",
                  "Pval_combined_Holm"
                  )]
            filtered.fisher.results <-
              filtered.fisher.results$fishresults[,
                c("pathwayName",
                  "Pval_combined_Holm"
                  )]

            Filt.Test <-
              max(filtered.fisher.results[, 2]
                             )
            Fish.Test <-
              max(fisher.results[, 2]
                             )
            expect_true(
              Filt.Test != Fish.Test)
            expect_true(
              !is.null(filtered.fisher.results))
          })



