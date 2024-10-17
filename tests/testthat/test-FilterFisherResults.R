test_that("Fisher test results does not equal filtered fisher test results, FilterFisherResults",
          {

            analytes <-
                c(
                  "ensembl:ENSG00000135679",
                  "hmdb:HMDB0000064",
                  "hmdb:HMDB0000148",
                  "ensembl:ENSG00000141510"
                )

            fisher.results <-
              runCombinedFisherTest(db = rampDB, analytes = analytes)

            filtered.fisher.results <-
              FilterFishersResults(fisher.results, pValType='holm', pValCutoff  = 0.05)

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



