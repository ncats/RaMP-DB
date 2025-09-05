for (rampDB in test_databases) {
  test_that("getPathwayFromAnalyte,Returns table for analytes not NULL",
  {

    analytes <- c(
      'hmdb:HMDB0000056',
      'hmdb:HMDB0000439',
      'hmdb:HMDB0000479'
    )

    Pathways <- getPathwayFromAnalyte(db = rampDB, analytes = analytes)

    expect_true(
      !is.null(Pathways))
    expect_true(
      NROW(Pathways) != 0)
  })
}