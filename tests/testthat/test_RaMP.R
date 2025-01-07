for (rampDB in test_databases) {
  test_that("listRaMPVersions works", {
    res <- listRaMPVersions(local = TRUE)
    if (length(res))
      expect_true(is(res, "character"))
    res <- listRaMPVersions()
    expect_true(is(res, "character"))
    expect_true(length(res) > 0)
  })

  test_that(".version_from_db_file works", {
    ch <- "RaMP-SQLite_v4.5.2.sqlite.gz"
    res <- .version_from_db_file(ch)
    expect_equal(res, "4.5.2")

    ch <- "some_other_string"
    res <- .version_from_db_file(ch)
    expect_equal(res, "")
  })

  test_that(".valid_ramp_database works", {
    expect_error(.valid_ramp_database("error", TRUE), "database connection")
    dbc <- dbConnect(SQLite(), tempfile())
    expect_error(.valid_ramp_database(dbc, TRUE), "lacks required database")
  })

  test_that(".RaMP and accessors work", {
    res <- .RaMP(dbname = rampDB@dbname)
    expect_s4_class(res, "RaMP")

    expect_output(show(res), "RaMP")

    con <- .dbcon(res)
    expect_s4_class(con, "SQLiteConnection")
    dbDisconnect(con)
  })

  test_that(".get_local_db_version_list works", {
    vList <- .get_local_db_version_list()
    testthat::expect_true(length(vList) > 0)
  })


  test_that("get_remote_db_version_list works", {
    vList <- get_remote_db_version_list()
    testthat::expect_true(length(vList) > 0)
  })
}
