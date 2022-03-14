Sys.setenv("R_TESTS"="")
library(testthat)
library(RaMP)

pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2", username = "root", conpass = "",
    host = "localhost",
)

test_check("RaMP")

