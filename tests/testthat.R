Sys.setenv("R_TESTS"="")
library(testthat)
library(RaMP)

pkg.globals <- setConnectionToRaMP(
    dbname = "ramp2", username = "root", conpass = "mysql123",
    host = "localhost",
    socket = paste0(
        "/lscratch/",
        Sys.getenv("SLURM_JOB_ID"),
        "/mysql/mysql.sock"
    )
)

test_check("RaMP")
