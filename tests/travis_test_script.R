print(R.version)
.libPaths("/ramp-db/lib")
install.packages(c('filelock'))
install.packages('https://www.bioconductor.org/packages/3.18/bioc/src/contrib/Archive/BiocFileCache/BiocFileCache_2.10.1.tar.gz', repos=NULL, method='libcurl', dependencies=TRUE)
install.packages('RMariaDB')
install.packages('highcharter')
library(highcharter)
library(RMariaDB)
remotes::install_deps()

# Define a function to run tests with a given database type
run_tests <- function() {

  # Run tests
  test_results <- devtools::test()

  # Get indices where the condition is true
  indices <- sapply(1:length(test_results), function(i) {
    class_value <- class(test_results[[i]][[7]][[1]])[[1]]
    class_value == 'expectation_failure'
  })

  # Filter test_results based on the condition
  failures <- test_results[indices]

  # Check if any tests failed and exit with an error if true
  if (length(failures) > 0) {
    return(failures)
  } else {
    return(NULL)
  }
}

# Run tests with sqlite database
failures <- run_tests()

if (!is.null(failures)) {
    message("Some tests failed. Exiting with an error.")
    print(failures)
    q("no", status = 1)
}

message("All tests passed successfully.")
q("yes", status = 0)