print(R.version)
.libPaths("/ramp-db/lib")
install.packages(c('filelock'))
install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/BiocFileCache_2.10.1.tar.gz', repos=NULL, method='libcurl', dependencies=TRUE)
install.packages('RMariaDB')
install.packages('highcharter')
library(highcharter)
library(RMariaDB)
remotes::install_deps()

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
  message("Some tests failed using the sqlite database. Exiting with an error.")
  print(failures)
  q("no", status = 1)
}

Sys.setenv(MYSQL_TEST = "true")
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
  message("Some tests failed using the mysql database. Exiting with an error.")
  print(failures)
  q("no", status = 1)
} else {
  message("All tests passed successfully.")
  q("yes", status = 0)
}