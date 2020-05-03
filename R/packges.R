## R packages 
#install.packages("devtools")
library(devtools)
#install.packages("roxygen2")
library(roxygen2)
devtools::install_github("hadley/devtools")

#my_function <- function(x) { x*2}

devtools::create("mypackage") 
