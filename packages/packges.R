install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)


setwd("packages")
create("neo_pack")
