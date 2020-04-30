# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...


view(my_data_clean_aug)

colnames(my_data_clean_aug)
#scatterplot_function

scatterplot_function <- function(x,y) {  my_data_clean_aug %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_point(aes(color=response, alpha = response))+
    scale_y_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    scale_x_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    labs(title= "Elution rank score of neoepitope vs WT epitope", x= "Neoepitope elution rank ", y="WT epitope elution rank")
  }


#1.

scatterplot_function('mut_mhcrank_el','norm_mhcrank_el')

#2.

scatterplot_function('mut_mhcrank_ba', 'norm_mhcrank_ba')

#3.

scatterplot_function('mut_mhcrank_el', 'expression_level')















