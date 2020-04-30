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
my_data <- read_tsv(file = "data/01_my_data.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
my_data_clean <- my_data %>% 
  # Remove "...17" column (useless - added from excel import); "HLA_allele" (same as HLA); identifier
  select(-`...17`, -HLA_allele, -identifier) %>% 
  # Remove control peptides that has introduced all NAs (present in barracoda file but not in Mupexi file)
  drop_na(Norm_peptide) %>% 
  # Round so that all numeric columns have 3 decimals
  mutate_if(is.numeric, round, 3)

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean,
          path = "data/02_my_data_clean.tsv")
