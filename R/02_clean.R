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
  # add response column
  mutate(response = case_when(
          # yes response if logfold & pvalue conditions match & there is input in all three triplicates 
          log_fold_change >= 2 & p <= 0.05 & input.1&input.2&input.3 != 0 & count.1 > input.1  ~ "yes",
          # everything else does not match any of the previous criterias as is labelled no
          TRUE ~ "no")) %>% 
  # add organ column
  mutate(organ = case_when(
    str_detect(sample, "TU$") ~ "tumor",
    str_detect(sample, "SP$") | str_detect(sample, "^CT26") ~ "spleen")) %>% 
  # add treatment column
  mutate(treatment = case_when(str_detect(sample, "^4T1_19") | str_detect(sample,"^4T1_23") | str_detect(sample,"^4T1_20") | str_detect(sample,"^4T1_16") | sample == "CT26_C1" | sample == "CT26_D1" | sample == "CT26_D2" ~ "yes",
    str_detect(sample,"^4T1_22") | str_detect(sample,"^4T1_17") | str_detect(sample,"^4T1_18") | sample == "CT26_C3" | sample == "CT26_C4" | sample == "CT26_D4" ~ "no"))

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean,
          path = "data/02_my_data_clean.tsv")
