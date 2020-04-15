# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library(readxl)

# Define functions
# ------------------------------------------------------------------------------
file.choose()
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
path_ct26 <- "data/_raw/barracoda_output_CT26.xlsx"

all_ct26_barracoda_raw <- path_ct26 %>% 
  # function to import all sheets
  excel_sheets() %>% 
  # give names to each sheet
  set_names() %>% 
  # apply read_excel to each sheet, and add the number to the colum sheet
  map_df(~ read_excel(path = path_ct26, sheet = .x), .id = "sheet") 

path_4t1 <- "data/_raw/barracoda_output_4T1.xlsx"

all_4t1_barracoda_raw <- path_4t1 %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path_4t1, sheet = .x), .id = "scolheet") %>% 
  # remove introduced column "...17" because of weird excel format
  select(-"...17")

# Wrangle data
# ------------------------------------------------------------------------------
# all barracoda files together 
my_data <- full_join(all_ct26_barracoda_raw, all_4t1_barracoda_raw)

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data,
          path = "data/01_my_data.tsv")
