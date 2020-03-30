# Script explanation
# This scripts aims to pre-processed the raw barracoda output by:
## - Adding a response == yes / no column.
## - Adding a treatment == yes / no column.

# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)

# Load data
# -----------------------------------------------------------------------------
## The raw barracoda output comprises one sheet per sample. We will import all sheets into a unique dataframe.
path_ct26 <- "data/raw/barracoda_output_CT26.xlsx"

all_ct26_barracoda_raw <- path_ct26 %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path_ct26, sheet = .x), .id = "sheet")

path_4t1 <- "data/raw/barracoda_output_4T1.xlsx"

all_4t1_barracoda_raw <- path_4t1 %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path_4t1, sheet = .x), .id = "sheet")

# Wrangle data
# ------------------------------------------------------------------------------
## Add response column based on cut-off: 
all_ct26_barracoda <- all_ct26_barracoda_raw %>% 
  mutate(response = case_when(
    # counts of a response are enriched compared to input (barcode added to sample)
    count.1 < input.1 ~ "count<input", 
    # yes response if logfold & pvalue conditions match & there is input in all three triplicates 
    log_fold_change >= 2 & p <= 0.05 & input.1&input.2&input.3 != 0   ~ "yes",
    # everything else does not match any of the previous criterias as is labelled no
    TRUE ~ "no")) %>% 
  mutate(treatment = case_when(
    sample == "CT26_C1" | sample == "CT26_D1" | sample == "CT26_D2" ~ "yes",
    sample == "CT26_C3" | sample == "CT26_C4" | sample == "CT26_D4" ~ "no"
    )) %>% 
  mutate(organ="spleen") %>% 
  mutate(identifier = paste(HLA, Peptide.name, sep="_"))

all_ct26_barracoda %>% 
  filter(response=="yes") %>% 
  distinct(Peptide.name, .keep_all = T) #19 unique peptides triggered responses 

all_4t1_barracoda <- all_4t1_barracoda_raw %>% 
  select(-`...17`) %>% 
  mutate(response = case_when(
    # counts of a response are enriched compared to input (barcode added to sample)
    count.1 < input.1 ~ "count<input", 
    # yes response if logfold & pvalue conditions match & there is input in all three triplicates 
    log_fold_change >= 2 & p <= 0.05 & input.1&input.2&input.3 != 0   ~ "yes",
    # everything else does not match any of the previous criterias as is labelled no
    TRUE ~ "no")) %>% 
  mutate(treatment = case_when(
    str_detect(sample, "^4T1_19") | str_detect(sample,"^4T1_23") | str_detect(sample,"^4T1_20") | str_detect(sample,"^4T1_16") ~ "yes",
    str_detect(sample,"^4T1_22") | str_detect(sample,"^4T1_17") | str_detect(sample,"^4T1_18") ~ "no",
    TRUE ~ as.character(sample))) %>% 
  mutate(organ = case_when(
    str_detect(sample, "TU$") ~ "tumor",
    str_detect(sample, "SP$") ~ "spleen"
    )) %>% 
  mutate(identifier = paste(HLA, Peptide.name, sep="_"))

all_4t1_barracoda %>% 
  filter(response=="yes") %>% 
  distinct(Peptide.name, .keep_all = T) #14 unique peptides triggered responses 

