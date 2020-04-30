# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library(readxl)

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean <- read_tsv(file = "data/02_my_data_clean.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
my_data_clean_aug <- my_data_clean %>% 
  # add response column
  # yes response if logfold & pvalue conditions match & there is input in all three triplicates 
  mutate(response = case_when(log_fold_change >= 2 & p_value <= 0.01 & input_1&input_2&input_3 != 0 & count > input_1  ~ "yes",
                              # everything else does not match any of the previous criterias as is labelled no
                              TRUE ~ "no")) %>% 
  
  # add organ column
  mutate(organ = case_when(str_detect(sample, "TU$") ~ "tumor",
                           str_detect(sample, "SP$") | str_detect(sample, "^CT26") ~ "spleen")) %>% 
  
  # add treatment column
  mutate(treatment = case_when(str_detect(sample, "^4T1_19") | str_detect(sample,"^4T1_23") | str_detect(sample,"^4T1_20") | 
                                 str_detect(sample,"^4T1_16") | sample == "CT26_C1" | sample == "CT26_D1" | sample == "CT26_D2" ~ "yes",
                               str_detect(sample,"^4T1_22") | str_detect(sample,"^4T1_17") | str_detect(sample,"^4T1_18") | 
                                 sample == "CT26_C3" | sample == "CT26_C4" | sample == "CT26_D4" ~ "no")) %>% 
  # add cell_line column
  mutate(cell_line = case_when(str_detect(sample, "^4T1") ~ "4T1", str_detect(sample, "^CT26") ~ "CT26")) %>% 
  
  # add percent_count.fraction column = barcode_count / sum(barcode_counts) for a particular sample * 100
  group_by(sample) %>% 
  mutate(percent_count_fraction = round(count/sum(count)*100, 3)) %>% 
  
  # add estimated_frequency column
  mutate(estimated_frequency = percent_pe*percent_count_fraction/100) %>% 
  
  # add a sum of the normalized counts column for the significant epitopes (response == yes) to calculate the normalized estimated freq
  mutate(count_norm_signif = case_when(response == "yes" ~ sum(count_normalised_edger))) %>% 
  
  # add estimated_frequency_normalized column 
  mutate(estimated_frequency_norm =(count_norm_signif*percent_pe/count_norm_signif)) %>% 
  
  # add identifier column
  mutate(identifier = paste(neoepitope_sequence, hla, sep = "_"),
         identifier = paste(identifier, cell_line, sep = "_"))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          path = "data/03_my_data_clean_aug.tsv")

