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
  mutate(response = case_when(log_fold_change >= 2 & p_value <= 0.01 & 
                                input_1 & input_2 & input_3 != 0 & count > input_1  ~ "yes",
                              # everything else does not match any of the previous criterias as is labelled no
                              TRUE ~ "no")) %>% 
  # add organ column
  mutate(organ = case_when(str_detect(sample, "TU$") ~ "tumor",
                           str_detect(sample, "SP$") | str_detect(sample, "^CT26") ~ "spleen")) %>% 
  
  # remove _SP or _TU from sample
  mutate(sample = case_when(str_detect(sample, "TU$") ~ str_replace(sample, ".{3}$",  " ") ,
                           str_detect(sample, "SP$") ~ str_replace(sample, ".{3}$", " "),
                           TRUE ~ as.character(sample))) %>% 
  
  # add treatment column
  mutate(treatment = case_when(str_detect(sample, "^4T1_19") | str_detect(sample,"^4T1_23") | str_detect(sample,"^4T1_20") | 
                               str_detect(sample,"^4T1_16") | sample == "CT26_C1" | sample == "CT26_D1" | sample == "CT26_D2" ~ "CPI",
                               str_detect(sample,"^4T1_22") | str_detect(sample,"^4T1_17") | str_detect(sample,"^4T1_18") | 
                               sample == "CT26_C3" | sample == "CT26_C4" | sample == "CT26_D4" ~ "Isotype control")) %>% 
  
  # add cell_line column
  mutate(cell_line = case_when(str_detect(sample, "^4T1") ~ "4T1", 
                               str_detect(sample, "^CT26") ~ "CT26")) %>% 
  # add percent_count.fraction column = barcode_count / sum(barcode_counts) for a particular sample * 100
  group_by(sample) %>% 
  mutate(percent_count_fraction = count/sum(count)*100) %>% 
  # add estimated_frequency column
  mutate(estimated_frequency = percent_pe*percent_count_fraction/100) %>% 
  
  # add a sum of the normalized counts column for the significant 
  #epitopes (response == yes) to calculate the normalized estimated freq
  mutate(identifier = paste(neoepitope_sequence, hla, cell_line, sep = "_"))
     #    identifier = paste(identifier, cell_line, sep = "_"))


count_norm_signif <- my_data_clean_aug %>%
  group_by(sample) %>%
  filter(response == "yes") %>% 
  summarise(count_norm_signif = sum(count_normalised_edger))


my_data_clean_aug <- full_join(my_data_clean_aug, count_norm_signif) %>% 
  # add estimated_frequency_normalized column 
  mutate(estimated_frequency_norm = (count_norm_signif*percent_pe/count_norm_signif)) %>% 
  
  # change estimated_frequency_norm of non-response peptides to 0, as this measure is not relevant
  mutate(estimated_frequency_norm = case_when(response == "no" ~ 0,
                                              TRUE ~ as.numeric(estimated_frequency_norm)))


# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          path = "data/03_my_data_clean_aug.tsv")
