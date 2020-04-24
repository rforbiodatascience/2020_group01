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
  mutate(response = case_when(log_fold_change >= 2 & p <= 0.01 & input.1&input.2&input.3 != 0 & count.1 > input.1  ~ "yes",
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
  # add percent_count.fraction column = barcode_count / sum(barcode_counts) for a particular sample * 100
  group_by(sample) %>% 
  mutate(percent_count.fraction = round(count.1/sum(count.1)*100, 3)) %>% 
  
  # add estimated_frequency column
  mutate(estimated_frequency = percent_PE*percent_count.fraction/100) 


my#########################################################################################
# To add the estimated_frequency_normalized column = `count.normalised (edgeR)`*percent_PE/sum(count.norm.signif),
# we need calculate the sum of the normalized counts for the peptides with significant responses. 

## I only know how to do it by getting the sum in a summarised dataset and merging it into the origianl dataframe

h <- my_data_clean_aug %>% group_by(sample) %>% 
  filter(response == "yes") %>% 
  summarise(count.norm.signif = sum(`count.normalised (edgeR)`))


my_data_clean_aug <- full_join(my_data_clean_aug, h) %>% 
  # add estimated_frequency_normalized column 
  mutate(estimated_frequency_norm =(`count.normalised (edgeR)`*percent_PE/count.norm.signif))

### summary tables to see how est.freq and est.freq.norm correlate
#########################################################################################
my_data_clean_aug %>% select(sample, estimated_frequency, estimated_frequency_norm, response, count.1, input.1, input.2, input.3, log_fold_change, p) %>% 
group_by(sample) %>% arrange(-estimated_frequency) %>% head(10)
my_data_clean_aug %>% select(sample, estimated_frequency, estimated_frequency_norm, response, count.1, input.1, input.2, input.3, log_fold_change, p) %>% 
group_by(sample) %>% arrange(-estimated_frequency_norm) %>% head(10)

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          path = "data/03_my_data_clean_aug.tsv")

