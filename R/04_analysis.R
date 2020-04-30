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
my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% 
  ggplot(., aes(Peptide.name, estimated_frequency, colour=sample)) +
  geom_point() +
  geom_text(my_data_clean_aug %>% 
              filter(response == "yes"), mapping = aes(label = Peptide.name))
  #facet_grid(vars(response))

# need to split axes until 1.5 to compare with estimated_frequency
my_data_clean_aug %>% 
  ggplot(., aes(Peptide.name, estimated_frequency_norm, colour=response)) +
  geom_point() +
  geom_text(my_data_clean_aug %>% 
            filter(response == "yes"), mapping = aes(label = Peptide.name))

## new thing for bar plot missense mutatio 
my_data_clean_aug %>%
  filter(str_length(Mut_peptide)==9,Mutation_Consequence=="M") %>% 
  ggplot(aes(x=peptide_position)) + 
  geom_bar(aes(fill = response), stat = "count")+
  scale_y_log10() + 
  theme_bw()
  


# Write data
# ------------------------------------------------------------------------------
write_tsv(...)
ggsave(...)