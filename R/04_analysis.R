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
  geom_text(my_data_clean_aug[my_data_clean_aug$response == "yes",], mapping = aes(label = Peptide.name))+
  facet_grid(vars(response))

# need to split axes until 1.5 to compare with estimated_frequency
my_data_clean_aug %>% 
  ggplot(., aes(Peptide.name, estimated_frequency_norm, colour=response)) +
  geom_point() +
  geom_text(my_data_clean_aug[my_data_clean_aug$response == "yes",], mapping = aes(label = Peptide.name))


######


my_data_clean_aug %>% filter(Mutation_Consequence=="M", str_length(Mut_peptide)==9) %>% 
  ggplot(aes(x=peptide_position))+
  geom_bar(aes(fill = response),stat = "count") +
  scale_y_log10()


# Write data
# ------------------------------------------------------------------------------
write_tsv(...)
ggsave(...)