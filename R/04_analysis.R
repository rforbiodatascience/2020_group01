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
my_data_clean_aug<- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(., identifier, .keep_all = TRUE)

# Model data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% filter(cell_line=="CT26" & response == "yes") %>% 
  ggplot(., aes(peptide_name, estimated_frequency, colour=sample)) +
  geom_point() +
  geom_text(mapping = aes(label = peptide_name)) +
  facet_grid(vars(hla), scales = "free_y")

my_data_clean_aug %>% filter(cell_line=="CT26") %>% 
  ggplot(., aes(peptide_name, estimated_frequency_norm, colour=sample)) +
  geom_point() +
  geom_text(my_data_clean_aug %>% filter(cell_line=="CT26") %>% 
              filter(response == "yes"), mapping = aes(label = peptide_name)) +
  facet_grid(vars(hla), scales = "free_y")

my_data_clean_aug %>% filter(cell_line=="4T1") %>% 
  ggplot(., aes(peptide_name, estimated_frequency, colour=sample)) +
  geom_point() +
  geom_text(my_data_clean_aug %>% filter(cell_line=="4T1")  %>% 
              filter(response == "yes"), mapping = aes(label = peptide_name)) +
  facet_grid(vars(hla), scales = "free_y")


## new thing for bar plot missense mutatio 
my_data_clean_aug %>%
  filter(str_length(mut_peptide)==9,mutation_consequence=="M") %>% 
  ggplot(aes(x=peptide_position)) + 
  geom_bar(aes(fill = response), stat = "count")+
  scale_y_log10() + 
  theme_bw()

#### Trial linear model
library(modelr)
options(na.action = na.warn)

mod_diamond <- lm(mut_mhcrank_el ~ expression_score, data = my_data_clean_aug)

my_data_clean_aug %>% 
  ggplot() +
  geom_point(aes(x = mut_mhcrank_el, y = expression_score))


# Write data
# ------------------------------------------------------------------------------
write_tsv(...)
ggsave(...)