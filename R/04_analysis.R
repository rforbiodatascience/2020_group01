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


# Model data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
# No. mutations by HLA faceted by cell_line. This plot was made by Sara because of curiosity.
# It shows in what MHC allele do the response peptides bind. Most of the responsive peptides
# bind H2-Kd allele in both cell lines 
data_single_peptides %>% 
  filter(response == "yes") %>% 
  ggplot(aes(x = hla)) +
  geom_bar(aes(fill = hla), stat = "count") +
  facet_grid(vars(cell_line)) +
  theme_bw()


## new thing for bar plot missense mutation 
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


# I make some nice code :) 

my_data_clean_aug %>%
  filter(str_length(mut_peptide)==10,mutation_consequence=="M") %>% 
  ggplot(aes(x=peptide_position)) + 
  geom_bar(aes(fill = response), stat = "count")+
  scale_y_log10() + 
  theme_bw()





# Write data
# ------------------------------------------------------------------------------
write_tsv(...)
ggsave(...)