# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
#source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
treatment_names <- list(
  'no'="Treatment: No",
  'yes'="Treatment: Yes"
)

treatment_labeller <- function(variable,value){
  return(treatment_names[value])
}

##Geom_boxplot self-similarity on y and x = response yes and no, facet by treatment 
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = response, y = Self_Similarity, fill = response)) +
  geom_boxplot() +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Self-Similarity vs Response", x = "Response", y = "Self-Similarity")

##Geom_boxplot Mut_rank_el on y and x = response yes and no, facet by treatment
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = response, y = Mut_MHCrank_EL)) +
  geom_boxplot(fill="steelblue") +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Mutation Elution (Rank score) vs Response", x = "Response", y = "Mutated Elution (Rank score)")

##Bar plot illustrating the frequency of mutant position with use of column peptide_position
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = peptide_position)) +
  geom_bar(fill="steelblue") +
  labs(title ="Bar plot of Peptide Position", x = "Peptide Position", y = "Count")

##Bar plot illustrating the number of F(frameshift), M(missense), I(inframe insertion),  
## D(inframe deletion). Mutations are under "Mutation_Consequence" column.
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = Mutation_Consequence)) +
  geom_bar(fill="steelblue") +
  labs(title ="Bar plot of Mutation Consequence", x = "Mutation Consequence", y = "Count")

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)