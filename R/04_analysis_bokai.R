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
#my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...


view(my_data_clean_aug)
col(my_data_clean_aug)

colnames(my_data_clean_aug)
#scatterplot_function

scatterplot_function <- function(x,y) {  my_data_clean_aug %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_point(aes(color=response, alpha = response))+
    scale_y_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    scale_x_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))
  }


#1.

scatterplot_function('mut_mhcrank_el','norm_mhcrank_el')+
  labs(title= "Elution rank score of neoepitope vs WT epitope", x= "Neoepitope elution rank ", y="WT epitope elution rank")

#2.

scatterplot_function('mut_mhcrank_ba', 'norm_mhcrank_ba')+
  labs(title= "Barcode rank score of neoepitope vs WT epitope", x= "Neoepitope barcode rank ", y="WT epitope barcode rank")

#3.

scatterplot_function('mut_mhcrank_el', 'expression_level')+
  labs(title= "Elution rank score of neoepitope vs Expression level", x= "Neoepitope elution rank ", y="Expression level")





# create box_function for boxplot
box_function <- function(x,y) {  
  my_data_clean_aug %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_boxplot(aes(fill = response)) }

# boxplot of mutated elution vs response, facet by cell line
box_function('response','mut_mhcrank_el') +
  facet_wrap(cell_line~.) +
  labs(title ="Box plot of Mutated Elution vs Response (by Cell Line)", x = "Response", y = "Mutated Elution (Rank score)")

# boxplot of self-similarity vs response, facet by cell line
box_function('response','self_similarity') +
  facet_wrap(cell_line~.) +
  labs(title ="Box plot of Self Similarity vs Response (by Cell Line)", x = "Response", y = "Self Similarity")

# ------------------------------------------------------------------------------

treatment_names <- list(
  'no'="Treatment: No",
  'yes'="Treatment: Yes"
)

treatment_labeller <- function(variable,value){
  return(treatment_names[value])
}

# boxplot of mutated elution vs response, facet by treatment
box_function('response','mut_mhcrank_el') +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Mutated Elution vs Response (by Treatment)", x = "Response", y = "Mutated Elution (Rank score)")

# boxplot of self-similarity vs response, facet by treatment
box_function('response','self_similarity') +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Self Similarity vs Response (by Treatment)", x = "Response", y = "Self Similarity")

# ------------------------------------------------------------------------------

# create barplot_func for barplot
barplot_func <- function(num) {
  my_data_clean_aug %>%
    filter(str_length(mut_peptide)==num, mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() + 
    theme_bw()
  
}

# barplot of count vs peptide position, facet by cell line; str_length=10
barplot_func(10) +
  facet_wrap(cell_line~.) +
  labs(title ="Bar plot of Peptide Position (by Cell Line)", x = "Peptide Position", y = "Count")

# barplot of mutation consequence to illustrate FMID, facet by cell line
# FMID, i.e. F(frameshift), M(missense), I(inframe insertion), D(inframe deletion)
my_data_clean_aug %>% 
  ggplot(aes(x = mutation_consequence)) +
  geom_bar(aes(fill = response), stat = "count") +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(cell_line~.) +
  labs(title ="Bar plot of Mutation Consequence (by Cell Line)", x = "Mutation Consequence", y = "Count")








