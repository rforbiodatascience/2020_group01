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

mupexi_ct26 <- read_xlsx(path = "data/_raw/ct26_library_mupexi.xlsx") %>% 
  # remove extra columns from previous handling
  select(-identifier, -Mut_peptide.y, -Allele) %>% 
  # convert Mut_MHCrank_EL and Expression level to numeric so we can join both files
  mutate(Mut_MHCrank_EL = as.numeric(Mut_MHCrank_EL),
         Expression_Level = as.numeric(Expression_Level))
mupexi_4t1 <- read_xlsx(path = "data/_raw/4T1_library_mupexi.xlsx") %>% 
  select(-identifier, -"...1")

# Wrangle data
# ------------------------------------------------------------------------------
mupexi_all <- full_join(mupexi_4t1, mupexi_ct26) %>% 
  # identifier column to merge with barracoda - HLA_peptidename
  mutate(identifier = paste(HLA_allele, Mut_peptide, sep = "_"))

my_data_clean <- my_data_clean %>% mutate(identifier = paste(HLA, Sequence, sep = "_")) %>% 
  mutate(HLA = str_replace(HLA, "^H-2", "H2-"))
  

my_data_clean_aug <- left_join(my_data_clean, mupexi_all, by = "identifier") %>% 
  group_by(response) %>% 
  # as a peptide could give (or not) a response in multiple samples, we will just keep one entry on each response group
  distinct(identifier, .keep_all = T) %>% 
  # select columns of interest to plot peptide characteristics
  select(c(,15:53)) %>% 
  # remove NAs (those are peptides that were controls for the barcoding experiment but not predicted by Mupexi)
  drop_na(HLA_allele)

h <- my_data_clean_aug %>% filter(duplicated(identifier))

duplicated(my_data_clean_aug$identifier)

# if a peptide gave a response in one sample but it did not in another, it is duplicated!

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          path = "data/03_my_data_clean_aug.tsv")