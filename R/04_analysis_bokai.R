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
my_data_clean_aug %>% ...


#1.	


ggplot(data =my_data_clean_aug,
       mapping = aes(x=Mut_MHCrank_EL,y= Norm_MHCrank_EL) )+
  geom_point(aes(color=response))+
  labs( titile= "" x= "")
#2.

ggplot(data =my_data_clean_aug,
       mapping = aes(x=Mut_MHCrank_BA,y= Norm_MHCrank_BA) )+
  geom_point(aes(color=response))+
  labs( titile= "", x= "")


#3.
ggplot(data =my_data_clean_aug,
       mapping = aes(x=Mut_MHCrank_EL,y= Expression_Level) )+
  geom_point(aes(color=response))+
  labs( titile= "",x= "")

# Write data
# ------------------------------------------------------------------------------




