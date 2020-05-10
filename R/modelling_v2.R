# modelling v2

#### modelling regression model 
# leave one out method 


# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
# packages to install 
#install.packages("devtools")
#devtools::install_github("tidyverse/broom")
#remotes::install_github("dariyasydykova/tidyroc")
library(tidyverse)
library(tidyroc)
library(broom)
library(yardstick)
library(cowplot)
# set the second level of truth to positive in yardstick
options(yardstick.event_first = FALSE)



#### decisions 
# Generalised linear models, e.g. stats::glm().
# Linear models assume that the response is continuous 
# and the error has a normal distribution. Generalised
# linear models extend linear models to include non-continuous
# responses (e.g. binary data or counts). They work by defining 
# a distance metric based on the statistical idea of likelihood.


# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")
# Wrangle data
# ------------------------------------------------------------------------------
# make yes no varrable as numeric and select unique peptides
data_single_peptides <- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(identifier, .keep_all = T) %>% 
  ungroup() %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el)) %>%
  mutate(label = case_when(response=="yes" ~ 1,
                           response=="no" ~ 0))

#### downsample negative data 
X_neg <- data_single_peptides %>% 
  filter(label==0) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,allele_frequency,priority_score,new_score,label) %>% 
  sample_n(40)
X_pos <- data_single_peptides %>% 
  filter(label==1) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,allele_frequency,priority_score,new_score,label)
# bindpos and negative dataset 
X_selected <- bind_rows(X_pos,X_neg) %>% 
  sample_n(70)




# Logistic regression with all chossen parameters 
glm_out1 <- glm(
  formula = label ~ mut_mhcrank_el+
  self_similarity+
  expression_level,
  family = binomial,
  data = X_selected) %>%
  augment() %>%
  mutate(model = "all") 

# logistic regression with self similarity
glm_out2 <- glm(label ~ self_similarity,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "self_sim") # name the model

# logistic regression with mutant rank el 
glm_out3 <- glm(label ~ mut_mhcrank_el,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "mut_mhcrank_el") # name the model

# logistic regression with expression score 
glm_out4 <- glm(label ~ expression_level,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "expression_level") # name the model


# combine the two datasets to make an ROC curve for each model
glm_out <- bind_rows(glm_out1, glm_out2,glm_out3,glm_out4)
glm_out$label <- as.factor(glm_out$label)
# plot ROC curves
glm_out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(truth = label, .fitted) %>% # get values to plot an ROC curve
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
 # scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() 
 






