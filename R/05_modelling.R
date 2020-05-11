#### modelling regression model 
# leave one out method 



# Clear workspace ---------------------------------------------------------
rm(list = ls())


# Load libraries ----------------------------------------------------------

# libraries for modelling 
# packages to install 
#install.packages("devtools")
#devtools::install_github("tidyverse/broom")
#remotes::install_github("dariyasydykova/tidyroc")
library(tidyverse)
library(tidyroc)
library(broom)
library(yardstick)
library(cowplot)

# Load data ---------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")


# # Wrangle data ----------------------------------------------------------
# make yes no varrable as numeric and select unique peptides
data_single_peptides <- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(identifier, .keep_all = T) %>% 
  ungroup() %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el)) %>%
  mutate(response_binary = case_when(response=="yes" ~ 1,
                           response=="no" ~ 0))



# Roc_curve ---------------------------------------------------------------
#### downsample negative data 
Roc_plot_list <- list()
for (i in 1:4) {
X_neg <- data_single_peptides %>% 
  filter(response_binary==0) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,
         allele_frequency,priority_score,new_score,response_binary) %>% 
  sample_n(40)
X_pos <- data_single_peptides %>% 
  filter(response_binary==1) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,
         allele_frequency,priority_score,new_score,response_binary)
# bind pos and negative dataset 
X_selected <- bind_rows(X_pos,X_neg) %>% 
  sample_n(70)


# Logistic regression with all chossen parameters 
glm_out1 <- glm(
  formula = response_binary ~ mut_mhcrank_el+
    self_similarity+
    expression_level,
  family = binomial,
  data = X_selected) %>%
  augment() %>%
  mutate(model = "all") 

# logistic regression with self similarity
glm_out2 <- glm(response_binary ~ self_similarity,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "self_similarity")

# logistic regression with mutant rank el 
glm_out3 <- glm(response_binary ~ mut_mhcrank_el,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "mut_mhcrank_el") 

# logistic regression with expression score 
glm_out4 <- glm(response_binary ~ expression_level,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "expression_level")


# combine all the data set 
glm_out <- bind_rows(glm_out1, glm_out2,glm_out3,glm_out4)
# make binary variable as factor 
glm_out <- glm_out %>% 
  mutate(response_binary = as.factor(response_binary))

# plot ROC curves
Roc_curves <- glm_out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(truth = response_binary, .fitted) %>% # get values to plot an ROC curve
  ggplot(
    aes(
      x = 1 - specificity,  # equal to FPR(False positive rate)
      y = sensitivity, # equal to TPR ( True Positive Rate)
      color = model
    )
  ) + # plot roc curves for each models 
  geom_line(size = 1.1) +
  labs( x = "FPR", y = "TPR") +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(labels = c("All", "Self similarity","Neoepitope EL (% rank)", "Expression Level"), 
                     values = c("#d7191c", "#fdae61","#abdda4","#2b83ba"), 
                     breaks = c("all", "self_similarity","mut_mhcrank_el","expression_level")) +
  coord_fixed() +
  theme_bw()
# reomve legedn on curves and save only legend for one curve 
if (i<4) {
  Roc_curves <- Roc_curves + theme(legend.position = "none")
}
if (i==4) {
  roc_legend <- get_legend(Roc_curves)
  Roc_curves <- Roc_curves + theme(legend.position = "none")
}
# add roc curve to plot list
Roc_plot_list[[i]] <- Roc_curves

AUC_values <- glm_out %>%
  group_by(model) %>% # group to get individual AUC value for each model
  roc_auc(truth = response_binary, .fitted)

}

# Combine plots
png(file = "Results/05_Roc_curves.png", width = 700, height = 500)
ggdraw() +
  draw_plot(Roc_plot_list[[1]], 0, .49, .5, .5) +
  draw_plot(Roc_plot_list[[2]], 0, .0, .5, .5) +
  draw_plot(Roc_plot_list[[3]], .38, .49, .5, .5) +
  draw_plot(Roc_plot_list[[4]], .38, .0, .5, .5) +
  draw_plot(roc_legend, .68, .33, .45, .4)
dev.off()






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

# Leave one out method 
# -----------------------------------------------------------------------------

# select predictive varriable 
y <- X_selected %>% 
  select(label) %>% 
  pull(.) # make as vector
# select  X data 
X_selected <- X_selected %>% 
  select(mut_mhcrank_el,self_similarity,expression_level)
# numer of rows

N <- X_selected %>% 
  nrow()
attributeNames <- colnames(X_selected ) %>% as.vector()
set.seed(1234) # for reproducibility
# subset dataset by random 
#classNames <- c('1','0')
# Initialize variables
tb <- as.tibble(c("true pos","true neg", "false pos", "false neg"))

Error_LogReg = matrix(rep(NA, times=N*length(attributeNames)), nrow=N)
#Error_LogReg  = cbind(Error_LogReg,y) %>% as.tibble()
colnames(Error_LogReg) <- c(attributeNames,"True_val")                     
# For each crossvalidation fold
# make n for counting col possition
#val <- tibble()
n=0
for (at in attributeNames) {
  X <- X_selected %>%  select(at)
  n =  n + 1
  for(k in 1:N){
    print(paste('Crossvalidation fold ', k, '/', N, sep=''));
  # Extract the training and test set
    X_train <- X[-k,]
    y_train <- y[-k] 
    X_test <-   X[k,]
    y_test <- y[k] 

  (fmla <- as.formula(paste("y_train ~ ", at)))
  # Fit logistic regression model to predict the response
  model = glm(fmla,family=binomial(link="logit"), data=X_train)
#  model = lm(fmla, data=X_train)
  p = model %>% predict(X_test,type="response") 

 # Error_LogReg[k,n]  <- p                          
  Error_LogReg[k,n]  <- (round(p,0)==y_test)
  }
}

Error_rate <-  (table(Error_LogReg)[1]/N)
Error_LogReg %>% gather(., key = "var", value = "T/F")

################### model merge all varriables 
Error_LogReg_merge = rep(NA, times=N)
for(k in 1:N){
    print(paste('Leave one out number: ', k, '/', N, sep=''));
    # Extract the training and test set
    X_train <- X_selected[-k,]
    y_train <- y[-k] 
    X_test <-   X_selected[k,]
    y_test <- y[k] 
    
    (fmla <- as.formula(paste("y_train ~ ", paste(attributeNames, collapse= "+"))))
   
     # Fit logistic regression model to predict the response
    model = glm(fmla,family=binomial(link="logit"), data=X_train)
  #  model = lm(fmla, data=X_train)
    p = predict(model, newdata=X_test,type="response")
   
     Error_LogReg_merge[k] = (round(p,0)==y_test)
}

Error_rate <-  (table(Error_LogReg_merge)[1]/N)
# 0.4 <- lm
# 0.39 <- glm  
# lm model 
Error_LogReg_merge %>% sum(FALSE)/length(Error_LogReg_merge)


###################### leave-one-out is done ####################


