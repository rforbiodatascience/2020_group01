#### modelling regression model 
# leave one out method 


# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library(tidyverse)  # for data manipulation
#library(dlstats)    # for package download stats
#library(pkgsearch)
library(ROCR)
#library(ggplotify) # to make as.ggplot 
#library(gridExtra)
#library(grid)
#library(cvTools)
#library(R.matlab)
#library(modelr)

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

print(attributeNames,Error_LogReg[k+1,])

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


# old code 


## ROC Curves 

png(file = "Results/ROC_curve.png", width = 1000, height = 600)
par(mfrow=c(2,3)) 
pred <- prediction(X_selected$allele_frequency, X_selected$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "allele_frequency")

pred <- prediction(X_selected$expression_level, X_selected$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "expression_level")

pred <- prediction(X_selected$mut_mhcrank_el, X_selected$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "mut_mhcrank_el")

pred <- prediction(X_selected$priority_score, X_selected$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "priority_score")

pred <- prediction(X_selected$new_score, X_selected$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "expression_level/mut_mhcrank_el")

pred <- prediction(X_selected$self_similarity, X_selected$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "self_similarity")

dev.off()

