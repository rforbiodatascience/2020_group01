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
#library(ROCR)
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
  mutate(label = case_when(response=="yes" ~ 1,
                           response=="no" ~ 0))
#### downsample negative data 
X_neg <- data_single_peptides %>% 
  filter(label==0) %>% select(mut_mhcrank_el,self_similarity,expression_level,label) %>% 
  sample_n(40)
X_pos <- data_single_peptides %>% 
  filter(label==1) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,label)
# bindpos and negative dataset 
X_selected <- bind_rows(X_pos,X_neg) %>% 
  sample_n(70)

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
Error_LogReg  = cbind(Error_LogReg,y) %>% as.tibble()
colnames(Error_LogReg) <- c(attributeNames,"True_val")                     
# For each crossvalidation fold
# make n for counting col possition
val <- tibble()
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

  Error_LogReg[k,n]  <- p                          
      # Error_LogReg[k,n]  <- (round(p,0)==y_test)
  }
}


Error_LogReg %>% gather(., key = "var", value = "T/F")

print(attributeNames,Error_LogReg[k+1,])

################### model merge all varriables 
Error_LogReg_merge = rep(NA, times=N+1)
for(k in 1:N){
    print(paste('Crossvalidation fold ', k, '/', N, sep=''));
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

library(tidyverse)
theme_set(theme_minimal())


roc <- Error_LogReg %>% gather(., key = variable, value  = error, - True_val) %>% 
  mutate(confs_var = 
           case_when(error > 0.5 & True_val==1 ~ "TP",
                     error > 0.5 & True_val==0 ~ "FP",
                      error < 0.5 & True_val==0 ~ "TN",
                      error < 0.5 & True_val==1 ~ "FN")) %>% 
  group_by(variable,confs_var) %>% 
  tally() %>% spread(., confs_var ,n) %>% 
  na.omit() %>% 
  mutate(TPR = TP/(TP+FN),
         FPR = FP/(FP+TN)) 



ggplot(roc , aes(x = FPR, y = TPR ))+
  geom_line(aes(color = variable  ))
  


library(tidyverse)
library(broom)

glm(am ~ disp,
    family = binomial,
    data = mtcars
) %>% # fit a binary classification model
  augment() %>% # get fitted values
  foo(predictor = .fitted, known_class = am) %>% # foo() is the function we want
  ggplot(aes(x = fpr, y = tpr)) + # plot an ROC curve
  geom_line()



roc <- Error_LogReg %>% as_tibble() %>% 
  gather(., key = variable, value  = error) %>% 
  group_by(variable,error) %>% 
  tally() %>% spread(., error,n)
  mutate(TPR = TP/(TP+FN),
         FPR = TN/(TN+FP)) %>% 
    

  
  
  roc %>% group_by(Metric) %>%
    summarise(AUC = sum(diff(FPR) * na.omit(lead(TPR) + TPR)) / 2)
  
  
  ggplot(roc, aes(FPR, TPR, color = Metric, frame = FPR)) +
    geom_line() +
    geom_abline(lty = 2)
  
  

roc  %>% tally()

roc <- iris %>% 
  gather(Metric, Value, -Species) %>%
  mutate(Positive = Species == "virginica") %>%
  group_by(Metric, Value) %>%
  summarise(Positive = sum(Positive),
            Negative = n() - sum(Positive)) %>%
  arrange(-Value) %>%
  mutate(TPR = cumsum(Positive) / sum(Positive),
         FPR = cumsum(Negative) / sum(Negative))

roc %>% group_by(Metric) %>%
  summarise(AUC = sum(diff(FPR) * na.omit(lead(TPR) + TPR)) / 2)


ggplot(roc, aes(FPR, TPR, color = Metric, frame = FPR)) +
  geom_line() +
  geom_abline(lty = 2)














library(tidyverse)
library(caret)      # for confusion matrix
library(ROCR)  
library(AUC)



devtools::install_github("drsimonj/pipelearner")
library(pipelearner)
pl <- pipelearner(X_selected, lm, label ~ mut_mhcrank_el)

#How cross validation is done is handled by learn_cvpairs(). For leave-one-out, specify k = number of rows:
  
  pl <- learn_cvpairs(pl, k = nrow(mtcars))
#Finally, learn() the model on all folds:
  
  pl <- learn(pl)
#This can all be written in a pipeline:
  
  pl <- pipelearner(mtcars, lm, hp ~ .) %>% 
  learn_cvpairs(k = nrow(mtcars)) %>% 
  learn()





r <- roc(churn$predictions,churn$labels)
 
# without leave one out 
(fmla <- as.formula(paste("y_train ~ ", paste(attributeNames, collapse= "+"))))

# Fit logistic regression model to predict the response
model = glm(label ~ mut_mhcrank_el, data = X_selected)
#  model = lm(fmla, data=X_train)
p = predict(model, data = X_selected,type="response")

p <- if_else(p > 0.5 , 1, 0)
confusionMatrix(data = p, 
                reference = X_selected$label, 
                positive = "1")






library(cvTools)

# Load data
library(R.matlab)
(fmla <- as.formula(paste("y ~ ", paste(attributeNames, collapse= "+"))))

# Fit logistic regression model to predict the type of wine
model = glm(fmla,family=binomial(link="logit"), data=X_selected)

# Evaluate the logistic regression on the test data
p = predict(model, newdata=X_selected, type="response")

## Plot receiver operating characteristic (ROC) curve
rocplot(p, y_test);
rocplot(p)


library(tidyverse)
library(broom)
library(tidyroc)

glm_out %>%
  group_by(model) %>% # group to get individual AUC values for each ROC curve
  make_roc(predictor = .fitted, known_class = outcome) %>%
  summarise(auc = calc_auc(x = fpr, y = tpr))
library(Deducer)
modelfit <- glm(label ~ mut_mhcrank_el, data = X_selected) 
rocplot(modelfit)









glm(am ~ disp, 
    family = binomial,
    data = mtcars
) %>%
  augment() %>%
  make_roc(predictor = fitted.values, known_class = y) %>%
  ggplot(aes(x = fpr, y = tpr)) + 
  geom_line()







p <- glm(label ~ mut_mhcrank_el, data = X_selected) %>% 
  predict(.,type=c("response"))
library(pROC)
roc(mut_mhcrank_el ~ p, data = X_selected)


mylogit <- glm(admit ~ gre, data = mydata, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
mydata$prob=prob
library(pROC)
g <- roc(admit ~ prob, data = mydata)
plot(g) 



  tidyroc::
  ggplot(aes(x = fpr, y = tpr)) + 
  geom_line()


glm(am ~ disp,
    family = binomial,
    data = mtcars
) %>% # fit a binary classification model
  augment() %>% # get fitted values
  roc_curve(truth = outcome, .fitted)  %>% # foo() is the function we want
  ggplot(aes(x = fpr, y = tpr)) + # plot an ROC curve
  geom_line()



mod2 <- lm(y ~ x, data = sim2)




X_selected


mod2 <- lm(label ~ self_similarity , data = X_selected)
mod3 <- lm(label ~ expression_level, data = X_selected)

mod1 <- glm(label ~ mut_mhcrank_el, data = X_selected)
grid <- X_selected %>% 
  data_grid(mut_mhcrank_el) %>% 
  add_predictions(mod1)



grid <- X_selected %>% 
#  data_grid(x = seq_range(x, n = 50, expand = 0.1)) %>% 
  gather_predictions(mod1, mod2, mod3, .pred = "label")

ggplot(X_selected, aes(x, y)) + 
  geom_point() +
  geom_line(data = grid, colour = "red") +
  facet_wrap(~ model)


roc <- Error_LogReg %>% as_tibble() %>% 
  gather(., key = variable, value  = error) %>% 
  group_by(variable,error) %>% 
  summarise(Positive = sum(error), 
            Negative = n() - sum(error)) %>% 
  mutate(TPR = cumsum(Positive) / sum(Positive), 
         FPR = cumsum(Negative) / sum(Negative))


roc %>% 
  group_by(variable) %>% 
  summarise(AUC = sum(diff(FPR) + na.omit(lead(TPR) + TPR)) / 2)

roc %>% 
  ggplot(aes( x= FPR, y = TPR, color = variable)) +
  geom_line() +
  geom_abline(lty = 2) +
  xlab("False positive rate (1-specificity)") + 
  ylab("True positive rate (sensitivity)") +
  ggtitle("ROC at predicting Virginica iris species")



roc <- iris %>% 
  gather(Metric, Value, -Species) %>% 
  mutate(Positive = Species == "virginica") %>% 
  group_by(Metric, Value) %>% 
  summarise(Positive = sum(Positive), 
            Negative = n() - sum(Positive)) %>% 
  arrange(-Value) %>% 
  mutate(TPR = cumsum(Positive) / sum(Positive), 
         FPR = cumsum(Negative) / sum(Negative))



################### model merge 


library(modelr)
library(tidyverse)
library(gapminder)
library(dlstats)    # for package download stats
library(pkgsearch)  # for searching packages






##### ROC Curves 

model_1 <- glm(fmla,family=binomial(link="logit"), data=X_train) %>% augment()
library(tidyverse)
library(broom)
library(tidyroc)
remotes::install_github("dariyasydykova/tidyroc")

model_1 %>% # group to get individual ROC curve for each model
  make_roc(truth = y_train, .fitted) %>% # get values to plot an ROC curve
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot()

library(tidyverse)
library(broom)
library(tidyroc)

model_1 %>%
  augment() %>%
  make_roc(predictor = .fitted, known_class = am) %>%
  ggplot(aes(x = fpr, y = tpr)) + 
  geom_line()
)


glm(outcome ~ clump_thickness,
    family = binomial,
    data = biopsy
) %>%
  augment() %>% # group to get individual ROC curve for each model
  roc_curve(truth = outcome, .fitted) %>% # get values to plot an ROC curve
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot()


glm(am ~ disp,
    family = binomial,
    data = mtcars
) %>% # fit a binary classification model
  augment() %>% # get fitted values
  foo(predictor = .fitted, known_class = am) %>% # foo() is the function we want
  ggplot(aes(x = fpr, y = tpr)) + # plot an ROC curve
  geom_line()



r <- roc(X_selected$expression_level, X_selected$label)

################### roc curves ####################
data_single_peptides <- data_single_peptides %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el))

X_neg <- data_single_peptides %>% 
  filter(label==0) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,label,allele_frequency,priority_score,new_score ) %>% 
  sample_n(40)
X_pos <- data_single_peptides %>% 
  filter(label==1) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,label,allele_frequency,priority_score,new_score ) 
# bindpos and negative dataset 
X_selected <- bind_rows(X_pos,X_neg) %>% 
  sample_n(70)



pdf(file = "Results/ROC_curve.png", width = 1000, height = 600)
par(mfrow=c(2,3)) 
pred <- prediction(data_single_peptides$allele_frequency, my_data_clean_aug)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "expression_level")

pred <- prediction(my_data_clean_aug$expression_level, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "expression_level")

pred <- prediction(my_data_clean_aug$mut_mhcrank_el, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "mut_mhcrank_el")

pred <- prediction(my_data_clean_aug$priority_score, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "priority_score")

pred <- prediction(my_data_clean_aug$new_score, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "expression_level/mut_mhcrank_el")

pred <- prediction(data_single_peptides$self_similarity, data_single_peptides$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "self_similarity")

dev.off()


###########################################3


















############# old code 







# linear model 
mod_diamond <- lm(mut_mhcrank_el ~ expression_score, data = my_data_clean_aug)



#### Trial linear model

options(na.action = na.warn)

mod_diamond <- lm(mut_mhcrank_el ~ expression_score, data = my_data_clean_aug)

my_data_clean_aug %>% 
  ggplot() +
  geom_point(aes(x = mut_mhcrank_el, y = expression_score))



my_data_clean_aug <- my_data_clean_aug %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el))

my_data_clean_aug %>% 
  ggplot(aes( x = priority_score )) + 
  geom_line(stat = "count") + 
  facet_grid(.~response) 

my_data_clean_aug  %>% 
  ggplot(aes( x = new_score, y = self_similarity )) + 
  geom_point()+
  facet_grid(.~response) +
  scale_x_log10() +  
  scale_y_log10() + 
  geom_smooth(method='lm')

## rocplot 

rocplot <- function(p, y){
  # ROCPLOT Plots the receiver operating characteristic (ROC)  curve and
  # calculates the area under the curve (AUC).
  # Notice: The method requires the package caTools to be installed!
  #
  # Usage:
  #   rocplot(p, y);
  #   res = rocplot(p, y);
  # 
  # Input: 
  #   p: Estimated probability of class 1. (Between 0 and 1.)
  #   y: True class indices. (Equal to 0 or 1.)
  #
  # Output:
  #    list containing:
  #   AUC: The area under the ROC curve
  #   TPR: True positive rate
  #   FPR: False positive rate
  #
  # Author: Laura Frølich, lff@imm.dtu.dk
  
  ## old code assuming values of p are distinct
  #res <- sort(p, decreasing=FALSE, index.return=TRUE);
  #val <- res$x
  #ind <- res$ix
  #x = y[ind];
  #FNR = cumsum(x==1)/sum(x==1);
  #TPR = 1-FNR;
  #TNR = cumsum(x==0)/sum(x==0);
  #FPR = 1-TNR;
  #TPR = c(1, TPR);
  #FPR = c(1, FPR);
  #AUC = t(-diff(FPR)) %*% (TPR[1:(length(TPR)-1)]+TPR[2:length(TPR)])/2;
  
  ## new code that does not require values of p to be distinct
  
  
  res <- sort(p, decreasing=FALSE, index.return=TRUE);
  val <- res$x
  ind <- res$ix
  x = y[ind];
  
  N0=sum(1-x);
  N1=sum(x);
  FNR=c(rep(0, times = length(x)), 1) # false negative rate
  TNR=c(rep(0, times = length(x)), 1) # true negative rate
  N_true=x[1];
  N_false=1-x[1];
  t=1;
  for(k in 2:length(val)){
    if(val[k-1]!=val[k]){
      t=t+1;
      FNR[t]=N_true/N1;
      TNR[t]=N_false/N0;
    }
    N_true=N_true+x[k];
    N_false=N_false+(1-x[k]);            
  }
  FNR[t+1]=1;
  FNR = FNR[1:(t+1)]
  TNR[t+1]=1;
  TNR = TNR[1:(t+1)]
  TPR = 1-FNR;
  FPR = 1-TNR;
  if(require(caTools)){
    AUC = -trapz(FPR,TPR);
  } else {
    error('rocplot.R: The package caTools is required to compute the AUC (Area Under Curve).')
  }
  
  plot(c(0, 1), c(0, 1), col='black', type='l', xlab='False positive rate (1-Specificity)', ylab='True positive rate (Sensitivity)', main='Receiver operating characteristic (ROC)', yaxt='n', xaxt='n')
  ticks <- seq(from=0, to=1, by=0.1)
  axis(1, at=ticks)
  axis(2, at=ticks)
  mtext(paste('AUC =', round(AUC, digits=3)))
  lines(FPR, TPR, col='red')
  grid(nx=length(ticks), ny=length(ticks), lwd=2)
  
  res <- list()
  res$AUC <- AUC
  res$TPR <- TPR
  res$FPR <- FPR
  
  res
}

