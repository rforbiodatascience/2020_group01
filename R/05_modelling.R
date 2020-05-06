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

Error_LogReg = matrix(rep(NA, times=N*length(attributeNames)), nrow=N+1)
colnames(Error_LogReg) <- attributeNames
# For each crossvalidation fold
# make n for counting col possition
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
  p = predict(model, newdata=X_test,type="response")
  # 
   Error_LogReg[k,n] = (round(p,0)==y_test)
  }
  Error_LogReg[k+1,n] <- (table(Error_LogReg[,n])[1]/N)
}

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

#glm model 







################### model merge 


#### downsample negative data 
X_neg <- data_single_peptides %>% filter(label==0) %>% select(mut_mhcrank_el,self_similarity,expression_level,label) %>% sample_n(40)
X_pos <- data_single_peptides %>% filter(label==1) %>% select(mut_mhcrank_el,self_similarity,expression_level,label)
X_selected <- rbind(X_pos,X_neg) 
X_all <- X_selected %>% 
  select(mut_mhcrank_el,self_similarity,expression_level)

y <- X_selected$label
N <- nrow(X_all)
attributeNames <- as.vector(colnames(X_all))
set.seed(1234) # for reproducibility
# subset dataset by random 
classNames <- c('1','0')
# Initialize variables
##
attributeNames <- colnames(X_all)
Error_LogReg_merge = rep(NA, times=N)
#Error_LogReg = rep(NA,times = K*length(attributeNames))
# For each crossvalidation fold
#n=1k
n=0
for(k in 1:N){
  print(paste('Crossvalidation fold ', k, '/', N, sep=''));
  # Extract the training and test set
  X_train <- X_all [-k,]#X[CV$subsets[CV$which!=k], ]
  y_train <- y[-k] #y[CV$subsets[CV$which!=k]];
  X_test <-   X_all [k,]# X[CV$subsets[CV$which==k], ]
  y_test <- y[k] #y[CV$subsets[CV$which==k]]
  (fmla <- as.formula(paste("y_train ~ ", paste(attributeNames, collapse= "+"))))
  # Fit logistic regression model to predict the response
  w_est = glm(fmla,family=binomial(link="logit"), data=X_train)
  p = predict(w_est, newdata=X_test,type="response")
  Error_LogReg_merge[k] = (round(p,0)==y_test)
}






















############# old code 

table(Error_LogReg_merge)[1]/70











my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")
my_data_clean_aug$label <-  ifelse(my_data_clean_aug$response=="yes", 1,0)
my_data_clean_aug <- my_data_clean_aug %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el))

pdf(file = "Results/ROC_curve.pdf", width = 10, height = 6)
par(mfrow=c(2,3)) 
pred <- prediction(my_data_clean_aug$allele_frequency, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "allele_freq")

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

pred <- prediction(my_data_clean_aug$self_similarity, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE, main = "self_similarity")
dev.off()




# linear model 
mod_diamond <- lm(mut_mhcrank_el ~ expression_score, data = my_data_clean_aug)



#### Trial linear model

options(na.action = na.warn)

mod_diamond <- lm(mut_mhcrank_el ~ expression_score, data = my_data_clean_aug)

my_data_clean_aug %>% 
  ggplot() +
  geom_point(aes(x = mut_mhcrank_el, y = expression_score))












errors <- data.frame(Error_LogReg*100)
colnames(errors) <- c('Logistic regression', 'KNN')
pdf(file="Comparison Logistic regression_KNN.pdf", width = 10, height = 10)
boxplot(errors, ylab="Error rate ()%")
dev.off()
mean(errors$`Logistic regression`)
sd(errors$`Logistic regression`)
mean(errors$KNN)
sd(errors$KNN)
z <- Error_LogReg-Error_KNN[,22];
zb <- mean(z);
nu <- K-1;
sig <-  sd(z-zb) / sqrt(K);
alpha <- 0.05;
zL <- zb + sig * qt(alpha/2, nu);
zH <- zb + sig * qt(1-alpha/2, nu);


#testresult <- t.test(Error_LogReg, Error_DecTree, paired = TRUE)
if(zL <= 0 && zH >= 0){
  print('Classifiers are NOT significantly different');  
}else{
  print('Classifiers are significantly different');
  
}

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

