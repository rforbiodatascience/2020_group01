
#### modelling 
my_data_clean_aug <- my_data_clean_aug %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el+self_similarity))

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

library(tidyverse)  # for data manipulation
library(dlstats)    # for package download stats
library(pkgsearch)
library(ROCR)

my_data_clean_aug$label <-  ifelse(my_data_clean_aug$response=="yes", 1,0)
pred <- prediction(my_data_clean_aug$allele_frequency, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
ROC_allele_freq <- plot(perf,colorize=TRUE)



