library(ggplot2)
#install.packages(ggseqlogo)
library(ggseqlogo)
library(tidyverse)
library(cowplot)
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv") 


# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# plot seq log, there only exist responses in length 9 to 10 so that is the only one illustarted 
plot_List <- list()
i=1
for (num in 9:10) {
  for(r in c("yes","no"))  {
  p <- seqloggo_generator(my_data_clean_aug, num, r) 
  print(p)
  plot_List[[i]] <- p
  i <- i+1
  }
}

pdf(file = "Results/GGseq_plot.pdf", width = 12, height = 6)
ggdraw() +
  draw_plot(plot_List[[1]], 0, .51, .45, .4) +
  draw_plot(plot_List[[2]], 0, .0, .45, .4) +
  draw_plot(plot_List[[3]], .47, .51, .45, .4) +
  draw_plot(plot_List[[4]], .47, .0, .45, .4) +
  draw_plot_label(c("Responses", "No responses"), c(-0.05,  -0.06), c(1,  .48), size = 22)

dev.off()


# trash 


eightmer_no <- seqloggo_generator(my_data_clean_aug, 8, "no", "M", "CT26") 
ninemer_yes <- seqloggo_generator(my_data_clean_aug, 11, "yes", "M", "CT26")
ninemer_no <- seqloggo_generator(my_data_clean_aug, 9, "no", "M", "CT26")

tenmer_yes <- seqloggo_generator(my_data_clean_aug, 10, "yes", "M", "CT26")
tenmer_no <- seqloggo_generator(my_data_clean_aug, 10, "no", "M", "CT26")

elevenmer_yes <- seqloggo_generator(my_data_clean_aug, 11, "yes", "M", "CT26")
elevenmer_no <- seqloggo_generator(my_data_clean_aug, 11, "no", "M", "CT26")






my_data_clean_aug %>% 
  filter(str_length(mut_peptide)==10,response=="yes",mutation_consequence=="M") %>% 
  select(mut_peptide) %>% 
  ggseqlogo()

my_data_clean_aug %>% 
  filter(str_length(mut_peptide)==9,response=="no",mutation_consequence=="M") %>% 
  select(mut_peptide) %>% 
  ggseqlogo()





my_data_clean_aug %>% 
  filter(str_length(mut_peptide)==9,response=="yes",mutation_consequence=="M") %>% 
  select(mut_peptide) %>% 
  ggseqlogo()

hola <- my_data_clean_aug %>% select(mut_peptide, cell_line, peptide_length, response)

h <- hola %>% group_split(cell_line, peptide_length, response, keep=T)

ggplot() + geom_logo(seqs_dna) + theme_logo() + 
  facet_wrap(~seq_group, ncol=4, scales='free_x') 

# Dont know how to strapolate the faceting to data with a different structure
ggplot() + geom_logo(h) + theme_logo() + 
  facet_wrap(~seq_group, scales='free_x') 



#### modelling 
my_data_clean_aug <- my_data_clean_aug %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el+self_similarity))
max(my_data_clean_aug$mut_mhcrank_el)
my_data_clean_aug <- my_data_clean_aug %>% 
  mutate(new_score  = norm_mhcrank_el/(mut_mhcrank_el))
my_data_clean_aug %>% 
  ggplot(aes( x = priority_score )) + 
  geom_line(stat = "count") + 
  facet_grid(.~response) 

TEMPT  %>% 
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

yes_df <- TEMPT %>% filter(response=="yes")
no_df <- TEMPT %>% filter(response=="no") %>% sample_n(50)

df <- rbind(yes_df,no_df )

TEMPT_without_dup <- subset(TEMPT, TEMPT$cell_line=="CT26")

table(my_data_clean_aug_moddeling$label)

my_data_clean_aug$label <-  ifelse(my_data_clean_aug$response=="yes", 1,0)
my_data_clean_aug_moddeling <- my_data_clean_aug %>% filter(new_score < 3 )
pred <- prediction(my_data_clean_aug$allele_frequency, my_data_clean_aug$label)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)




