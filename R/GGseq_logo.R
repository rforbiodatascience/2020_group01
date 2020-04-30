library(ggplot2)
#install.packages(ggseqlogo)
library(ggseqlogo)

my_data_clean_aug %>% 
  filter(str_length(Mut_peptide)==9,response=="yes",Mutation_Consequence=="M") %>% 
  select(Mut_peptide) %>% 
  ggseqlogo()

my_data_clean_aug %>% 
  filter(str_length(Mut_peptide)==9,response=="no",Mutation_Consequence=="M") %>% sample_n(20) %>% 
  select(Mut_peptide) %>% 
  ggseqlogo()
