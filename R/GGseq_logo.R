library(ggplot2)
#install.packages(ggseqlogo)
library(ggseqlogo)

data(ggseqlogo_sample)

my_data_clean_aug %>% filter(str_length(Mut_peptide)==9) %>% 
  ggseqlogo(., Mut_peptide )


seqs_dna %>% select(MA0001.1) %>% 
ggseqlogo(.)



ggplot(.) + geom_logo(my_data_clean_aug$Mut_peptide) + theme_logo()


geom_logo(my_data_clean_aug$Mut_peptide)

geom_logo(data = NULL, method = "bits", seq_type = "auto",
          namespace = NULL, font = "roboto_medium", stack_width = 0.95,
          rev_stack_order = F, col_scheme = "auto", low_col = "black",
          high_col = "yellow", na_col = "grey20", plot = T, ...)