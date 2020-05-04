# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggplot2)
#install.packages(ggseqlogo)
library(ggseqlogo)
library(ggseqlogo)
library(cowplot)

# Define functions
# ------------------------------------------------------------------------------
wd = getwd()
source(file = paste0(wd,"/R/99_project_functions.R"))

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_tsv(paste0(wd,file = "/data/03_my_data_clean_aug.tsv"))


# Wrangle data
# ------------------------------------------------------------------------------
# select unique peptides 
data_single_peptides <- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(identifier, .keep_all = T)

# Model data
# ------------------------------------------------------------------------------
# my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------

P2 <- my_data_clean_aug %>% filter(cell_line == "CT26") %>% 
  ggplot(., aes(mut_mhcrank_el, expression_level)) +
  geom_point(aes(colour = response, alpha = response, size = estimated_frequency_norm))+
  facet_grid(vars(treatment)) +
  geom_text_repel(my_data_clean_aug %>% 
            filter(response == "yes", cell_line=="CT26"), 
            mapping = aes(label = peptide_name))+
  theme_bw()


## Bar plot for mutaion possition only missense mutaions 

# all nine mer 
P3 <- bar_plot_func(9) +
  facet_grid(vars(cell_line))+
  labs(title ="Bar plot of Peptide Position (by Cell Line)", 
       x = "Peptide Position", 
       y = "Count")

# all 10 mer 
P4 <- bar_plot_func(10) +
  facet_grid(vars(cell_line))+
  labs(title ="Bar plot of Peptide Position (by Cell Line)",
       x = "Peptide Position", 
       y = "Count")



# No. mutations by HLA faceted by cell_line. This plot was made by Sara because of curiosity.
# It shows in what MHC allele do the response peptides bind. Most of the responsive peptides
# bind H2-Kd allele in both cell lines 
P5 <- data_single_peptides %>% 
  filter(response == "yes") %>% 
  ggplot(aes(x = hla)) +
  geom_bar(aes(fill = hla), stat = "count") +
  facet_grid(vars(cell_line)) +
  theme_bw()



#1.

P6 <-scatterplot_function('mut_mhcrank_el','norm_mhcrank_el')+
  labs(title= "Elution rank score of neoepitope vs WT epitope", 
       x= "Neoepitope elution rank ", y="WT epitope elution rank")

#2.

P7 <- scatterplot_function('mut_mhcrank_ba', 'norm_mhcrank_ba')+
  labs(title= "Barcode rank score of neoepitope vs WT epitope",
       x= "Neoepitope binding affinity rank ", y="WT epitope barcode rank")

#3.

P8 <- scatterplot_function('mut_mhcrank_el', 'expression_level')+
  labs(title= "Elution rank score of neoepitope vs Expression level", 
       x= "Neoepitope elution rank ", y="Expression level")



# boxplot of mutated elution vs response, facet by cell line
P9 <-  box_function('response','mut_mhcrank_el') +
  facet_wrap(cell_line~.) +
  labs(title ="Box plot of Mutated Elution vs Response (by Cell Line)", x = "Response", y = "Mutated Elution (Rank score)")


# boxplot of self-similarity vs response, facet by cell line
P10 <- box_function('response','self_similarity') +
  facet_wrap(cell_line~.) +
  labs(title ="Box plot of Self Similarity vs Response (by Cell Line)", 
       x = "Response", 
       y = "Self Similarity") 


treatment_names <- list(
  'no'="Treatment: No",
  'yes'="Treatment: Yes"
)

treatment_labeller <- function(variable,value){
  return(treatment_names[value])
}



# boxplot of mutated elution vs response, facet by treatment
P11 <- box_function('response','mut_mhcrank_el') +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Mutated Elution vs Response (by Treatment)", 
       x = "Response", 
       y = "Mutated Elution (Rank score)")

# boxplot of self-similarity vs response, facet by treatment
P12 <- box_function('response','self_similarity') +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Self Similarity vs Response (by Treatment)", 
       x = "Response", 
       y = "Self Similarity")


## missense mutations 
P13 <- my_data_clean_aug %>% 
  ggplot(aes(x = mutation_consequence)) +
  geom_bar(aes(fill = response), stat = "count") +
  scale_y_log10() +
  theme_bw() +
  facet_wrap(cell_line~.) +
  labs(title ="Bar plot of Mutation Consequence (by Cell Line)", 
       x = "Mutation Consequence", 
       y = "Count")



################## GGseq logo 
# # ------------------------------------------------------------------------------------------
# #plot seq log, there only exist responses in length 9 to 10 so that is the only one illustarted
# plot_List <- list()
# i=1
# for (num in 9:10) {
#   for(r in c("yes","no"))  {
#     p <- seqloggo_generator(my_data_clean_aug, num, r)
#     print(p)
#     plot_List[[i]] <- p
#     i <- i+1
#   }
# }
# 
# # plot together
# pdf(file = "Results/GGseq_plot.pdf", width = 12, height = 6)
# ggdraw() +
#    draw_plot(plot_List[[1]], 0, .51, .45, .4) +
#    draw_plot(plot_List[[2]], 0, .0, .45, .4) +
#    draw_plot(plot_List[[3]], .47, .51, .45, .4) +
#    draw_plot(plot_List[[4]], .47, .0, .45, .4) +
#    draw_plot_label(c("Responses", "No responses"), c(-0.05,  -0.06), c(1,  .48), size = 22)
# dev.off()
# 
# # exploring data with shiny app 
# # ------------------------------------------------------------------------------
# source(file = "R/shiny_app.R")
# Exploring_data_shiny(my_data_clean_aug)
# 

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)