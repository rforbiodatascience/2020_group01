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
library(cowplot)

# Define functions
# ------------------------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv") 


# Wrangle data
# ------------------------------------------------------------------------------

# Model data
# ------------------------------------------------------------------------------
my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
data_peptides <- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(identifier, .keep_all = T)


my_data_clean_aug %>% filter(cell_line == "CT26") %>% 
  ggplot(., aes(mut_mhcrank_el, expression_level)) +
  geom_point(aes(colour = response, alpha = response, size = estimated_frequency_norm))+
  facet_grid(vars(treatment)) +
  geom_text_repel(my_data_clean_aug %>% 
            filter(response == "yes", cell_line=="CT26"), 
            mapping = aes(label = peptide_name))+
  theme_bw()


## new thing for bar plot missense mutatio 
bar_plot_func <- function(num) {
  my_data_clean_aug %>%
    filter(str_length(mut_peptide)==num,mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() + 
    theme_bw()
  
}

bar_plot_func(9) +
  facet_grid(vars(cell_line))




# No. mutations by HLA faceted by cell_line. This plot was made by Sara because of curiosity.
# It shows in what MHC allele do the response peptides bind. Most of the responsive peptides
# bind H2-Kd allele in both cell lines 
data_single_peptides %>% 
  filter(response == "yes") %>% 
  ggplot(aes(x = hla)) +
  geom_bar(aes(fill = hla), stat = "count") +
  facet_grid(vars(cell_line)) +
  theme_bw()


## new thing for bar plot missense mutation 

my_data_clean_aug %>%
  filter(str_length(mut_peptide)==9,mutation_consequence=="M") %>% 
  ggplot(aes(x=peptide_position)) + 
  geom_bar(aes(fill = response), stat = "count")+
  scale_y_log10() + 
  theme_bw()

#### Trial linear model
library(modelr)
options(na.action = na.warn)

mod_diamond <- lm(mut_mhcrank_el ~ expression_score, data = my_data_clean_aug)

my_data_clean_aug %>% 
  ggplot() +
  geom_point(aes(x = mut_mhcrank_el, y = expression_score))


# I make some nice code :) 

my_data_clean_aug %>%
  filter(str_length(mut_peptide)==10,mutation_consequence=="M") %>% 
  ggplot(aes(x=peptide_position)) + 
  geom_bar(aes(fill = response), stat = "count")+
  scale_y_log10() + 
  theme_bw()

################## GGseq logo 
# ------------------------------------------------------------------------------------------
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

# plot together 
pdf(file = "Results/GGseq_plot.pdf", width = 12, height = 6)
ggdraw() +
  draw_plot(plot_List[[1]], 0, .51, .45, .4) +
  draw_plot(plot_List[[2]], 0, .0, .45, .4) +
  draw_plot(plot_List[[3]], .47, .51, .45, .4) +
  draw_plot(plot_List[[4]], .47, .0, .45, .4) +
  draw_plot_label(c("Responses", "No responses"), c(-0.05,  -0.06), c(1,  .48), size = 22)

dev.off()



# Write data
# ------------------------------------------------------------------------------
write_tsv(...)
ggsave(...)