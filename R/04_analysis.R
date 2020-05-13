
# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggplot2)
#install.packages(ggseqlogo)
library(ggseqlogo)
library(cowplot)
library(ggbeeswarm)


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")


# Wrangle data ------------------------------------------------------------
# arrange data so "no" is plotted first 
my_data_clean_aug <- my_data_clean_aug %>%  
  arrange(response)
# select unique peptides 
data_single_peptides <- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(identifier, .keep_all = T) %>% 
  ungroup()


# Visualise data ----------------------------------------------------------
# 1) Barracoda characteristics --------------------------------------------
p1_CT26 <- barc_resp(data  = my_data_clean_aug, 
                     mouce_cell_line =  "CT26") +
  labs(title = "log fold change of CT26-tumor cell line neoepitope screen")

p1_4T1 <- barc_resp(data  = my_data_clean_aug, 
                    mouce_cell_line =  "4T1") +
  labs(title = "log fold change of 4T1-tumor cell line neoepitope screen")

png(file = "Results/04_fig1_barracoda.png", width = 1500, height = 1000)
ggdraw() +
  draw_plot(p1_4T1, .05, .52, .95, .45) +
  draw_plot(p1_CT26, .05, .0, .95, .45) +
  #draw_plot(barracoda_legend, .95, .4, .2, .2) +
  draw_plot_label(c("A", "B"), c(.0, .0), c(.95, .45), size = 22)
dev.off()

# pvalue and log fold cahnge of baracoda 
p2 <- data_single_peptides %>% 
  ggplot(aes(p_value, log_fold_change)) +
  geom_point(aes(color=response, size  = estimated_frequency_norm)) +
  scale_color_manual(values = respond_cols) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0)) +
  scale_y_continuous(breaks = c(-10, -5, 0, 2, 5)) +
  theme_bw() +
  geom_hline(yintercept = 2, linetype = "dotted") +
  geom_vline(xintercept = 0.01, linetype = "dotted") +
  facet_grid(cell_line~.) +
  labs(x = "p-value",
       y = "log fold change",
       size = "Nomalized estimated frequency")
ggsave(p2, filename ="Results/04_p2_barracoda_selection.png", width = 12, height = 7)

# 2) Boxplots response vs non-response  ---------------------------------------------
# Elution rank ------------------------------------------------------------
p3 <-  box_function(data = data_single_peptides, 
                    x = 'response',
                    y = 'mut_mhcrank_el') +
  labs(title ="Eluted ligand %rank",
       x = "Response", 
       y = "Neoepitopes eluted ligand %rank")

# Binding affinity --------------------------------------------------------
p4 <-  box_function(data = data_single_peptides, 
                    x = 'response',
                    y= 'mut_mhcrank_ba') +
  labs(title ="Binding affinity %rank",
       x = "Response", 
       y = "Neoepitopes binding affinity %rank")

# Self-similarity ---------------------------------------------------------
p5 <- box_function(data = data_single_peptides, 
                   x = 'response',
                   y= 'self_similarity',
                   no_legend = FALSE) +
  labs(title = "Self-similarity ", 
       x = "Response", 
       y = "Self Similarity",
       color = "Response")
# get legend
mut_bar_legend <- get_legend(p5)
p5 <- p5 + theme(legend.position = "none")

# Combine plots 
png(file = "Results/04_fig2_el_ba_similarity.png", width = 1200, height = 400)
ggdraw() +
  draw_plot(p3, .0, .0, .3, .95) +
  draw_plot(p4, .3, .0, .3, .95) +
  draw_plot(p5, .6, .0, .3, .95) +
  draw_plot(mut_bar_legend, .85, .4, .2, .2) +
  draw_plot_label(c("A", "B", "C"), c(.0,  .3, .6), c(.97, .97, .97), size = 22)
dev.off()

# 3) Barplots mutation position -------------------------------------------
# 9mer --------------------------------------------------------------------
p6 <- bar_plot_func(data = data_single_peptides,
                      pep_length = 9)

# 10mer -------------------------------------------------------------------
p7 <- bar_plot_func(data = my_data_clean_aug,
                    pep_length = 10,
                    no_legend = FALSE)
# get legend
mut_bar_legend <- get_legend(p7)
p7 <- p7 + theme(legend.position = "none")

# Combine plots 
png(file = "Results/04_fig_3_mutation_position.png", width = 800, height = 400)
ggdraw() +
  draw_plot(p6, 0, .0, .45, .95) +
  draw_plot(p7, .45, .0, .45, .95) +
  draw_plot(mut_bar_legend, .85, .4, .2, .2) +
  draw_plot_label(c("A", "B"), c(0, 0.45), c(1, 1), size = 22)
dev.off()

# Neoepitope characteristics and est.freq --------------------------------
# WT vs Neo El and BA -----------------------------------------------------
p8 <-scatterplot_function(data = data_single_peptides,
                          x = 'mut_mhcrank_el',
                          y = 'norm_mhcrank_el',
                          no_legend = FALSE)+
  labs(title= "Eluted ligand of neoepitope vs WT epitope", 
       x= "Neoepitope eluted ligand %Rank ", 
       y="WT epitope eluted ligand %Rank")
ggsave(p8, filename ="Results/04_fig4_wt_neo_el.png", width = 10, height = 10)


# 4) Expression level vs rank ----------------------------------
p10 <-  scatterplot_function(data = data_single_peptides,
                     x = 'expression_level', 
                     y= 'mut_mhcrank_el')+
  labs(y = "Neoepitope eluted ligand %Rank", 
       x =  "Expression level") +
  theme(plot.title = element_text(hjust = 0.5)) 
ggsave(p10, filename ="Results/04_expression_rank.png", width = 12, height = 7)

# GGseq logo plot ---------------------------------------------------------
#plot seq log, there only exist responses in length 9 to 10 so that is the only one illustarted
plot_list <- list()
i=1
for (num in 9:10) {
  for(r in c("yes","no"))  {
    p <- seqloggo_generator(data  = data_single_peptides, 
                            len = num,
                            resp = r)
    plot_List[[i]] <- p
    i <- i+1
  }
}

# Combine plots
png(file = "Results/04_GGseq_plot.png", width = 1200, height = 600)
ggdraw() +
   draw_plot(plot_list[[1]], 0, .51, .45, .4) +
   draw_plot(plot_list[[2]], 0, .0, .45, .4) +
   draw_plot(plot_list[[3]], .47, .51, .45, .4) +
   draw_plot(plot_list[[4]], .47, .0, .45, .4) +
   draw_plot_label(c("Response", "No response"), c(-0.03,  -0.03), c(1,  .48), size = 22)
dev.off()
