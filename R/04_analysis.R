
# Clear workspace ---------------------------------------------------------
rm(list = ls())

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggplot2)
#install.packages(ggseqlogo)
library(ggseqlogo)
library(ggseqlogo)
library(cowplot)
library(ggbeeswarm)


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")


# Wrangle data ------------------------------------------------------------
# arrange data so "no" is plotted first 
my_data_clean_aug <- my_data_clean_aug %>%  arrange(response)
# select unique peptides 
data_single_peptides <- my_data_clean_aug %>% 
  group_by(response) %>% 
  distinct(identifier, .keep_all = T) %>% 
  ungroup()


# Visualise data ----------------------------------------------------------
# 1) Barracoda characteristics --------------------------------------------
p1_CT26 <- barc_resp(my_data_clean_aug, "CT26") +
  labs(title = "log fold change of CT26-tumor cell line neoepitope screen")

p1_4T1 <- barc_resp(my_data_clean_aug, "4T1") +
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
  facet_grid(vars(cell_line)) +
  labs(x = "p-value",
       y = "log fold change",
       size = "Nomalized estimated frequency")
ggsave(p2, filename ="Results/04_p2_barracoda_selection.png", width = 12, height = 7)

# 2) Boxplots response vs non-response  ---------------------------------------------
# Elution rank ------------------------------------------------------------
p3 <-  box_function(data = data_single_peptides, 
                    x = 'response',
                    y = 'mut_mhcrank_el') +
  labs(title ="Elution (EL) rank score",
       x = "Response", 
       y = "Neoepitopes elution (%rank score)")

# Binding affinity --------------------------------------------------------
p4 <-  box_function(data = data_single_peptides, 
                    x = 'response',
                    y= 'mut_mhcrank_ba') +
  labs(title ="Binding affinity (BA) rank score",
       x = "Response", 
       y = "Neoepitopes binding affinity (%rank score)")

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
                          y = 'norm_mhcrank_el')+
  labs(title= "Elution (EL) of neoepitope vs WT epitope", 
       x= "Neoepitope EL %Rank ", 
       y="WT epitope EL %Rank")

p9 <- scatterplot_function(data = data_single_peptides,
                           x = 'mut_mhcrank_ba', 
                           y= 'norm_mhcrank_ba',
                           no_legend = FALSE)+
  labs(title= "Binding affinity (BA) of neoepitope vs WT epitope",
       x= "Neoepitope BA %Rank ",
       y="WT epitope BA %Rank")

bind_legend <- get_legend(p9)
p9 <- p9 + theme(legend.position = "none")

png(file = "Results/04_fig4_wt_neo_ba_el.png", width = 800, height = 400)
ggdraw() +
  draw_plot(p8, 0, .0, .43, .95) +
  draw_plot(p9, .44, .0, .43, .95) +
  draw_plot(bind_legend, .88, .4, .2, .2) +
  draw_plot_label(c("A", "B"), c(0,  0.45), c(1, 1), size = 20)
dev.off()

# 4) Expression level vs rank ----------------------------------
p10 <-  scatterplot_function(data = data_single_peptides,
                     x = 'expression_level', 
                     y= 'mut_mhcrank_el')+
  labs(y = "Neoepitoe elution (% rank score)", 
       x =  "Expression level") +
  theme(plot.title = element_text(hjust = 0.5)) 
ggsave(p10, filename ="Results/04_expression_rank.png", width = 12, height = 7)



# Mutaions overview -------------------------------------------------------
p11 <- my_data_clean_aug %>% 
  ggplot(aes(x = mutation_consequence)) +
  geom_bar(aes(fill = response), stat = "count") +
  scale_y_log10() +
  theme_bw() +
  scale_fill_manual(values = respond_cols) +
  facet_wrap(cell_line~.) +
  labs(title ="Bar plot of Mutation Consequence (by Cell Line)", 
       x = "Mutation Consequence", 
       y = "Aomunt of neoeptiopes",
       fill = "Response")
ggsave(p11, filename ="Results/04_Mutation_consequence.png", width = 10, height = 8)


# GGseq logo plot ---------------------------------------------------------
#plot seq log, there only exist responses in length 9 to 10 so that is the only one illustarted
plot_List <- list()
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
   draw_plot(plot_List[[1]], 0, .51, .45, .4) +
   draw_plot(plot_List[[2]], 0, .0, .45, .4) +
   draw_plot(plot_List[[3]], .47, .51, .45, .4) +
   draw_plot(plot_List[[4]], .47, .0, .45, .4) +
   draw_plot_label(c("Response", "No response"), c(-0.03,  -0.03), c(1,  .48), size = 22)
dev.off()




######## data not to use 

# No. mutations by HLA faceted by cell_line. This plot was made by Sara because of curiosity.
# It shows in what MHC allele do the response peptides bind. Most of the responsive peptides
# bind H2-Kd allele in both cell lines 
# hla <- data_single_peptides %>% 
#   filter(response == "yes") %>% 
#   ggplot(aes(x = hla)) +
#   geom_bar(aes(fill = hla), stat = "count") +
#   facet_grid(vars(cell_line)) +
#   theme_bw()
# ggsave(hla, filename ="Results/hla.png", width = 12, height = 7)

# p10 <- data_single_peptides %>% 
#   ggplot(., aes(expression_level, mut_mhcrank_el)) +
#   geom_point(aes(color = response, alpha = response, size = estimated_frequency_norm))+
#   scale_x_log10() +
#   facet_grid(vars(cell_line), scales = "free") +
#   # geom_text_repel(my_data_clean_aug %>%
#   #                   group_by(cell_line) %>%
#   #                   filter(response == "yes"), stat = "identity", mapping = aes(label = peptide_name))+
#   theme_bw() + 
#   scale_color_manual(values = respond_cols) +
#   labs(size = "Estimated frequency normalized", 
#        color = "Respond", 
#        alpha = "Respond", 
#        y = "Neoepitoe elution (% rank score)", 
#        x =  "Expression level") + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(p10, filename ="Results/p10.png", width = 12, height = 7)



