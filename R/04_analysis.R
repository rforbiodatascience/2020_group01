
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

# libraries for modelling 
# packages to install 
#install.packages("devtools")
#devtools::install_github("tidyverse/broom")
#remotes::install_github("dariyasydykova/tidyroc")
library(tidyverse)
library(tidyroc)
library(broom)
library(yardstick)
library(cowplot)



# Define functions --------------------------------------------------------
wd = getwd()
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
   
p1_CT26


p1_4T1 <- barc_resp(my_data_clean_aug, "4T1") +
  labs(title = "log fold change of 4T1-tumor cell line neoepitope screen")
p1_4T1 

png(file = "Results/fig1_barracoda.png", width = 1500, height = 1000)
ggdraw() +
  draw_plot(p1_4T1, .05, .52, .95, .45) +
  draw_plot(p1_CT26, .05, .0, .95, .45) +
  #draw_plot(barracoda_legend, .95, .4, .2, .2) +
  draw_plot_label(c("A", "B"), c(.0, .0), c(.95, .45), size = 22)
dev.off()

p2 <- data_single_peptides %>% 
  ggplot(aes(p_value, log_fold_change)) +
  geom_point(aes(color=response, size  = estimated_frequency_norm)) +
  scale_color_manual(values = respond_cols) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0)) +
  scale_y_continuous(breaks = c(-10, -5, 0, 2, 5)) +
  theme_bw() +
  geom_hline(yintercept = 2, linetype = "dotted") +
  geom_vline(xintercept = 0.01, linetype = "dotted") +
  facet_grid(vars(cell_line))
ggsave(p2, filename ="Results/p2_barracoda_selection.png", width = 12, height = 7)

# 2) Boxplots response vs non-response  ---------------------------------------------
# Elution rank ------------------------------------------------------------
p3 <-  box_function('response','mut_mhcrank_el') +
  labs(title ="Eslution (EL) rank score",
       x = "Response", 
       y = "Neoepitopes elution (%rank score)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# Binding affinity --------------------------------------------------------
p4 <-  box_function('response','mut_mhcrank_ba') +
  labs(title ="Binding affinity (BA) rank score",
       x = "Response", 
       y = "Neoepitopes binding affinity (%rank score)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


# Self-similarity ---------------------------------------------------------
p5 <- box_function('response','self_similarity') +
  labs(title = "Self-similarity ", 
       x = "Response", 
       y = "Self Similarity")

mut_bar_legend <- get_legend(p5)
p5 <- p5 + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

png(file = "Results/fig2_el_ba_similarity.png", width = 1200, height = 400)
ggdraw() +
  draw_plot(p3, .0, .0, .3, .95) +
  draw_plot(p4, .3, .0, .3, .95) +
  draw_plot(p5, .6, .0, .3, .95) +
  draw_plot(mut_bar_legend, .85, .4, .2, .2) +
  draw_plot_label(c("A", "B", "C"), c(.0,  .3, .6), c(.97, .97, .97), size = 22)
dev.off()

# 3) Barplots mutation position -------------------------------------------
# 9mer --------------------------------------------------------------------
p6 <- bar_plot_func(9) +
  facet_grid(vars(cell_line))+
  labs(x = "Peptide Position", #title ="Bar plot of Peptide Position (by Cell Line)", 
       y = "Count") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


# 10mer -------------------------------------------------------------------
p7 <- bar_plot_func(10) +
  facet_grid(vars(cell_line))+
  labs( x = "Peptide Position", #title ="Bar plot of Peptide Position (by Cell Line)",
        y = "Count" ,
        fill = "Respond")

mut_bar_legend <- get_legend(P4)
p7 <- p7 + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


png(file = "Results/fig_3_mutation_position.png", width = 800, height = 400)
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
       x= "Neoepitope EL %Rank ", y="WT epitope EL %Rank") + 
  facet_grid(vars(cell_line)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p9 <- scatterplot_function(data = data_single_peptides,
                           x = 'mut_mhcrank_ba', 
                           y= 'norm_mhcrank_ba')+
  labs(title= "Binding affinity (BA) of neoepitope vs WT epitope",
       x= "Neoepitope BA %Rank ", y="WT epitope BA %Rank") +
  facet_grid(vars(cell_line)) +
  theme(plot.title = element_text(hjust = 0.5))

bind_legend <- get_legend(p9)
p9 <- p9 + theme(legend.position = "none")

png(file = "Results/fig4_wt_neo_ba_el.png", width = 800, height = 400)
ggdraw() +
  draw_plot(p8, 0, .0, .43, .95) +
  draw_plot(p9, .44, .0, .43, .95) +
  draw_plot(bind_legend, .88, .4, .2, .2) +
  draw_plot_label(c("A", "B"), c(0,  0.45), c(1, 1), size = 20)
dev.off()

# 4) Expression level vs rank ----------------------------------
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

# why not with function ??
p10 <-  scatterplot_function(data = data_single_peptides,
                     x = 'expression_level', 
                     y= 'Neoepitope EL (% rank)')+
  labs(size = "Estimated frequency normalized", 
       color = "Respond", 
       alpha = "Respond", 
       y = "Neoepitoe elution (% rank score)", 
       x =  "Expression level") +
  facet_grid(vars(cell_line), scales = "free")+
  theme(plot.title = element_text(hjust = 0.5)) 
ggsave(p10, filename ="Results/p10.png", width = 12, height = 7)

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
       y = "Count")
ggsave(p11, filename ="Results/p11.png", width = 12, height = 7)

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

# plot together
png(file = "Results/GGseq_plot.png", width = 1200, height = 600)
ggdraw() +
   draw_plot(plot_List[[1]], 0, .51, .45, .4) +
   draw_plot(plot_List[[2]], 0, .0, .45, .4) +
   draw_plot(plot_List[[3]], .47, .51, .45, .4) +
   draw_plot(plot_List[[4]], .47, .0, .45, .4) +
   draw_plot_label(c("Responses", "No responses"), c(-0.05,  -0.06), c(1,  .48), size = 22)
dev.off()


# Modelling ---------------------------------------------------------------

# Roc curves
data_single_peptides <- data_single_peptides %>% 
  mutate(new_score  = expression_level/(mut_mhcrank_el)) %>%
  mutate(response_binary = case_when(response=="yes" ~ 1,
                           response=="no" ~ 0))

#### downsample negative data 
X_neg <- data_single_peptides %>% 
  filter(response_binary==0) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,
         allele_frequency,priority_score,new_score,response_binary) %>% 
  sample_n(40)
X_pos <- data_single_peptides %>% 
  filter(response_binary==1) %>% 
  select(mut_mhcrank_el,self_similarity,expression_level,
         allele_frequency,priority_score,new_score,response_binary)
# bind pos and negative dataset 
X_selected <- bind_rows(X_pos,X_neg) %>% 
  sample_n(70)

# set the second level of truth to positive in yardstick
options(yardstick.event_first = FALSE)

# Logistic regression with all chossen parameters 
glm_out1 <- glm(
  formula = response_binary ~ mut_mhcrank_el+
    self_similarity+
    expression_level,
  family = binomial,
  data = X_selected) %>%
  augment() %>%
  mutate(model = "all") 

# logistic regression with self similarity
glm_out2 <- glm(response_binary ~ self_similarity,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "self_similarity")

# logistic regression with mutant rank el 
glm_out3 <- glm(response_binary ~ mut_mhcrank_el,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "mut_mhcrank_el") 

# logistic regression with expression score 
glm_out4 <- glm(response_binary ~ expression_level,
                family = binomial,
                data = X_selected) %>%
  augment() %>%
  mutate(model = "expression_level")


# combine all the data set 
glm_out <- bind_rows(glm_out1, glm_out2,glm_out3,glm_out4)
# make binary variable as factor 
glm_out <- glm_out %>% 
  mutate(response_binary = as.factor(response_binary))

# plot ROC curves
Roc_curves <- glm_out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(truth = response_binary, .fitted) %>% # get values to plot an ROC curve
  ggplot(
    aes(
      x = 1 - specificity,  # equal to FPR(False positive rate)
      y = sensitivity, # equal to TPR ( True Positive Rate)
      color = model
    )
  ) + # plot roc curves for each models 
  geom_line(size = 1.1) +
  labs( x = "FPR", y = "TPR") +
geom_abline(slope = 1, intercept = 0, size = 0.4) +
   scale_color_manual(labels = c("All", "Self similarity","Neoepitope EL (% rank)", "Expression Level"), 
  values = c("#d7191c", "#fdae61","#abdda4","#2b83ba"), 
  breaks = c("all", "self_similarity","mut_mhcrank_el","expression_level")) +
  coord_fixed() +
  theme_bw()
ggsave(Roc_curves, file = "Results/Roc_curves.png")


AUC_values <- glm_out %>%
  group_by(model) %>% # group to get individual AUC value for each model
  roc_auc(truth = response_binary, .fitted)











