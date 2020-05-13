# Define project functions

# ------------------------------------------------------------------------------
# define some colors 
respond_cols <- c("#91bfdb","#ef8a62")

# Baracoda respond function -----------------------------------------------
barc_resp <- function(data = my_data_clean_aug, 
                      mouce_cell_line = "CT26" ){
  
  p <- data %>% filter(cell_line == mouce_cell_line) %>% 
    ggplot(., aes(peptide_name, log_fold_change)) +
    geom_point(aes(color = sample, shape = organ, 
                   alpha = response, size = estimated_frequency_norm)) +
    geom_text_repel(data  %>%
                      filter(cell_line == mouce_cell_line, response == "yes"),
                    mapping = aes(label = peptide_name, size = 14)) +
    facet_grid(treatment ~.) +
    labs(size = "Normalized estimated frequency %",
         shape = "Organ", 
         color = "Sample", 
         alpha = "Response", 
         y = "logFC") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16), 
          axis.title = element_text(size = 16),
          axis.title.x = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           alpha = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))
  return(p)
}

# GGseq logo function -----------------------------------------------------
seqloggo_generator <-  function(data = my_data_clean_aug ,
                                len = 9,
                                resp = c("yes","no")){
 p <- my_data_clean_aug %>%
  filter(str_length(mut_peptide)==len,response==resp) %>% 
  select(mut_peptide) %>%
  ggseqlogo()
return(p)
}


# Bar plot function -------------------------------------------------------

bar_plot_func <- function(data = my_data_clean_aug,
                            pep_length = 9 ,
                          no_legend = TRUE) {
 p <- data %>%
    filter(str_length(mut_peptide)==pep_length,mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
   # scale_y_log10() + 
    scale_x_discrete(limits = factor(1:pep_length)) +
    scale_fill_manual(values = respond_cols) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))+
    facet_grid(cell_line~.)+
    labs(x = "Peptide Position", 
         y = "Amount neopeptides",
         fill = "Response")
 
  # Determine legend 
  if (no_legend == TRUE) p <- p + theme(legend.position = 'none') else NULL 

 return(p)
}


# Scatterplot function ----------------------------------------------------

scatterplot_function <- function(data = data_single_peptides,
                                 x = 'mut_mhcrank_el',
                                 y = 'expression_level',
                                 no_legend = TRUE) {
   p <- data %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_point(aes(color=response, alpha = response, size  = estimated_frequency_norm))+
    scale_y_log10(breaks = c(0.01, 0.10, 0.5, 1.00, 2.00, 10))+
    scale_x_log10(breaks = c(0.01, 0.10, 0.5, 1.00, 2.00, 10))+
    theme_bw() + 
    scale_alpha_manual(breaks = c("no","yes"),labels = c("no","yes"),values = c(0.3,0.9))+
    scale_color_manual(values = respond_cols) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    facet_grid(cell_line~., scales = "free") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))+
    labs(size = "Normalized estimated frequency %",
         color = "Response", 
         alpha = "Response")
  # Determine legend 
  if (no_legend == TRUE) p <- p + theme(legend.position = 'none') else NULL 
  
  return(p)
  
}


# Box plot function -------------------------------------------------------

box_function <- function(data = data_single_peptides,
                         x ='response',
                         y= 'mut_mhcrank_el',
                         no_legend = TRUE) {  
   p <-  data %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_quasirandom(aes(color = response),size = 2) + 
    geom_boxplot(aes(fill = response), 
                 alpha = .5, outlier.shape = NA, colour = '#525252') +
    facet_grid(cell_line~., scales = "free") +
    theme_bw() +
    theme(plot.title = element_text(size = 14,hjust = 0.5), 
           axis.title = element_text(size = 14),
           legend.text = element_text(size = 12),
           legend.title = element_text(size = 14))+
    scale_fill_manual(values = respond_cols) +
    scale_color_manual(values = respond_cols) +
    guides(fill = FALSE, color = guide_legend(override.aes = list(size = 4)))
   
   # Determine legend 
   if (no_legend == TRUE) p <- p + theme(legend.position = 'none') else NULL 
   
   return(p)
  }



