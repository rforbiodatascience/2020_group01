# Define project functions

# ------------------------------------------------------------------------------
# define some colors 
respond_cols <- c("#91bfdb","#ef8a62")


# Baracoda respond function -----------------------------------------------

barc_resp <- function(d, c){
  d <- d %>% filter(cell_line == c) %>% 
    ggplot(., aes(peptide_name, log_fold_change)) +
    geom_point(aes(color = sample, shape = organ, alpha = response, size = estimated_frequency_norm)) +
    geom_text_repel(d %>%
                      filter(cell_line == c, response == "yes"), mapping = aes(label = neoepitope_sequence, size = 14)) +
    facet_grid(vars(treatment)) +
    labs(size = "Normalized estimated frequency",
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
}

# barc_freq <- function(d, c, p){
#   d <- d %>% filter(cell_line == c) %>% 
#     ggplot(., aes_string("peptide_name", p)) +
#     geom_point(aes(color = organ, shape = sample), size = 3) +
#     #geom_text_repel(d %>%
#     #                 filter(response == "yes"), stat = "identity", mapping = aes(label = peptide_name)) +
#     scale_fill_manual(values = respond_cols) +
#     facet_grid(vars(treatment)) +
#     labs(#shape = "Sample", 
#       color = "Sample", 
#       #alpha = "Response", 
#       y = "Normalized estimated frequency", 
#       x = "Neoepitopes")+
#     theme_bw() +
#     guides(color = guide_legend(override.aes = list(size = 2)))
# }



# GGseq logo function -----------------------------------------------------
seqloggo_generator <-  function(data = my_data_clean_aug ,
                                len = 9,
                                resp = c("yes","no"))
{ p <- my_data_clean_aug %>%
  filter(str_length(mut_peptide)==len,response==resp) %>% #,mutation_consequence==cons,cell_line == mouse
  select(mut_peptide) %>%
  ggseqlogo()
return(p)
}


# Bar plot function -------------------------------------------------------

bar_plot_func <- function(data = my_data_clean_aug,
                            pep_length = 9 ) {
  data %>%
    filter(str_length(mut_peptide)==pep_length,mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() + 
    scale_x_discrete(limits = factor(1:pep_length)) +
    scale_fill_manual(values = respond_cols) +
    theme_bw()
}


# Scatterplot function ----------------------------------------------------

scatterplot_function <- function(data = data_single_peptides,
                                 x = 'mut_mhcrank_el',
                                 y = 'expression_level') 
  {  data %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_point(aes(color=response, alpha = response, size  = estimated_frequency_norm))+
    scale_y_log10(breaks = c(0.01, 0.10, 0.5, 1.00, 2.00, 10))+
    scale_x_log10(breaks = c(0.01, 0.10, 0.5, 1.00, 2.00, 10))+
    theme_bw() + 
    scale_alpha_manual(breaks = c("no","yes"),labels = c("no","yes"),values = c(0.3,0.9))+
    scale_color_manual(values = respond_cols) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(size = "Estimated frequency normalized",
         color = "Respond", 
         alpha = "Respond")
}


# Box plot function -------------------------------------------------------

box_function <- function(data = data_single_peptides,
                         x ='response',
                         y= 'mut_mhcrank_el') {  
    data %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_quasirandom(aes(color = response),size = 2) + 
    geom_boxplot(aes(fill = response), alpha  = 0.8 ) +
    facet_grid(vars(cell_line), scales = "free") +
    theme_bw() +
    scale_fill_manual(values = respond_cols) +
    scale_color_manual(values = respond_cols) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  }



