# Define project functions

# ------------------------------------------------------------------------------
# define some colors 
respond_cols <- c("#91bfdb","#ef8a62")

### seq logo generater
seqloggo_generator <-  function(data = my_data_clean_aug ,
                                len = 9,
                                resp = c("yes","no"),
                                cons = c("F","M","I","D"),
                                mouse = c("CT26","4T1"))
{ p <- my_data_clean_aug %>%
    filter(str_length(mut_peptide)==len,response==resp,mutation_consequence==cons,cell_line == mouse) %>%
    select(mut_peptide) %>%
    ggseqlogo()
return(p)
}
#


## Bar plot function 
bar_plot_func <- function(num) {
  my_data_clean_aug %>%
    filter(str_length(mut_peptide)==num,mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() + 
    scale_x_discrete(limits = c(1:num)) +
    scale_fill_manual(values = respond_cols) +
    theme_bw()
  
}


scatterplot_function <- function(data = my_data_clean_aug,
                                 x = 'mut_mhcrank_el',
                                 y = 'expression_level') 
  {  data %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_point(aes(color=response, alpha = response, size  = estimated_frequency))+
    scale_y_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    scale_x_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    theme_bw() + 
    scale_alpha_manual(breaks = c("no","yes"),labels = c("no","yes"),values = c(0.3,0.9))+
    scale_color_manual(values = respond_cols) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    labs(size = "Estimated frequency",
         color = "Respond", 
         alpha = "Respond")
}


# create box_function for boxplot
box_function <- function(x,y) {  
  my_data_clean_aug %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_boxplot(aes(fill = response)) }



