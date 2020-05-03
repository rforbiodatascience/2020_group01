# Define project functions
# ------------------------------------------------------------------------------


### seq logo generater
# seqloggo_generator <-  function(data = my_data_clean_aug ,
#                                 len = 9,
#                                 resp = c("yes","no"),
#                                 cons = c("F","M","I","D"),
#                                 mouse = c("CT26","4T1"))
# { p <- my_data_clean_aug %>%
#     filter(str_length(mut_peptide)==len,response==resp,mutation_consequence==cons,cell_line == mouse) %>%
#     select(mut_peptide) %>%
#     ggseqlogo()
# return(p)
# }

## Bar plot function 
bar_plot_func <- function(num) {
  my_data_clean_aug %>%
    filter(str_length(mut_peptide)==num,mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() + 
    scale_x_discrete(limits = c(1:num)) +
    theme_bw()
  
}



scatterplot_function <- function(x,y) {  my_data_clean_aug %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_point(aes(color=response, alpha = response))+
    scale_y_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    scale_x_log10(breaks = c(0.01, 0.10, 1.00, 2.00, 10))+
    theme_bw()
}


# create box_function for boxplot
box_function <- function(x,y) {  
  my_data_clean_aug %>% 
    ggplot(mapping = aes_string(x = x, y = y)) +
    geom_boxplot(aes(fill = response)) }



