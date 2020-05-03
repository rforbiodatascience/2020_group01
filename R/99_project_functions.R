# Define project functions
# ------------------------------------------------------------------------------


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

## Bar plot function 
bar_plot_func <- function(num) {
  my_data_clean_aug %>%
    filter(str_length(mut_peptide)==num,mutation_consequence=="M") %>% 
    ggplot(aes(x=peptide_position)) + 
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() + 
    theme_bw()
  
}


# 

