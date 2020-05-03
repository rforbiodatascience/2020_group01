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




# 

