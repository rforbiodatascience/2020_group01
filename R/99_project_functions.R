# Define project functions
# ------------------------------------------------------------------------------


### seq logo generater 
seqloggo_generator <-  function(data =my_data_clean_aug ,
                                len = 9,
                                resp = "yes",
                                cons = "M",
                                mouse = c("CT26","4T1")) 
{ my_data_clean_aug %>% 
    filter(str_length(mut_peptide)==len,response==resp,mutation_consequence==cons,cell_line == mouse) %>% 
    select(mut_peptide) %>% 
    ggseqlogo()}

# 

