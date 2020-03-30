library(tidyverse)
library(readxl)
library(writexl)

file.choose()
## Load multiple sheets into one single dataframe
path <- "data/raw/barracoda_output_CT26.xlsx"

all_ct26_barracoda <- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path, sheet = .x), .id = "sheet")


## Add response column based on cut-off: Log2>2 & p-value<0.05
ct26_barracoda <- all_ct26_barracoda %>% mutate(response = case_when(
    log_fold_change >= 2 & p <= 0.05 & count.1 >= input.1 ~ "yes",
    TRUE ~ "no")) %>% 
  mutate(identifier = paste(HLA, Peptide.name, sep="_")) 

responses <- all_ct26_barracoda %>% filter(response=="yes") %>% distinct(identifier, .keep_all = T)

# Sample de-selected because when choosing uniques identifiers, all come from C1. 
ct26_barracoda <- ct26_barracoda %>% 
  group_by(sample) %>% 
  distinct(identifier, .keep_all = T)

write_xlsx(ct26_responses, "data/ct26_responses_barracoda.xlsx")

ct26_responses %>% group_by(sample, response) %>% count() 
