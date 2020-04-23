# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")

# Define functions
# ------------------------------------------------------------------------------
#source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv")

# Wrangle data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Model data
# ------------------------------------------------------------------------------
#my_data_clean_aug %>% ...

# Visualise data
# ------------------------------------------------------------------------------
treatment_names <- list(
  'no'="Treatment: No",
  'yes'="Treatment: Yes"
)

treatment_labeller <- function(variable,value){
  return(treatment_names[value])
}

##Geom_boxplot self-similarity on y and x = response yes and no, facet by treatment 
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = response, y = Self_Similarity, fill = response)) +
  geom_boxplot() +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Self-Similarity vs Response", x = "Response", y = "Self-Similarity")

##Geom_boxplot Mut_rank_el on y and x = response yes and no, facet by treatment
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = response, y = Mut_MHCrank_EL)) +
  geom_boxplot(fill="steelblue") +
  facet_wrap(treatment~., labeller=treatment_labeller) +
  labs(title ="Box plot of Mutation Elution (Rank score) vs Response", x = "Response", y = "Mutated Elution (Rank score)")

##Bar plot illustrating the frequency of mutant position with use of column peptide_position
lib <- matrix(nrow = 1, ncol = 11)
colnames(lib) <- c("pos1","pos2","pos3","pos4","pos5","pos6","pos7","pos8","pos9","pos10","pos11")
lib[1, ] <- c(0,0,0,0,0,0,0,0,0,0,0)
lib <- as.data.frame(lib)

for (pos in 1:length(lib)) {
  for (position in my_data_clean_aug$peptide_position) {
    if (grepl(":", position)==T) {
      pos1 <-  as.numeric(str_split(position, ':')[[1]][1]  )
      pos2 <-  as.numeric( str_split(position, ':')[[1]][2]  )
      posistions <- pos1:pos2
      if (pos %in% posistions) {
        lib[,pos] = lib[,pos]+1}
    }
    else if (grepl(":", position)==F) {
      print(position)
      position <- as.numeric(position)
      if (pos %in% position) {
        lib[,pos] = lib[,pos]+1 }
    }
  }
  
}

lib <- as.data.frame(t(lib))
lib$pos <- rownames(lib)

ggplot(lib, aes(x=factor(pos,
                         levels = c("pos1" , "pos2" , "pos3" , "pos4" , "pos5" , "pos6" ,"pos7" , "pos8" , "pos9" , "pos10", "pos11")),
                y=V1)) +
  geom_bar(stat = "identity", fill="steelblue") +
  labs(title ="Bar plot of Peptide Position", x = "Peptide Position", y = "Number Changes")

##Bar plot illustrating the number of F(frameshift), M(missense), I(inframe insertion),  
## D(inframe deletion). Mutations are under "Mutation_Consequence" column.
my_data_clean_aug %>% 
  ggplot(mapping = aes(x = Mutation_Consequence)) +
  geom_bar(fill="steelblue") +
  labs(title ="Bar plot of Mutation Consequence", x = "Mutation Consequence", y = "Count")

# Write data
# ------------------------------------------------------------------------------
#write_tsv(...)
#ggsave(...)