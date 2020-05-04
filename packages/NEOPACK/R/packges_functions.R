# packges functions
library(tidyverse)

#' Multiply function
#'
#' This function can make par plots
#' @param num
#' @keywords number
#' @export
#' @examples
#' bar_plot_func()
bar_plot_func <- function(data,num) {
  data() %>%
    filter(str_length(mut_peptide)==num,mutation_consequence=="M") %>%
    ggplot(aes(x=peptide_position)) +
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() +
    scale_x_discrete(limits = c(1:num)) +
    theme_bw()

}

bar_plot_func(9)



