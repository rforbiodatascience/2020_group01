# packges functions


#' Multiply function
#'
#' This function can make par plots
#' @param num
#' @keywords number
#' @export
#' @examples
#' cat_function()
bar_plot_func <- function(num) {
  my_data_clean_aug %>%
    filter(str_length(mut_peptide)==num,mutation_consequence=="M") %>%
    ggplot(aes(x=peptide_position)) +
    geom_bar(aes(fill = response), stat = "count")+
    scale_y_log10() +
    scale_x_discrete(limits = c(1:num)) +
    theme_bw()

}
bar_plot_func(9)

setwd("/Users/annbor/Documents/courses/R_course_22100/2020_group01_project/packages/NEOPACK")
document()


