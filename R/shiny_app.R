#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(openxlsx)
library(shiny)
library(readxl)
library(ggplot2)
library(tidyverse)

Plotting_data <- read_tsv(file = "data/03_my_data_clean_aug.tsv")

# Define UI for application 
ui <- fluidPage(
  titlePanel("Ploting away :) "),
  checkboxInput("logarithmicX", "show x-axis in log10", FALSE),
  checkboxInput("logarithmicY", "show y-axis in log10", FALSE),
  checkboxInput("selected_dat", "Plot only responses", FALSE),
  br(),
  # Sidebar layout with a input and output definitions
  sidebarLayout(
    # Inputs
    sidebarPanel(
      # Select variable for y-axis
      selectInput(inputId = "y", label = "Y-axis:",
                  choices = colnames(Plotting_data[,c("Expression_Level","Mut_MHCrank_EL","Norm_MHCrank_EL" ,"Self_Similarity")])
      ),
      # Select variable for x-axis
      selectInput(inputId = "x", label = "X-axis:",
                  choices = c(colnames(Plotting_data[,c("Expression_Level","Mut_MHCrank_EL", "Norm_MHCrank_EL","Self_Similarity","response","HLA_allele","Mutation_Consequence","cell.line")]))
      ),
      selectInput(inputId = "ColorVar", label = "Color stuff",
                  choices = c(colnames(Plotting_data[,c("HLA_allele","Mutation_Consequence","response","organ","cell.line")])),
                  
      ),
   #   selectInput(inputId = "facet", label = "facet",
  #                choices = colnames(Plotting_data[,c("Mutation_Consequence","response","HLA","organ")])
  #    ),
      selectInput("plot.type","Plot Type:",
                  list(boxplot = "boxplot", dotplot = "dotplot")
      ),
      
      
      
      width = 6
    ),
    
    # Output:
    mainPanel(
      # Create a container for tab panels
      tabsetPanel(
        
        tabPanel(
          title = "Explore the Data",
          # Show scatterplot
          plotOutput(outputId = "plots")
        )
        
      ),
      
      width = 6  
    )
  )
)

# Define server function required to create the scatterplot
server <- function(input, output) {
  
  # Create scatterplot object the plotOutput function is expecting
  output$plots <- renderPlot({
    
    if(input$selected_dat)
      Plotting_data <- Plotting_data %>% filter(response=="yes")
    
    renderText({
      switch(input$plot.type,
             "boxplot" 	= 	"Boxplot",
             "dotplot" =	"dotplot")
    })
    plot.type<-switch(input$plot.type,
                      "boxplot" 	= Plotting_data %>% ggplot(aes_string(x = input$x, y = input$y)) +
                        geom_boxplot(aes_string(color = input$ColorVar), alpha=0.5) + 
                        theme_bw(),
                    
                      "dotplot" = Plotting_data %>%	ggplot(aes_string(x = input$x, y = input$y)) +
                        geom_point(aes_string(color = input$ColorVar))) +
                        theme_bw()
    
    if(input$logarithmicX)
      plot.type <- plot.type + scale_x_log10()
    
    if(input$logarithmicY)
      plot.type <- plot.type + scale_y_log10()
    
    return(plot.type)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


###########################################
