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
                  choices = c(colnames(Plotting_data[,c("expression_level","mut_mhcscore_el","norm_mhcrank_el" ,"self_similarity")]))
      ),
      # Select variable for x-axis
      selectInput(inputId = "x", label = "X-axis:",
                  choices = c(colnames(Plotting_data[,c("expression_level","mut_mhcscore_el","norm_mhcrank_el" ,"self_similarity","response","hla","mutation_consequence")]))
      ),
      selectInput(inputId = "ColorVar", label = "Color stuff",
                  choices = c(colnames(Plotting_data[,c("hla","mutation_consequence","response","organ")])),
                  
      ),
      selectInput(inputId = "alpha", label = "alpha",
                  choices = c(colnames(Plotting_data[,c("response")]))
                  ),
      selectInput(inputId = "facet", label = "facet",
                  choices = c("none",colnames(Plotting_data[,c("response")]))
                  , selected = "response"),
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
    if(input$facet!="none")
    plot.type<-switch(input$plot.type,
                      "boxplot" 	= Plotting_data %>% ggplot(aes_string(x = input$x, y = input$y)) +
                        geom_boxplot(aes_string(color = input$ColorVar), alpha=0.5) + 
                        theme_bw() + 
                        facet_grid(~get(input$facet)),
                    
                      "dotplot" = Plotting_data %>%	ggplot(aes_string(x = input$x, y = input$y)) +
                        geom_point(aes_string(color = input$ColorVar, alpha= input$alpha)) +
                        theme_bw() + 
      facet_grid(~get(input$facet)))
    
    if(input$facet=="none")
      plot.type<-switch(input$plot.type,
                        "boxplot" 	= Plotting_data %>% ggplot(aes_string(x = input$x, y = input$y)) +
                          geom_boxplot(aes_string(color = input$ColorVar), alpha=0.5) + 
                          theme_bw() ,
                        
                        "dotplot" = Plotting_data %>%	ggplot(aes_string(x = input$x, y = input$y)) +
                          geom_point(aes_string(color = input$ColorVar, alpha= input$alpha)) +
                          theme_bw() )
      
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



my_data_clean_aug %>% group_split(cell_line) %>% select(sample)





