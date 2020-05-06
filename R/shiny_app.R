# This is a shiny app for exploring data in mupexi files with screened neoepitope responses 
# The app is en Beta version 
# The app is only for exploring data 

library(openxlsx)
library(shiny)
library(readxl)
library(ggplot2)
library(tidyverse)

my_data_clean_aug <- read_tsv(file = "data/03_my_data_clean_aug.tsv") 

Exploring_data_shiny <- function(Plotting_data) {
# Define UI for application 
ui <- fluidPage(
  titlePanel("Ploting away :) "),
  # checkboxes for log scale 
  checkboxInput("logarithmicX", "show x-axis in log10", FALSE),
  checkboxInput("logarithmicY", "show y-axis in log10", FALSE),
  
  # check boxes for plotting with ONLY responses 
  checkboxInput("selected_dat", "Plot only responses", FALSE),
  br(),
  # Sidebar layout with a input and output definitions
  sidebarLayout(
    # Inputs
    sidebarPanel(
      # Select choices to put in the y-axis
      selectInput(inputId = "y", label = "Y-axis:",
                  choices = Plotting_data %>% 
                    select(expression_level,mut_mhcscore_el,norm_mhcrank_el ,self_similarity) %>% 
                    colnames() 
                  ),
               
      # Select choices to put in the x-axis
      selectInput(inputId = "x", label = "X-axis:",
                  choices = Plotting_data %>% 
                    select(expression_level,mut_mhcscore_el,norm_mhcrank_el ,self_similarity,response,hla,mutation_consequence) %>% 
                    colnames() 
                  ),
      # select choices to color stuff
      selectInput(inputId = "ColorVar", label = "Color stuff",
                  Plotting_data %>% 
                    select(hla,mutation_consequence,response,organ) %>% 
                    colnames() 
                  ),
     # select choises for the alpha scale 
      selectInput(inputId = "alpha", label = "alpha",
                  choices = Plotting_data %>% 
                    select(response,expression_level,mut_mhcscore_el,norm_mhcrank_el ,self_similarity) %>% 
                    colnames() 
                  ),
      # select choises for facet 
      selectInput(inputId = "facet", label = "facet",
                  choices = c("none",colnames(Plotting_data[,c("response")]))
                  , selected = "response"),
     # select choises for plot type 
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
    
    # if statement for plotting only responses
    if(input$selected_dat)
      Plotting_data <- Plotting_data %>% filter(response=="yes")
    
    # switch for plot type 
    renderText({
      switch(input$plot.type,
             "boxplot" 	= 	"Boxplot",
             "dotplot" =	"dotplot")
    })
    
    # if statement for duing facet or not 
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
    # if statement for making log scales   
    if(input$logarithmicX)
      plot.type <- plot.type + scale_x_log10()
    
    if(input$logarithmicY)
      plot.type <- plot.type + scale_y_log10()
    
    return(plot.type)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
}

###########################################







