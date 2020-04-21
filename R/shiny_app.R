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

Mupexi_file <- read_xlsx("/Users/annbor/Documents/courses/R_course_22100/2020_group01_project/data/raw/ct26_library_mupexi.xlsx")


Mupexi_file[c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")] <- sapply(Mupexi_file[c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")],as.numeric)
# make my own function 
#vchoices <- Mupexi_file[c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")]
#names(vchoices) <- names(Mupexi_file[c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")])
# Define UI for application that draws a histogram
Mupexi_file$Mutation_Consequence
#data <- as.data.frame(Mupexi_file[,c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity","HLA_allele","Mutation_Consequence")])
data <- Mupexi_file
library(shiny)


# Define UI for application 
ui <- fluidPage(
  titlePanel("Scatterplot"),
  checkboxInput("logarithmicX", "show x-axis in log10", FALSE),
  checkboxInput("logarithmicY", "show y-axis in log10", FALSE),
  
  br(),
  
  # Sidebar layout with a input and output definitions
  sidebarLayout(
    # Inputs
    sidebarPanel(
      # Select variable for y-axis
      selectInput(inputId = "y", label = "Y-axis:",
                  choices = colnames(data[,c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")])
      ),
      # Select variable for x-axis
      selectInput(inputId = "x", label = "X-axis:",
                  choices = colnames(data[,c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")])
      ),
      selectInput(inputId = "ColorVar", label = "Color stuff",
                  choices = colnames(data[c("HLA_allele","Mutation_Consequence")])
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
          plotOutput(outputId = "scatterplot")
        )
        
      ),
      
      width = 6  
    )
  )
)

# Define server function required to create the scatterplot
server <- function(input, output) {
  
  # Create scatterplot object the plotOutput function is expecting
  output$scatterplot <- renderPlot({
    
    
    p <-    ggplot(data = data, aes_string(x = input$x, y = input$y)) +
      geom_point(aes_string(color = input$ColorVar))
    if(input$logarithmicX)
      p <- p + scale_x_log10()
    
    if(input$logarithmicY)
      p <- p + scale_y_log10()
    
    return(p)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
