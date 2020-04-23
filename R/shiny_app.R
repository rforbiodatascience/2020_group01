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
  titlePanel("plots"),
  checkboxInput("logarithmicX", "show x-axis in log10", FALSE),
  checkboxInput("logarithmicY", "show y-axis in log10", FALSE),
  
  br(),
  #[,c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")]
  # Sidebar layout with a input and output definitions
  sidebarLayout(
    # Inputs
    sidebarPanel(
      # Select variable for y-axis
      selectInput(inputId = "y", label = "Y-axis:",
                  choices = colnames(Plotting_data[,c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity")])
      ),
      # Select variable for x-axis
      selectInput(inputId = "x", label = "X-axis:",
                  choices = colnames(Plotting_data[,c("Expression_Level","Mut_MHCrank_EL", "Self_Similarity","response")])
      ),
      selectInput(inputId = "ColorVar", label = "Color stuff",
                  choices = colnames(Plotting_data[,c("HLA_allele","Mutation_Consequence","response","HLA","organ")])
      ),
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
    
    renderText({
      switch(input$plot.type,
             "boxplot" 	= 	"Boxplot",
             "dotplot" =	"dotplot")
    })
    plot.type<-switch(input$plot.type,
                      "boxplot" 	= ggplot(data = Plotting_data, aes_string(x = input$x, y = input$y)) +
                        geom_boxplot(aes_string(color = input$ColorVar), alpha=0.5) +
                        facet_grid(.~input$facet),
                    
                      "dotplot" =	ggplot(data = Plotting_data, aes_string(x = input$x, y = input$y)) +
                        geom_point(aes_string(color = input$ColorVar))) +
                        facet_grid(.~input$facet)
    
    
   # p <-    ggplot(data = Plotting_data, aes_string(x = input$x, y = input$y)) +
  #    geom_point(aes_string(color = input$ColorVar))
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


library(shiny)
library(ggplot2)
library(purrr)
library(dplyr)


#example data
# data(iris)
# 
# #make some factors
# #easier to let ggplot2 control plotting (color, fill) based on type
# data(mtcars)
# uvals<-sapply(mtcars,function(x){length(unique(x))})
# mtcars<-map_if(mtcars,uvals<4,as.factor) %>%
#   as.data.frame()


#plotting theme for ggplot2
.theme<- theme(
  axis.line = element_line(colour = 'gray', size = .75),
  panel.background = element_blank(),
  plot.background = element_blank()
)


# UI for app
ui<-(pageWithSidebar(
  # title
  headerPanel("Select Options"),
  
  #input
  sidebarPanel
  (
    # Input: Select a file ----
    
  #  fileInput("file1", "Choose CSV File",
    #          multiple = TRUE,
     #         accept = c("text/csv",
      #                   "text/comma-separated-values,text/plain",
       #                  ".csv")),
    # Input: Checkbox if file has header ----
    #checkboxInput("header", "Header", TRUE),
    
    # Input: Select separator ----
  #  radioButtons("sep", "Separator",
  #               choices = c(Semicolon = ";",
  #                           Comma = ",",
  #                           Tab = "\t"),
  #               selected = ","),
    # Horizontal line ----
  #  tags$hr(),
    
    
    # Input: Select what to display
  #  selectInput("dataset","Data:",
  #              choices =list(iris = "iris", mtcars = "mtcars",
   #                           uploaded_file = "inFile"), selected=NULL),
    selectInput("variable","Variable:", choices = NULL),
    selectInput("group","Group:", choices = NULL),
    selectInput("plot.type","Plot Type:",
                list(boxplot = "boxplot", histogram = "histogram", density = "density", bar = "bar")
    ),
    checkboxInput("show.points", "show points", TRUE)
  ),
  
  # output
  mainPanel(
    h3(textOutput("caption")),
    #h3(htmlOutput("caption")),
    uiOutput("plot") # depends on input
  )
))


# shiny server side code for each call
server<-(function(input, output, session) {
  
  #update group and
  #variables based on the data
#  observe({
    #browser()
#    if(!exists(input$dataset)) return() #make sure upload exists
#    var.opts<-colnames(get(input$dataset))
#    updateSelectInput(session, "variable", choices = var.opts)
#    updateSelectInput(session, "group", choices = var.opts)
#  })
  
  output$caption<-renderText({
    switch(input$plot.type,
           "boxplot" 	= 	"Boxplot",
           "histogram" =	"Histogram",
           "density" 	=	"Density plot",
           "bar" 		=	"Bar graph")
  })
  
  
  output$plot <- renderUI({
    plotOutput("p")
  })
  
  #get data object
  get_data<-reactive({
    
 #   if(!exists(input$dataset)) return() # if no upload
    
#    check<-function(x){is.null(x) || x==""}
#    if(check(input$dataset)) return()
    
    obj<-list(data= Plotting_data,
              variable=Plotting_data$variable,
              group=Plotting_data$group
    )
    
    #require all to be set to proceed
#    if(any(sapply(obj,check))) return()
    #make sure choices had a chance to update
 #   check<-function(obj){
 #     !all(c(obj$variable,obj$group) %in% colnames(obj$data))
#    }
    
 #   if(check(obj)) return()
    
    
 #   obj
    
  })
  
  #plotting function using ggplot2
  output$p <- renderPlot({
    
    plot.obj<-get_data()
    
    #conditions for plotting
    if(is.null(plot.obj)) return()
    
    #make sure variable and group have loaded
    if(plot.obj$variable == "" | plot.obj$group =="") return()
    
    #plot types
    plot.type<-switch(input$plot.type,
                      "boxplot" 	= geom_boxplot(),
                      "histogram" =	geom_histogram(alpha=0.5,position="identity"),
                      "density" 	=	geom_density(alpha=.75),
                      "bar" 		=	geom_bar(position="dodge")
    )
    
    
    if(input$plot.type=="boxplot")	{		#control for 1D or 2D graphs
      p<-ggplot(plot.obj$data,
                aes_string(
                  x 		= plot.obj$group,
                  y 		= plot.obj$variable,
                  fill 	= plot.obj$group # let type determine plotting
                )
      ) + plot.type
      
      if(input$show.points==TRUE)
      {
        p<-p+ geom_point(color='black',alpha=0.5, position = 'jitter')
      }
      
    } else {
      
      p<-ggplot(plot.obj$data,
                aes_string(
                  x 		= plot.obj$variable,
                  fill 	= plot.obj$group,
                  group 	= plot.obj$group
                  #color 	= as.factor(plot.obj$group)
                )
      ) + plot.type
    }
    
    p<-p+labs(
      fill 	= input$group,
      x 		= "",
      y 		= input$variable
    )  +
      .theme
    print(p)
  })
  
  # set uploaded file
  # upload_data<-reactive({
  #   
  #   inFile <- input$file1
  #   
  #   if (is.null(inFile))
  #     return(NULL)
  #   
  #   #could also store in a reactiveValues
  #   read.csv(inFile$datapath,
  #            header = input$header,
  #            sep = input$sep)
  # })
  
  # observeEvent(input$file1,{
  #   inFile<<-upload_data()
  # })
  
  
})


# Create Shiny app ----
shinyApp(ui, server)


