library(shiny)
library(SAMURAI)
data(greentea)
data(Hpylori)
data(Fleiss1993)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Sensitivity Analysis of a Meta-Analytic Effect Size"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Choose dataset ----
      selectInput("dataset", "Choose a dataset:",
                  choices = c("greentea", "Hpylori")
      ),      
            
      # Copy the line below to make a select box 
      selectInput(
        "studies_outlook1", 
        label = h4("Select outlook of unpublished studies in Plot 1"), 
        choices = list("reset",
                       "-- --",
                       "very positive", "positive", "no effect", "negative", "very negative", 
                       "-- based on published studies --",
                       "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL"), 
        selected = 1
      ),
  
      selectInput(
        "studies_outlook2", 
        label = h4("Select outlook of unpublished studies in Plot 2"), 
        choices = list("reset",
                       "-- --",
                       "very positive", "positive", "no effect", "negative", "very negative", 
                       "-- based on published studies --",
                       "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL"), 
        selected = 1
      )
    
    # actionButton(inputId = "reset", label = "Reset")  
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot 1", plotOutput(outputId = "plot1")),
                  tabPanel("Plot 2", plotOutput(outputId = "plot2"))
      )      
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Reactive value for selected dataset ----
  inputData <- reactive({
    switch(input$dataset,
           "greentea" = greentea,
           "Hpylori" = Hpylori)
  })  
  
  qbinary <- reactive({
    switch(input$dataset,
           "greentea" = FALSE,
           "Hpylori" = TRUE)
  })  

  qmean.sd <- reactive({
    switch(input$dataset,
           "greentea" = TRUE,
           "Hpylori" = FALSE)
  })
  
  qhigher.is.better <- reactive({
    switch(input$dataset,
           "greentea" = FALSE,
           "Hpylori" = FALSE)
  })
  
  output$plot1 <- renderPlot({
    
    if (input$studies_outlook1 != "reset" && substring(input$studies_outlook1, 1, 1) != '-' ){
      uoutlook <- input$studies_outlook1
    }
    else {
      uoutlook <- NA
    }

    forestsens(inputData(), 
               binary = qbinary(), 
               mean.sd = qmean.sd(), 
               higher.is.better = qhigher.is.better(),
               outlook = uoutlook)

  })

  output$plot2 <- renderPlot({
    
    if (input$studies_outlook2 != "reset" && substring(input$studies_outlook2, 1, 1) != '-' ){
      uoutlook <- input$studies_outlook2
    }
    else{
      uoutlook <- NA
    }
    
    forestsens(inputData(), 
               binary = qbinary(), 
               mean.sd = qmean.sd(), 
               higher.is.better = qhigher.is.better(),
               outlook = uoutlook)
    
  })
  
    
  # observeEvent(
  #   input$reset,
  #   { output$plot <- renderPlot({forestsens(greentea,
  #                           binary=FALSE,
  #                           mean.sd=TRUE,
  #                           higher.is.better=FALSE,
  #                           outlook = NA)})
  #   }
  # )
  
}

shinyApp(ui = ui, server = server)