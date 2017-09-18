library(shiny)
library(SAMURAI)

## load data sets
data(greentea)
data(Hpylori)
data(Fleiss1993)

## Define UI for app
ui <- fluidPage(
  
  ## app title
  tags$h1("Sensitivity Analysis of a Meta-Analysis to Unpublished Studies"),
  
  fluidRow(

    # sidebar of inputs 
    column(4,
      wellPanel(
        
        ## select existing dataset
        selectInput("dataset", 
                    label = "Choose an existing dataset:",
                    choices = c("green tea", "Fleiss (1993)", "H.pylori"),
                    width = '100%'),       
  
        ## Choose dataset to download
        downloadButton("downloadData", "Download as a CSV file"),
        br(),
        tags$a(href='https://noory.shinyapps.io/samurai/', 'Click here to upload your own dataset (coming soon)')
      ),
        
      wellPanel(                
        ## select box 1
        selectInput("studies_outlook1", 
                    label = "Select outlook of all unpublished studies in Plot 1", 
                    choices = list("initial settings",
                                   "-- --",
                                   "very positive", "positive", "no effect", "negative", "very negative", 
                                   "-- based on published studies' effect size and CL --",
                                   "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL"), 
                    selected = "initial settings",
                    width = '100%'),
        
        ## select box 2
        selectInput("studies_outlook2", 
                    label = "Select outlook of all unpublished studies in Plot 2", 
                    choices = list("initial settings",
                                   "-- --",
                                   "very positive", "positive", "no effect", "negative", "very negative", 
                                   "-- based on published studies' effect size and CL --",
                                   "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL"), 
                    selected = "initial settings",
                    width = '100%')
      ),
      
      wellPanel(           
        ## References
        tags$h4('References:'),
        tags$a(href='https://doi.org/10.1186/2046-4053-3-27', 'Systematic Reviews 2014 3:27'),
        br(),
        tags$a(href='https://cran.r-project.org/web/packages/SAMURAI/index.html', 'SAMURAI, an R package')
      )
      
    ),
    
    ## display outputs 
    column(8,
        # Output: Tabset w/ plot, summary, and table ----
        tabsetPanel(type = "tabs",
                    tabPanel("Click here for Plot 1", plotOutput(outputId = "plot1")),
                    tabPanel("Click here for Plot 2", plotOutput(outputId = "plot2")) )      
    )  
    
  )
)

## Define server logic for app
server <- function(input, output) {
  
  ## Reactive values for dataset selected
  
  inputData <- reactive({
    switch(input$dataset,
           "green tea" = greentea,
           "H.pylori" = Hpylori,
           "Fleiss (1993)" = Fleiss1993)
  })  
  
  qbinary <- reactive({
    switch(input$dataset,
           "green tea" = FALSE,
           "H.pylori" = TRUE,
           "Fleiss (1993)" = TRUE)
  })  
  
  qmean.sd <- reactive({
    switch(input$dataset,
           "green tea" = TRUE,
           "H.pylori" = FALSE,
           "Fleiss (1993)" = FALSE)
  })
  
  qhigher.is.better <- reactive({
    switch(input$dataset,
           "green tea" = FALSE,
           "H.pylori" = FALSE,
           "Fleiss (1993)" = FALSE)
  })
  
  
  ## render plots
  
  output$plot1 <- renderPlot({
    
    if (input$studies_outlook1 != "initial settings" && substring(input$studies_outlook1, 1, 1) != '-' ){
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
    
    if (input$studies_outlook2 != "initial settings" && substring(input$studies_outlook2, 1, 1) != '-' ){
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
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(inputData(), file, row.names = FALSE)
    }
  )
  
  
  # observeEvent(
  #   input$initial settings,
  #   { output$plot <- renderPlot({forestsens(greentea,
  #                           binary=FALSE,
  #                           mean.sd=TRUE,
  #                           higher.is.better=FALSE,
  #                           outlook = NA)})
  #   }
  # )
  
}

shinyApp(ui = ui, server = server)