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

    ## left hand sidebar 
    column(4,
           
      wellPanel(
        
        fileInput("file1", "Choose data text file: ",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values, text/plain",
                             ".csv",
                             ".txt")),
        
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),

        tags$h6("Note that the file should have specific column headers, as seen in one of the data set examples (which you can download below).")
      ),      
           
      ## forestsens() parameters
      wellPanel(
        tags$h4("Plot settings: "),
        checkboxInput("binary", "binary (Outcome is dichotomous)", TRUE),
        checkboxInput("msd", "mean.sd (Data includes means and standard deviations)", TRUE),
        checkboxInput("higher_better", "higher.is.better (Higher outcome or counts desired)", TRUE),
        
        tags$h6("If you get an error message for the plots, first try changing these settings (so that they are appropriate for the input data set).")
      ),          
           
      ## examples of data sets
      wellPanel(
        tags$h4("Data set examples: "),
        
        ## select existing dataset
        selectInput("dataset", 
                    label = "",
                    choices = c("Fleiss1993", "greentea", "Hpylori"),
                    selected = "greentea",
                    width = '100%'),       
  
        ## Choose dataset to download
        downloadButton("downloadData", "Download as a *.txt file"),
        # Select separator 
        radioButtons("seperator", "Select separator for the text file to be downloaded: ",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ",")
      ),

      wellPanel(
        tags$a(href='https://noory.shinyapps.io/samurai_example', 'Click here to see a demo with the data set examples')
      ),
      
      ## References        
      wellPanel(           
        tags$h4('References:'),
        tags$a(href='https://doi.org/10.1186/2046-4053-3-27', 'Systematic Reviews (2014) 3:27'),
        br(),
        tags$a(href='https://cran.r-project.org/web/packages/SAMURAI/index.html', 'SAMURAI, an R package')
      )
      
    ),
    
    ## right hand side 
    column(8,
      tags$h4("Choose outlook of all unpublished studies"),
      column(6,
        ## select box 1
        selectInput("studies_outlook1", 
                   label = "in Plot 1", 
                   choices = list("initial settings",
                                  "-- --",
                                  "very positive", "positive", "no effect", "negative", "very negative", 
                                  "-- based on published studies' effect size and CL --",
                                  "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL"), 
                   selected = "no effect",
                   width = '100%')
      ),
      column(6,
        ## select box 2
        selectInput("studies_outlook2", 
                   label = "in Plot 2", 
                   choices = list("initial settings",
                                  "-- --",
                                  "very positive", "positive", "no effect", "negative", "very negative", 
                                  "-- based on published studies' effect size and CL --",
                                  "very positive CL", "positive CL", "current effect", "negative CL", "very negative CL"), 
                   selected = "current effect",
                   width = '100%')
        )
      ),  
  
      ## display outputs 
      column(8,
        tags$h4("Output"),           
        tabsetPanel(type = "tabs",
                    tabPanel("Data set", tableOutput(outputId = "dscontents")),
                    tabPanel("Plot 1", plotOutput(outputId = "plot1")),
                    tabPanel("Plot 2", plotOutput(outputId = "plot2"))
      )
      
    )
    
  )
)

## Define server logic for app
server <- function(input, output) {
  
  ## Reactive values for dataset selected
  
  #inputData <- reactive({ input$file1 })
  inputData <- reactive({     
    df <- read.table(input$file1$datapath,
                     header = TRUE,
                     sep = input$sep)
    return(df)
  })
  
  downloadData <- reactive({
    switch(input$dataset,
           "greentea" = greentea,
           "Hpylori" = Hpylori,
           "Fleiss1993" = Fleiss1993)
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
               binary = input$binary, 
               mean.sd = input$msd, 
               higher.is.better = input$higher_better,
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
               binary = input$binary, 
               mean.sd = input$msd, 
               higher.is.better = input$higher_better,
               outlook = uoutlook)
    
  })
  
  # Downloadable text file of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".txt", sep = "")
    },
    content = function(file) {
      write.table(downloadData(), file, row.names = FALSE, quote = FALSE, sep = input$seperator)
    }
  )
  
  output$dscontents <- renderTable({
    req(input$file1)
    
    df <- read.table(input$file1$datapath,
                   header = TRUE,
                   sep = input$sep
#                   quote = input$quote
                   )
    return(df)
  })
    
}

shinyApp(ui = ui, server = server)