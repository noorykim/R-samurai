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
                   choices = c("Fleiss (1993)", "green tea", "H. pylori"),
                   selected = "green tea",
                   width = '100%'),       
        
        ## Choose dataset to download
        downloadButton("downloadData", "Download as a *.txt file"),
        tags$br(),
        # Select separator 
        radioButtons("seperator", "Select separator for the text file to be downloaded: ",
                    choices = c(Comma = ",",
                                Semicolon = ";",
                                Tab = "\t"),
                    selected = ",")
      ),

      ## Select outlooks           
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
        tags$a(href='https://noory.shinyapps.io/samurai/', 'Click here to upload your own dataset (coming soon)')
      ),
      
      ## References
      wellPanel(           
        tags$h4('References for software and methodology:'),
        tags$a(href='https://doi.org/10.1186/2046-4053-3-27', 'Systematic Reviews (2014) 3:27'),
        br(),
        tags$a(href='https://cran.r-project.org/web/packages/SAMURAI/index.html', 'SAMURAI, an R package')
      )
      
    ),
    
    ## display outputs 
    column(8,
           
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                 tabPanel("Click here for Plot 1", plotOutput(outputId = "plot1")),
                 tabPanel("Plot 2", plotOutput(outputId = "plot2")),
                 tabPanel("Data set (with initial settings)", tableOutput(outputId = "dscontents")) 
      ),

      wellPanel(
        tags$h4('Data set description:'),
        textOutput(outputId = "ds_description")
      ),
      
      wellPanel(
        tags$h4('Syntax:'),
        textOutput(outputId = "syntax")
      ),
      
      
      wellPanel(
        tags$h4('Data set reference:'),
        textOutput(outputId = "ds_reference"),
        
        tags$h4('Webpage:'),
        textOutput(outputId = "ds_reference_link")
      )
      

    )  
    
  )
)

## Define server logic for app
server <- function(input, output) {
  
  ## Reactive values for dataset selected
  
  inputData <- reactive({
    switch(input$dataset,
           "green tea" = greentea,
           "H. pylori" = Hpylori,
           "Fleiss (1993)" = Fleiss1993)
  })  

  # inputDatasetName <- reactive({
  #   switch(input$dataset,
  #          "green tea" = "greentea",
  #          "H. pylori" = "Hpylori",
  #          "Fleiss (1993)" = "Fleiss1993")
  # })   
    
  qbinary <- reactive({
    switch(input$dataset,
           "green tea" = FALSE,
           "H. pylori" = TRUE,
           "Fleiss (1993)" = TRUE)
  })  
  
  qmean.sd <- reactive({
    switch(input$dataset,
           "green tea" = TRUE,
           "H. pylori" = FALSE,
           "Fleiss (1993)" = FALSE)
  })
  
  qhigher.is.better <- reactive({
    switch(input$dataset,
           "green tea" = FALSE,
           "H. pylori" = FALSE,
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
  
  # Downloadable text file of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".txt", sep = "")
    },
    content = function(file) {
      write.table(inputData(), file, row.names = FALSE, quote = FALSE, sep = input$seperator)
    }
  )
  
  ## to show data set
  output$dscontents <- renderTable({
    inputData()
  })  
  
  ## to show syntax 
  output$syntax <-  reactive({
    switch(input$dataset,
           "green tea" = "> forestsens(greentea, binary=FALSE, mean.sd=TRUE, higher.is.better=FALSE)",
           "H. pylori" = "> forestsens(table=Hpylori, binary=TRUE, mean.sd=FALSE, higher.is.better=FALSE)",
           "Fleiss (1993)" = "> forestsens(table=Fleiss1993, binary=TRUE, mean.sd=FALSE, higher.is.better=FALSE)" )
  })

  ## to show description of data set
  output$ds_description <-  reactive({
    switch(input$dataset,
           "green tea" = "Randomized clinical trials of at least 12 weeks duration assessing the effect of green tea consumption on weight loss. The standardized mean differences reflect changes in weight (kg).",
           "H. pylori" = "Randomized clinical trials comparing duodenal ulcer acute healing among (1) patients on ulcer healing drug + Helicobacter pylori eradication therapy vs. (2) patients ulcer healing drug alone. The event counts represent the numbers of patients not healed.",
           "Fleiss (1993)" = "This data set originally included 7 published placebo-controlled randomized studies on the effect of aspirin in preventing death after myocardial infarction. The defined binary outcome event is death. 2 (fictional) unpublished studies have been added." )
  })  

  ## to show description of data set
  output$ds_reference <-  reactive({
    switch(input$dataset,
           "green tea" = "Jurgens TM, Whelan AM, Killian L, Doucette S, Kirk S, Foy E. (2012) Cochrane Database of Systematic Reviews. Art. No.: CD008650. DOI: 10.1002/14651858.CD008650.pub2.",
           "H. pylori" = "Ford AC, Delaney B, Forman D, Moayyedi P. (2006) Cochrane Database of Systematic Reviews. Art No.: CD003840. DOI: 10.1002/14651858.CD003840.pub4.",
           "Fleiss (1993)" = "Fleiss, JL. (1993) Stat Methods Med Res. 2(2): 121-45." )
  })  
  
  output$ds_reference_link <-  reactive({
    switch(input$dataset,
           "green tea" = "http://www.cochrane.org/CD008650/ENDOC_green-tea-for-weight-loss-and-weight-maintenance-in-overweight-or-obese-adults",
           "H. pylori" = "http://onlinelibrary.wiley.com/doi/10.1002/14651858.CD003840.pub4/abstract",
           "Fleiss (1993)" = "http://journals.sagepub.com/doi/abs/10.1177/096228029300200202" )
  })    
  
  
      
}

shinyApp(ui = ui, server = server)