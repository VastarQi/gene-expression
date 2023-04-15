yeastdata=read.table("C:/Users/97481/Downloads/spellman.txt", header=T, row.names = 1)
cdcyeast=yeastdata[grepl("cdc15",colnames(yeastdata))]

ui<- fluidPage(
  sidebarLayout(
    sidebarPanel(selectInput('xcol', 'X Variable', dimnames(cdcyeast)[[2]]),
                 selectInput('ycol','y Variable', dimnames(cdcyeast)[[2]]),
                 selectInput('color','Point color',list("Red" = "Red", "Blue" = "Blue",
                                                        "Green" = "#70AD47"), selected = 1)
                        ),
    mainPanel(plotOutput("plot1"))
    )
  

    )
server <- function(input, output,session) {
     
    selectedData <- reactive({cdcyeast[, c(input$xcol, input$ycol)]}) 
   
    output$plot1 <- renderPlot({ par(mar=c(5.1, 4.1,  0, 1))    
      plot(selectedData(), col = input$color, pch = 20, cex 
           = 1) 
}) 
  }


shinyApp(ui=ui,server = server)



