library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)

server <- function(input, output) {
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  
  output$"table" <- renderDataTable({
    datatable(sample1)})
  
  output$text_density_plot <- renderText({
    paste("Density plot for", input$variable)})
  
  ### sequences min/max ###
  
  observeEvent(input$do, {
    colmn <- input$column
    n <- input$x
    colmn_quo = parse_quosure(colmn)
    minmax <- as.numeric(input$type)
    #    print(colmn)
    #    print(n)
    #    print(input$column)
    seq_column <- (select(sample1, c("seq", colmn)))
    #    print(seq_column)
    if (input$type==1){
      seq_sorted <- arrange(seq_column, !!colmn_quo)
    } else if (input$type==2){
      seq_sorted <- arrange(seq_column, desc(!!colmn_quo))
    }
    #    print(seq_sorted)
    
    output$"filtered_sequences" <- renderDataTable(
      head(seq_sorted, n)
    )
  }
  )
  
  
  ### density plot ###
  
  y_click <- reactive(input$plot_click$y)
  output$test <- renderText({
    y_click <- input$plot_click$y
    #  print(y_click())
    paste0(as.numeric(y_click))})
  
  
  ##### v.1 #####
  
  # output$density_plot <- reactive({
  #   if (is.null(y_click)==TRUE){
  #     renderPlot({ggplot(sample1, aes(x=ben_WT)) +
  #                 geom_density() + coord_flip()
  #                 + geom_vline(xintercept = 0)})
  #   } else {renderPlot({ggplot(sample1, aes(x=ben_WT)) +
  #         geom_density() + coord_flip()
  #       + geom_vline(xintercept = y_click() )})
  #   }
  #     })
  
  
  
  ##### v.2 #####
  
  # output$density_plot <- reactive({
  #     renderPlot({ggplot(sample1, aes(x=input$variable)) +
  #         geom_density() + coord_flip() +
  #       if (is.null(y_click)==TRUE){
  #       geom_vline(xintercept = 0)}
  #    else { geom_vline(xintercept = y_click())}
  #   })
  # })
  
  
  
  ##### v.3 #####
  
  if (is.null(y_click)==TRUE){
    output$density_plot <- renderPlot({ggplot(sample1, aes(x=ben_WT)) +
        geom_density() + coord_flip() + geom_vline(xintercept=0)}) }

  else {output$density_plot <- renderPlot({ggplot(sample1, aes(x=ben_WT)) +
      geom_density() + coord_flip() + geom_vline(xintercept=as.numeric(y_click()))})
  }
  
}