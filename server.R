library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)

server <- function(input, output) {
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  
  output[["table"]] <- renderDataTable({
    datatable(sample1)})
  
  output[["text_density_plot"]] <- renderText({
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
    
    output[["filtered_sequences"]] <- renderDataTable(
      head(seq_sorted, n)
    )
  }
  )
  
  
  ### density plot ###
  
  coord <- reactiveValues(x=0, y=0)
  
  observeEvent(input$plot_click,  {
    coord$y <- input$plot_click$y
    coord$x <- input$plot_click$x
    colm <- (parse_quosure(input$variable))
    seq_c <- select(sample1, c(seq, !!colm))
    seq_a <- arrange(seq_c, !!colm)
    if(input$thresh_button==1){
      sec_f <- filter(seq_a, !!colm<coord$y)
      output[["table"]] <- renderDataTable({datatable(sec_f)})
    }else if (input$thresh_button==2){
      sec_f <- filter(seq_a, !!colm>coord$y)
      output[["table"]] <- renderDataTable({datatable(sec_f)})
    }
  })
  
  output[["test"]] <- renderText({
    paste0("Filter position: ", as.numeric(coord$y))})
  
  output[["density_plot"]] <- renderPlot({ggplot(sample1, aes(x=sample1[, input$variable])) + labs(x="value", y="density") + 
      geom_density() + coord_flip() + geom_vline(xintercept=coord$y)})
  
  
}