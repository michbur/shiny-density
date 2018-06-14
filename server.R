library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)

server <- function(input, output) {
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")

    
  output[["table"]] <- renderDataTable({
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
    
    output[["filtered_sequences"]] <- renderDataTable(
       head(seq_sorted, n)
    )
  }
)
  
  
  ### density plot ###
  
  eventReactive(input$variable, {
    colmn_var = parse_quosure(input$variable)
  output$density_plot<- reactivePlot({
    ggplot(sample1, aes(x=(!!colmn_var))) + geom_density() + coord_flip()
   })
}
)
}
