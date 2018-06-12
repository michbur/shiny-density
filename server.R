library(shiny)
library(DT)
library(ggplot2)

server <- function(input, output) {
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  output[["table"]] <- renderDataTable({
    datatable(sample1)})
  
  output$text_density_plot <- renderText({
    paste("Density plot for", input$variable)})
  
  
  # data <- reactive(sample1$input$variable)
  # output[["density_plot"]] <- reactive({ggplot(sample1, 
  #                                              aes(x=data) + geom_density() + coord_flip())})
  
  
  
  # output[["density_plot"]] <- reactive({ggplot(sample1, 
  #                                              aes(x=input$variable) + geom_density() + coord_flip())})
  
  
  
  # if (input$variable == "ben_WT"){
  #   output[["density_plot"]] <- reactive({ggplot(sample1, 
  #                                                  aes(x=ben_WT) + geom_density() + coord_flip())})
  # } else if (input$variable == "ben_mal"){
  #     output[["density_plot"]] <- reactive({ggplot(sample1, 
  #                                                    aes(x=ben_mal)) + geom_density() + coord_flip()})
  # } else if (input$variable == "mal_WT"){
  #    output[["density_plot"]] <- reactive({ggplot(sample1, 
  #                                                   aes(x=mal_WT)) + geom_density() + coord_flip()})
  #   }
  
}