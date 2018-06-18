library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)

server <- function(input, output) {
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  output[["table"]] <- renderDataTable({
    datatable(sample1)
  })
  
  output[["text_density_plot"]] <- renderText({
    paste("Density plot for", 
          input[["variable"]])
  })
  
  output[["density_plot"]] <- renderPlot({
    ggplot(sample1, 
           aes(x = sample1[, input[["variable"]]])) + 
           labs(x = "value", y = "density") + 
           geom_density() + 
           coord_flip() + 
           geom_vline(xintercept = coord[["y"]])
    })
  
  ### density plot click ###
  
  coord <- reactiveValues(x = 0, y = 0)
  coord[["y"]] <- 0
  observeEvent(input[["plot_click"]], { 
    coord[["y"]] <- input[["plot_click"]][["y"]]
    coord[["x"]] <- input[["plot_click"]][["x"]]
    })
  
  output[["test"]] <- renderText({
    paste0("Threshold position: ", 
           as.numeric(coord[["y"]]))
    })   
  
  ### filtering ###
  
  observeEvent(input[["filter"]], {
    if (is.null(input[["method"]])) { 
      showModal(modalDialog(
        title = "Error:",
        "Please select the method of filtering.",
        easyClose = TRUE))
    } else {
      colmn <- input[["variable"]]
      colmn_quo = parse_quosure(colmn)
      seq_column <- (select(sample1, c(seq, !!colmn_quo)))
      seq_sorted <- arrange(seq_column, !!colmn_quo)
      colnames(seq_sorted) <- c("sequence", "values")
      
      if (input[["method"]] == "Density plot") { 
        if (input[["thresh_button"]] == 1) {
          seq_out <- filter(seq_sorted, 
                            values < coord[["y"]])
        } else if (input[["thresh_button"]] == 2) {
          seq_out <- filter(seq_sorted, 
                            values > coord[["y"]])
        }
        
      } else if (input[["method"]] == "Minimum/maximum values") {
        n <- input[["x"]]
        if (input[["type"]] == 1) {
          seq_out <- head(seq_sorted, n)
        } else if (input[["type"]] == 2) {
          seq_out <- head(arrange(seq_sorted, desc(values)), n)
        }
      }
      
      output[["table"]] <- renderDataTable({
        datatable(seq_out)
        })
    }
  })
  
  observeEvent(input[["reset"]], {
    output[["table"]] <- renderDataTable({
      datatable(sample1)
      })
  })
}