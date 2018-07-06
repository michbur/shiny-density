library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)

server <- function(input, output) {
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  variable = levels(sample1[["type"]])
  
  output[["table"]] <- renderDataTable({
    datatable(sample1)
  })
  
  output[["text_density_plot"]] <- renderText({
    paste("Density plot for", 
          input[["var"]])
  })
  
  observeEvent(input[["var"]], {
    sample_chosen <- filter(sample1, type == input[["var"]])
 
    output[["density_plot"]] <- renderPlot({
      ggplot(sample_chosen,
            aes(x = InROPE)) + 
            labs(x = "value", y = "density") + 
            geom_density() + 
            coord_flip() + 
            geom_vline(xintercept = coord[["y"]])
    }) 
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
      seq_filtered <- filter(sample1, type == input[["var"]])
      seq_selected <- select(seq_filtered, c(seq, prot_id, InROPE))
      seq_sorted <- arrange(seq_selected, InROPE)
      
      if (input[["method"]] == "Density plot") { 
        if (input[["thresh_button"]] == 1) {
          seq_out <- filter(seq_sorted, 
                            InROPE < coord[["y"]])
        } else if (input[["thresh_button"]] == 2) {
          seq_out <- filter(seq_sorted, 
                            InROPE > coord[["y"]])
        }
        
      } else if (input[["method"]] == "Minimum/maximum values") {
        n <- input[["x"]]
        if (input[["type"]] == 1) {
          seq_out <- head(seq_sorted, n)
        } else if (input[["type"]] == 2) {
          seq_out <- head(arrange(seq_sorted, desc(InROPE)), n)
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