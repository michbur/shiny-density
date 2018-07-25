library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)

server <- function(input, output) {
  
  counting <- function(dt) {
    prot_grouped <- arrange(group_by(dt, prot_id, type))
    counts <- (summarise(prot_grouped, count = n()))
    types <- levels(prot_grouped[["type"]])
    columns <- c("prot_id", types)
    data <- data.frame(prot_id = levels(prot_grouped[["prot_id"]]))
    for (i in types) {
      filtered <- filter(counts, type == i)
      vals <- select(filtered, c(prot_id, count))
      colnames(vals)[2] <- i
      df2 <- data.frame(vals)
      df3 <- merge.data.frame (data, df2, by = "prot_id", all = TRUE)
      data <- df3
    }
    good_data <- data[rowSums(is.na(data)) != (ncol(data)-1),]
    return(good_data)
  }
  
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  proteins <- counting(sample1)
  
  output[["table"]] <- renderDataTable({
    datatable(sample1)
  })
  
  output[["proteins"]] <- renderDataTable({
    datatable(proteins)
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
 #     seq_selected <- select(seq_filtered, c(seq, prot_id, InROPE))
      seq_sorted <- arrange(seq_filtered, InROPE)
      
      if (input[["method"]] == "Density plot") { 
        if (input[["thresh_button"]] == 1) {
          seq_out <- filter(seq_sorted,
                            InROPE < coord[["y"]])
          seq_prot <- filter(sample1,
                             InROPE < coord[["y"]])
          
        } else if (input[["thresh_button"]] == 2) {
          seq_out <- filter(seq_sorted, 
                            InROPE > coord[["y"]])
          seq_prot <- filter(sample1,
                             InROPE < coord[["y"]])
        }
        peptides <- counting(seq_prot)
        
      } else if (input[["method"]] == "Minimum/maximum values") {
        n <- input[["x"]]
        if (input[["type"]] == 1) {
          seq_out <- head(seq_sorted, n)
          seq_prot <- head(arrange(sample1, InROPE), n)
        } else if (input[["type"]] == 2) {
          seq_out <- head(arrange(seq_sorted, desc(InROPE)), n)
          seq_prot <- head(arrange(sample1, desc(InROPE)), n)
        }
        peptides <- counting(seq_prot)
      }
      
      output[["table"]] <- renderDataTable({
        datatable(seq_out)
      })
      output[["proteins"]] <- renderDataTable({
        datatable(peptides)
      })
    }
  })
  
  observeEvent(input[["reset"]], {
    output[["table"]] <- renderDataTable({
      datatable(sample1)
    })
    output[["proteins"]] <- renderDataTable({
      datatable(proteins)
    })
  })
}