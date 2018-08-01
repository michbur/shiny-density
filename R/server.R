library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)
library(reshape2)

options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

#' Datatable output
#'
#' Create a HTML datatable widget containing data with export buttons using library DT.
#' 
#' @param x Matrix or a data frame containing data.
#' @param ... Optional arguments to \code{datatable}.
#' @return HTML datatable widget displaying data.
#' @examples 
#' my_DT(df)
#' my_DT(matrix, colnames = FALSE)
#' @export
my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE)

#' Counting peptides
#' 
#' Count the occurences of relevant peptides for each protein.
#' 
#' @param dt A data table or data frame.
#' @return Data table containing counts of peptides for each protein.
#' @examples
#' counting(data)
#' counting(sample)
#' @export
counting <- function(dt) {
  dt %>% 
    group_by(prot_id, type) %>% 
    summarise(count = length(seq)) %>% 
    dcast(prot_id ~ type)
}

server <- function(input, output) {
  
  sample1 <- read.csv("https://raw.githubusercontent.com/michbur/shiny-density/master/data/sample1.csv")
  
  phenotypes <- levels(sample1[["type"]]) 
  
  proteins <- counting(sample1)
  
  sketch = htmltools::withTags(table(
    class = 'display',
    thead(
      tr(
        th(rowspan = 2, colspan = 1, 'Protein ID'),
        th(colspan = as.numeric(length(phenotypes)), 'Count of peptides')
      ),
      tr(
        lapply(rep(phenotypes, 1), th)
      )
    )
  ))
  
  output[["table"]] <- renderDataTable({
    my_DT(sample1)
  })
  
  output[["proteins"]] <- renderDataTable({
    my_DT(proteins, 
              container = sketch, 
              colnames = FALSE)
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
            labs(x = "value", y = "count") + 
            #geom_histogram(binwidth = 2.5) +
            geom_density(aes(y = ..count..)) +
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
      seq_sorted<- sample1 %>%
        filter(type == input[["var"]]) %>%
        arrange(InROPE)
      
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
        my_DT(seq_out)
      })
      output[["proteins"]] <- renderDataTable({
        my_DT(peptides, 
                  container = sketch, 
                  colnames = FALSE)
      })
    }
  })
  
  observeEvent(input[["reset"]], {
    output[["table"]] <- renderDataTable({
      my_DT(sample1)
    })
    output[["proteins"]] <- renderDataTable({
      my_DT(proteins,
                container = sketch, 
                colnames = FALSE)
    })
  })
}