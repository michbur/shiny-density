library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(rlang)
library(reshape2)

# DT output options
options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE)


# Function for counting peptides in proteins and phenotypes
counting <- function(dt) {
  dt %>% 
    group_by(prot_id, type) %>% 
    summarise(count = length(seq)) %>% 
    dcast(prot_id ~ type)
}

server <- function(input, output) {
  
  # Data input
  if(Sys.info()[["nodename"]] == "amyloid")
    sample1 <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
  if(Sys.info()[["nodename"]] == "LENOVO")
    sample1 <- read.csv("C:/Users/Kaede/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
  
  # Extract phenotypes from data
  phenotypes <- levels(sample1[["type"]]) 
  
  # Count peptides in proteins
  proteins <- counting(sample1)
  
  # Container for proteins data table
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
  
  # Initial outputs:
  # all input data
  output[["table"]] <- renderDataTable({
    my_DT(sample1)
  })
  # counts of peptides in proteins
  output[["proteins"]] <- renderDataTable({
    my_DT(proteins, 
              container = sketch, 
              colnames = FALSE)
  })
  # density plot for the first phenotype
  output[["text_density_plot"]] <- renderText({
    paste("Density plot for", 
          input[["var"]])
  })
  

  # Choice of a phenotype and density plot output 
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
  
  # Density plot click - set threshold and display its position
  coord <- reactiveValues(x = 0, y = 0)
  coord[["y"]] <- 0
  observeEvent(input[["plot_click"]], { 
    coord[["y"]] <- input[["plot_click"]][["y"]]
    coord[["x"]] <- input[["plot_click"]][["x"]]
  })
  
  output[["test"]] <- renderText({
    paste0("Threshold position: ", 
           round(as.numeric(coord[["y"]]), 4))
  })   
  
  ### Filtering ###
  
  observeEvent(input[["filter"]], {
    # None of the methods selected - display modal with error
    if (is.null(input[["method"]])) { 
      showModal(modalDialog(
        title = "Error:",
        "Please select the method of filtering.",
        easyClose = TRUE))
    } else {
      # Sorting selected data
      seq_sorted <- sample1 %>%
        filter(type == input[["var"]]) %>%
        arrange(InROPE)
      # Results of filtering by density plot method:
      # seq_out - selected data for 'Data' DT
      # seq_prot - selected data for peptide counting
      if (input[["method"]] == "Density plot") { 
        # values lower than threshold selected 
        if (input[["thresh_button"]] == 1) {
          seq_out <- filter(seq_sorted,
                            InROPE < coord[["y"]])
          seq_prot <- filter(sample1,
                             InROPE < coord[["y"]])
        # values higher than threshold selected  
        } else if (input[["thresh_button"]] == 2) {
          seq_out <- filter(seq_sorted, 
                            InROPE > coord[["y"]])
          seq_prot <- filter(sample1,
                             InROPE < coord[["y"]])
        }
        peptides <- counting(seq_prot)
        
        # Results of filtering by min/max values:
      } else if (input[["method"]] == "Minimum/maximum values") {
        # how many values
        n <- input[["x"]]
        # minimum values selected
        if (input[["type"]] == 1) {
          seq_out <- head(seq_sorted, n)
          seq_prot <- head(arrange(sample1, InROPE), n)
        # maximum values selected
        } else if (input[["type"]] == 2) {
          seq_out <- head(arrange(seq_sorted, desc(InROPE)), n)
          seq_prot <- head(arrange(sample1, desc(InROPE)), n)
        }
        peptides <- counting(seq_prot)

        # Results of filtering on protein level:
      } else if (input[["method"]] == "Peptides on protein level") {
        seq_prot <- filter(sample1, type == input[["pept_types"]])  
        peptides <- counting(sample1) %>%
                        select(prot_id, input[["pept_types"]])
        colnames(peptides) <- c("prot_id", "count of peptides")
        peptides_selected <- peptides %>%
                              filter(count >= input[["n_pept"]])
        seq_out <- inner_join(peptides_selected, seq_prot)
        }
      
      
      # Outputs with results of filtering:
      output[["table"]] <- renderDataTable({
        my_DT(seq_out)
      })
      output[["proteins"]] <- renderDataTable({
        my_DT(peptides, 
              container = sketch, 
              colnames = FALSE)
      })
      # Display number of selected peptides
      output[["n_selected"]] <- renderText({
        paste0("Selected ", nrow(seq_out), " peptides.")
      })
    }
  })
  
  # Reset filtering - return to initial outputs
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