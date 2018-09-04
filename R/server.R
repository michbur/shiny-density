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
    summarise(count = length(Sequence)) %>% 
    dcast(prot_id ~ type)
}

server <- function(input, output) {
  
  # Data input
  if(Sys.info()[["nodename"]] == "amyloid")
    load("/home/michal/Dropbox/PepArray_results/2018-06-07/full_data.RData")
    #full_data <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
  if(Sys.info()[["nodename"]] == "LENOVO")
    full_data <- read.csv("C:/Users/Kaede/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
  
  # Extract phenotypes from data
  phenotypes <- levels(full_data[["type"]]) 
  
  # Count peptides in proteins
  proteins <- counting(full_data)
  
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
    my_DT(full_data) %>% 
      formatRound(columns = c("muDiff", "effSz", "InROPE"),
                  digits = 4)
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
    sample_chosen <- filter(full_data, type == input[["var"]])
    
    output[["density_plot"]] <- renderPlot({
      ggplot(sample_chosen,
            aes(x = InROPE)) + 
            labs(x = "InROPE values", y = "count of peptides") + 
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
      seq_sorted <- full_data %>%
        filter(type == input[["var"]]) %>%
        arrange(InROPE)
      # Results of filtering by density plot method:
      # seq_out - selected data for 'Data' DT
      # seq_prot - selected data for peptide counting
      
      ### Density plot - peptides 
      
      if (input[["method"]] == "Density plot" & input[["level1"]] == "Peptides") { 
        # values lower than threshold selected 
        if (input[["thresh_button"]] == 1) {
          seq_out <- filter(seq_sorted,
                            InROPE < coord[["y"]])
          seq_prot <- filter(full_data,
                             InROPE < coord[["y"]])
          seq_map <- full_data %>%
            filter(Sequence %in% seq_out[["Sequence"]])
        # values higher than threshold selected  
        } else if (input[["thresh_button"]] == 2) {
          seq_out <- filter(seq_sorted, 
                            InROPE > coord[["y"]])
          seq_prot <- filter(full_data,
                             InROPE < coord[["y"]])
          seq_map <- full_data %>%
            filter(Sequence %in% seq_out[["Sequence"]])
        }
        peptides <- counting(seq_prot)
      
      ### Density plot - proteins
        
        } else if (input[["method"]] == "Density plot" & input[["level1"]] == "Proteins") {
      
        # values lower than threshold selected 
        if (input[["thresh_button"]] == 1) {
          seq_out <- filter(seq_sorted,
                            InROPE < coord[["y"]])
          seq_prot <- full_data %>%
            filter(InROPE < coord[["y"]])
          peptides_initial <- counting(seq_prot)
          peptides_selected <- peptides_initial %>%
            filter_(paste(input[["var"]], ">=", input[["n_pept1"]]))
          peptides <- peptides_initial %>%
            filter_at(.vars=2:ncol(peptides_initial), any_vars(. >= input[["n_pept1"]]))
          seq_map <- full_data %>%
            filter(prot_id %in% peptides_selected[["prot_id"]])
          seq_out <- seq_map %>%
            filter(type == input[["var"]])
       
        # values higher than threshold selected  
        } else if (input[["thresh_button"]] == 2){
          seq_out <- filter(seq_sorted,
                            InROPE > coord[["y"]])
          seq_prot <- full_data %>%
            filter(InROPE > coord[["y"]])
          peptides_initial <- counting(seq_prot)
          peptides_selected <- peptides_initial %>%
            filter_(paste(input[["var"]], ">=", input[["n_pept1"]]))
          peptides <- peptides_initial %>%
            filter_at(.vars=2:ncol(peptides_initial), any_vars(. >= input[["n_pept1"]]))
          seq_map <- full_data %>%
            filter(prot_id %in% peptides_selected[["prot_id"]])
          seq_out <- seq_map %>% 
            filter(type == input[["var"]])
        }

        # Results of filtering by min/max values:
        # Min/max - peptides  
          
      } else if (input[["method"]] == "Minimum/maximum values" & input[["level2"]] == "Peptides") {
        # how many values
        n <- input[["x"]]
        # minimum values selected
        if (input[["type"]] == 1) {
          seq_out <- head(arrange(full_data, InROPE), n)
          seq_prot <- full_data %>%
            filter(Sequence %in% seq_out[["Sequence"]])
          seq_map <- seq_prot
        # maximum values selected
        } else if (input[["type"]] == 2) {
          seq_out <- head(arrange(full_data, desc(InROPE)), n)
          seq_prot <- full_data %>%
            filter(Sequence %in% seq_out[["Sequence"]])
          seq_map <- seq_prot
          }
        peptides <- counting(seq_prot)
        
        # Min/max - proteins
        
      } else if (input[["method"]] == "Minimum/maximum values" & input[["level2"]] == "Proteins") {
        n <- input[["x"]]
        # minimum values selected
        if (input[["type"]] == 1) {
          seq_sorted <- arrange(full_data, InROPE)
          peptides_initial <- counting(seq_sorted)
          peptides <- peptides_initial %>%
            filter_at(.vars=2:ncol(peptides_initial), any_vars(. >= input[["n_pept2"]])) %>%
            group_by(prot_id) %>%
            head(n)
          seq_out <- seq_sorted %>%
            filter(prot_id %in% peptides[["prot_id"]])
          seq_map <- seq_out
        # maximum values selected
        } else if (input[["type"]] == 2) {
            seq_sorted <- arrange(full_data, desc(InROPE))
            peptides_initial <- counting(seq_sorted)
            peptides <- peptides_initial %>%
              filter_at(.vars=2:ncol(peptides_initial), any_vars(. >= input[["n_pept2"]])) %>%
              group_by(prot_id) %>%
              head(n)
            seq_out <- seq_sorted %>%
              filter(prot_id %in% peptides[["prot_id"]])
            seq_map <- seq_out
        }
      }
      # Outputs with results of filtering:
      output[["table"]] <- renderDataTable({
        my_DT(seq_out) %>% 
          formatRound(columns = c("muDiff", "effSz", "InROPE"),
                      digits = 4)
      })
      output[["proteins"]] <- renderDataTable({
        my_DT(peptides)
  #            container = sketch, 
  #            colnames = FALSE)
      })
      # Display number of selected peptides
      output[["n_selected"]] <- renderText({
        paste0("Selected ", nrow(seq_out), " peptides.")
      })
      # Render heatmap
      output[["heatmap"]] <- renderPlot({
        ggplot(seq_map,
               aes(x=Sequence, y=type)) +
          geom_tile(aes(fill=InROPE), color = "black") +
          coord_flip() 
      })
        
    }
  })
  
  
  # Reset filtering - return to initial outputs
  observeEvent(input[["reset"]], {
    output[["table"]] <- renderDataTable({
      my_DT(full_data) %>% 
        formatRound(columns = c("muDiff", "effSz", "InROPE"),
                    digits = 4)
    })
    output[["proteins"]] <- renderDataTable({
      my_DT(proteins,
            container = sketch,
            colnames = FALSE)
    })
  })
}