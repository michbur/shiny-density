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


# Add protein count to proteins table
prot_number <- function(data) {
  numbers <- data %>%
              group_by(prot_id) %>%
              summarise(prot_number = length(Sequence))
  data_out <- left_join(counting(data), numbers)
  return(data_out)
}

server <- function(input, output) {
  
  # Data input
  if(Sys.info()[["nodename"]] == "amyloid")
    load("/home/michal/Dropbox/PepArray_results/2018-06-07/full_data.RData")
    #full_data <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
  if(Sys.info()[["nodename"]] == "LENOVO")
    #full_data <- read.csv("C:/Users/Kaede/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
    load("C:/Users/Kaede/Dropbox/PepArray_results/2018-06-07/full_data.RData")
  
  # Divide full data
  data <- select(full_data, -gene_name) %>% 
    arrange(prot_id)
  genes <- select(full_data, c(prot_id, gene_name)) %>%
    unique
  
  
  # Define parameter for filtering
  param <- reactiveValues()
  observe({
    param[["parameter"]] <- (input[["parameter"]])
    param[["quoted"]] <- shQuote(param[["parameter"]])
   })

  
  # Choice of parameter for filtering and a phenotype and density plot output 
  observeEvent(input[["var"]], {
      sample_chosen <- filter(data, type == input[["var"]])
      
      output[["density_plot"]] <- renderPlot({
        ggplot(sample_chosen, aes(x = sample_chosen[[param[["parameter"]]]])) + 
          labs(x = "values", y = "count of peptides") + 
          #geom_histogram(binwidth = 2.5) +
          geom_density(aes(y = ..count..)) +
          geom_density() + 
          coord_flip() + 
          geom_vline(xintercept = coord[["y"]]) +
          theme_bw()
      })
    })
  
  # output[["testing"]] <- renderText({
  #   paste(param[["parameter"]], param[["quoted"]])
  # })
  
  # Density plot click - set threshold and display its position
  coord <- reactiveValues(x = 0, y = 0)
  observeEvent(input[["parameter"]], {
  if (input[["parameter"]] == "InROPE") {
  coord[["y"]] <- 50    
  } else if (input[["parameter"]] == "effSz") {
    coord[["y"]] <- 0
  }
  })

  observeEvent(input[["plot_click"]], { 
    coord[["y"]] <- input[["plot_click"]][["y"]]
    coord[["x"]] <- input[["plot_click"]][["x"]]
  })
  
  output[["test"]] <- renderText({
    paste0("Threshold position: ", 
           round(as.numeric(coord[["y"]]), 4))
  })   
  
  # Calculate medians of InROPE/effect size
  medians <- reactive({
    data %>%
    group_by(prot_id) %>% 
    summarise(median = median(case_when(
      input[["parameter"]] == "InROPE" ~ InROPE,
      input[["parameter"]] == "effSz" ~ effSz)
    ))
  })
  
  
  ### Filtering functions ------------------------------------------------------------
  
  ### Filtering by density plot 
  
  filter_density_plot <- function(data) {
    # lower or higher than threshold
    if (input[["thresh_button"]] == 1) {
      sign <- "<"
    } else if (input[["thresh_button"]] == 2) {
      sign <- ">"
    }
    
    # Peptide level
    if (input[["level1"]] == "Peptides") {
      seq_prot <- data %>%
        filter_(paste(param[["parameter"]], sign, coord[["y"]]))
      seq_out <- seq_prot %>% 
        filter(type == input[["var"]]) %>% 
        arrange(prot_id)
      seq_map <- data %>%
        filter(Sequence %in% seq_out[["Sequence"]])
      peptides <- counting(seq_prot)
    }
    # Protein level
    else if (input[["level1"]] == "Proteins") {
      seq_for_out <- data %>%
        filter(type == input[["var"]]) %>%
        filter_(paste(param[["parameter"]], sign, coord[["y"]])) # peptides below/above the threshold, for 'Data' table, only one phenotype
      seq_for_prot <- data %>% filter_(paste(param[["parameter"]], "<", coord[["y"]]))
      peptides_for_out <- counting(seq_for_prot) %>%
        filter_(paste(input[["var"]], ">=", as.numeric(input[["n_pept1"]]))) # selected peptides for 'Data' table, only proteins of one phenotype with given peptide count
      peptides_for_prot <- counting(seq_for_prot)
      peptides <- peptides_for_prot %>%
        filter_at(.vars=2:ncol(peptides_for_prot), any_vars(. >= input[["n_pept1"]])) # selected peptide counts for 'Proteins' table, all phenotypes with given peptide count
      seq_out <- filter(seq_for_out, prot_id %in% peptides_for_out[["prot_id"]]) %>% 
        arrange(prot_id)# output for 'Data' table, all selected peptide sequences from selected proteins and given phenotype
      seq_map <- filter(data, Sequence %in% seq_out[["Sequence"]]) # all the peptide sequences from selected proteins + the same sequences from other phenotypes
    }
    
    return(list(seq_out, seq_map, peptides))
  }
  
  ### Filtering by minimum/maximum values
  
  filter_min_max <- function(data) {
    n <- input[["x"]] # how many peptides/proteins should be displayed
    # minimum or maximum values?
    sorting <- function(data, data_column) {
      if (input[["type"]] == 1) { # minimum
        arrange(data, data[[data_column]])
      } else if (input[["type"]] == 2) { # maximum
        arrange(data, desc(data[[data_column]]))
      }
    }
    
    
    # Peptide level
    if (input[["level2"]] == "Peptides") {
      seq_out <- sorting(data, param[["parameter"]]) %>%
        head(n) %>% 
        arrange(prot_id) # select n peptides with the lowes/highest InROPE values
      seq_prot <- seq_map <- data %>%
        filter(Sequence %in% seq_out[["Sequence"]]) # select from full data all sequences from seq_out (include all phenotypes)
      peptides <- counting(seq_prot)
    } 
    # Protein level
    else if (input[["level2"]] == "Proteins") {
      peptides_initial <- prot_number(data) # 'Protein' table + peptide count in each protein (for all phenotypes together)
      peptides <- peptides_initial %>%
        filter(prot_number >= input[["n_pept2"]]) %>%
        inner_join(medians(), by="prot_id") %>%
        sorting("median") %>%
        head(n) %>% 
        select(-c(prot_number, median)) # select n of proteins with the lowest/highest InROPE median
      seq_out <- seq_map <- data %>%
        filter(prot_id %in% peptides[["prot_id"]]) %>% 
        arrange(prot_id) # output data for 'Data' table and heatmap, all sequences of peptides in selected proteins
    }
    return(list(seq_out, seq_map, peptides))
  }
  
  ###----------------------------------------------------------------------------------
  
  
  # Extract phenotypes from data
  phenotypes <- levels(data[["type"]]) 
  
  # Count peptides in proteins
  proteins <- reactive({ 
    counting(data) %>%
      inner_join(genes, by="prot_id") %>% 
      inner_join(medians(), by="prot_id")
  })
  
  # Initial outputs:
  # all input data
  output[["table"]] <- renderDataTable({
    my_DT(data) %>% 
      formatRound(columns = c("muDiff", "effSz", "InROPE"),
                  digits = 4)
  })
  # counts of peptides in proteins
  output[["proteins"]] <- renderDataTable({
    my_DT(proteins())
  })
  
  # density plot for the first phenotype
  output[["text_density_plot"]] <- renderText({
    paste("Density plot for", 
          input[["var"]])
  })
  
  
  ### Filtering ###
  
  observeEvent(input[["filter"]], {
     
      ### Density plot 
      
      if (input[["method"]] == "Density plot") { 
        
        seq_out <- filter_density_plot(data)[[1]]
        seq_map <- filter_density_plot(data)[[2]]
        peptides <- filter_density_plot(data)[[3]]
      } 
    
      ### Minimum/maximum values
    
        else if (input[["method"]] == "Minimum/maximum values") {
          
          seq_out <- filter_min_max(data)[[1]]
          seq_map <- filter_min_max(data)[[2]]
          peptides <- filter_min_max(data)[[3]]
        }
        

      # Join peptides with gene names and InROPE medians
      peptides_genes <- reactive({
        inner_join(peptides, genes, by="prot_id") %>%
        inner_join(medians(), by="prot_id")
      })
      
      # Outputs with results of filtering:
      output[["table"]] <- renderDataTable({
        my_DT(seq_out) %>% 
          formatRound(columns = c("muDiff", "effSz", "InROPE"),
                      digits = 4)
      })
      output[["proteins"]] <- renderDataTable({
        my_DT(peptides_genes())
      })
      # Display number of selected peptides
      output[["n_selected"]] <- renderText({
        paste0("Selected ", nrow(seq_out), " peptides.")
      })
      # Render heatmap
      output[["heatmap"]] <- renderPlot({
        ggplot(seq_map, aes(x = type, y = Sequence, fill = case_when(
          input[["parameter"]] == "InROPE" ~ InROPE,
          input[["parameter"]] == "effSz" ~ effSz))) +
          labs(fill = input[["parameter"]]) +
          geom_tile(color = "black") +
          scale_fill_gradient2(midpoint = case_when(
            input[["parameter"]] == "InROPE" ~ 50,
            input[["parameter"]] == "effSz" ~ 0)) +
          theme_bw()
      }, height = 360 + 50*nrow(seq_map))
  })
  
  
  # Reset filtering - return to initial outputs
  observeEvent(input[["reset"]], {
    output[["table"]] <- renderDataTable({
      my_DT(data) %>% 
        formatRound(columns = c("muDiff", "effSz", "InROPE"),
                    digits = 4)
    })
    output[["proteins"]] <- renderDataTable({
      my_DT(proteins())
    })
  })
}