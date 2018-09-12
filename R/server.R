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

# Compute InROPE medians
median_inrope <- function(data) {
  data %>%
    group_by(prot_id) %>%
    summarise(InROPE_median=median(InROPE))
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
        filter_(paste("InROPE", sign, coord[["y"]]))
      seq_out <- seq_prot %>% 
        filter(type == input[["var"]])
      seq_map <- data %>%
        filter(Sequence %in% seq_out[["Sequence"]])
      peptides <- counting(seq_prot)
    }
    # Protein level
    else if (input[["level1"]] == "Proteins") {
      seq_for_out <- data %>%
        filter(type == input[["var"]]) %>%
        filter_(paste("InROPE", sign, coord[["y"]])) # peptides below/above the threshold, for 'Data' table, only one phenotype
      seq_for_prot <- filter(data, InROPE < coord[["y"]])
      peptides_for_out <- counting(seq_for_prot) %>%
        filter_(paste(input[["var"]], ">=", as.numeric(input[["n_pept1"]]))) # selected peptides for 'Data' table, only proteins of one phenotype with given peptide count
      peptides_for_prot <- counting(seq_for_prot)
      peptides <- peptides_for_prot %>%
        filter_at(.vars=2:ncol(peptides_for_prot), any_vars(. >= input[["n_pept1"]])) # selected peptide counts for 'Proteins' table, all phenotypes with given peptide count
      seq_out <- filter(seq_for_out, prot_id %in% peptides_for_out[["prot_id"]]) # output for 'Data' table, all selected peptide sequences from selected proteins and given phenotype
      seq_map <- filter(data, Sequence %in% seq_out[["Sequence"]]) # all the peptide sequences from selected proteins + the same sequences from other phenotypes
    }
    
    return(list(seq_out, seq_map, peptides))
  }
  
  
  
  # Divide full data
  data <- select(full_data, -gene_name)
  genes <- select(full_data, c(prot_id, gene_name)) %>%
    unique
  medians <- median_inrope(data)
  
  # Extract phenotypes from data
  phenotypes <- levels(data[["type"]]) 
  
  # Count peptides in proteins
  proteins <- counting(data) %>%
    inner_join(genes, by="prot_id") %>% 
    inner_join(medians, by="prot_id")
  
  # Initial outputs:
  # all input data
  output[["table"]] <- renderDataTable({
    my_DT(data) %>% 
      formatRound(columns = c("muDiff", "effSz", "InROPE"),
                  digits = 4)
  })
  # counts of peptides in proteins
  output[["proteins"]] <- renderDataTable({
    my_DT(proteins)
  })
  
  # density plot for the first phenotype
  output[["text_density_plot"]] <- renderText({
    paste("Density plot for", 
          input[["var"]])
  })
  

  # Choice of a phenotype and density plot output 
  observeEvent(input[["var"]], {
    sample_chosen <- filter(data, type == input[["var"]])
    
    output[["density_plot"]] <- renderPlot({
      ggplot(sample_chosen, aes(x = InROPE)) + 
        labs(x = "InROPE values", y = "count of peptides") + 
        #geom_histogram(binwidth = 2.5) +
        geom_density(aes(y = ..count..)) +
        geom_density() + 
        coord_flip() + 
        geom_vline(xintercept = coord[["y"]]) +
        theme_bw()
    })
  })
  
  # Density plot click - set threshold and display its position
  coord <- reactiveValues(x = 0, y = 0)
  coord[["y"]] <- 100
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
     
      ### Density plot 
      
      if (input[["method"]] == "Density plot") { 
        
        seq_out <- filter_density_plot(data)[[1]]
        seq_map <- filter_density_plot(data)[[2]]
        peptides <- filter_density_plot(data)[[3]]
        

        # Results of filtering by min/max values:
        # Min/max - peptides  
          
      } else if (input[["method"]] == "Minimum/maximum values" & input[["level2"]] == "Peptides") {
        # how many values
        n <- input[["x"]]
        # minimum values selected
        if (input[["type"]] == 1) {
          seq_out <- head(arrange(data, InROPE), n)
          seq_prot <- data %>%
            filter(Sequence %in% seq_out[["Sequence"]])
          seq_map <- seq_prot
        # maximum values selected
        } else if (input[["type"]] == 2) {
          seq_out <- head(arrange(data, desc(InROPE)), n)
          seq_prot <- data %>%
            filter(Sequence %in% seq_out[["Sequence"]])
          seq_map <- seq_prot
          }
        peptides <- counting(seq_prot)
        
        # Min/max - proteins
        
      } else if (input[["method"]] == "Minimum/maximum values" & input[["level2"]] == "Proteins") {
        n <- input[["x"]]
        # minimum values selected
        if (input[["type"]] == 1) {
          peptides_initial <- prot_number(data)
          peptides <- peptides_initial %>%
            filter(prot_number >= input[["n_pept2"]]) %>%
            inner_join(medians, by="prot_id") %>%
            arrange(InROPE_median) %>%
            head(n) %>%
            select(-c(prot_number, InROPE_median))
          seq_out <- seq_map <- data %>%
            filter(prot_id %in% peptides[["prot_id"]])

        # maximum values selected
        } else if (input[["type"]] == 2) {
          peptides_initial <- prot_number(data)
          peptides <- peptides_initial %>%
            filter(prot_number >= input[["n_pept2"]]) %>%
            inner_join(medians, by="prot_id") %>%
            arrange(desc(InROPE_median)) %>%
            head(n) %>%
            select(-c(prot_number, InROPE_median))
          seq_out <- seq_map <- data %>%
            filter(prot_id %in% peptides[["prot_id"]])
        }
      }
      
      # Join peptides with gene names and InROPE medians
      peptides_genes <- inner_join(peptides, genes, by="prot_id") %>%
        inner_join(medians, by="prot_id")
      
      
      # Outputs with results of filtering:
      output[["table"]] <- renderDataTable({
        my_DT(seq_out) %>% 
          formatRound(columns = c("muDiff", "effSz", "InROPE"),
                      digits = 4)
      })
      output[["proteins"]] <- renderDataTable({
        my_DT(peptides_genes)
      })
      # Display number of selected peptides
      output[["n_selected"]] <- renderText({
        paste0("Selected ", nrow(seq_out), " peptides.")
      })
      # Render heatmap
      output[["heatmap"]] <- renderPlot({
        ggplot(seq_map, aes(x = type, y = Sequence, fill = InROPE)) +
          geom_tile(color = "black") +
          scale_fill_gradient2(midpoint = 50) +
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
      my_DT(proteins)
    })
  })
}