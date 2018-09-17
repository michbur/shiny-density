library(shiny)
library(DT)

 if(Sys.info()[["nodename"]] == "amyloid")
   sample1 <- read.csv("/home/michal/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")
 if(Sys.info()[["nodename"]] == "LENOVO")
   sample1 <- read.csv("C:/Users/Kaede/Dropbox/PepArray_results/2018-06-07/full_best_res.csv")

ui <-  (fluidPage(title = "PAAP",
  sidebarPanel(
    radioButtons("method",
                 label = "Method of filtering:",
                 choices = c("Density plot", "Minimum/maximum values")),
                 #selected = character(0)),
    # br(),
    conditionalPanel(
      condition = "input.method == 'Density plot'",
      selectInput("var", 
                  label = "Select a type:", 
                  choices = levels(sample1[["type"]])),
      # br(),
      textOutput("text_density_plot"),
      includeMarkdown("density_plot.md"),
      plotOutput("density_plot", 
                 click = "plot_click"),
      textOutput("test"),
      textOutput("n_selected"),
      radioButtons("thresh_button", 
                   "Select values lower or higher than threshold?",
                   choices = list("Lower" = 1, "Higher" = 2)),
      radioButtons("level1",
                   "Filter proteins or peptides?",
                   choices = list("Peptides", "Proteins")),
    conditionalPanel(
      condition = ("input.level1 == 'Proteins'"),
      numericInput("n_pept1",
                   label = "Minimal number of relevant peptides in a protein:",
                   value = 1,
                   min = 1,
                   max = 20))),
    
    conditionalPanel(
      condition = "input.method == 'Minimum/maximum values'",  
      radioButtons("type",
                   "Select the minimum or maximum values?",
                   choices = list("Minimum" = 1, "Maximum" = 2)),
      numericInput("x",
                   label = "Number of sequences/proteins to show:",
                   value = 0,
                   min = 0,
                   max = 100),
      radioButtons("level2",
                   "Filter proteins or peptides?",
                   choices = list("Peptides", "Proteins")),
      conditionalPanel(
        condition = "input.level2 == 'Proteins'",
        numericInput("n_pept2",
                     label = "Minimal number of peptides in a protein:",
                     value = 1,
                     min = 1,
                     max = 20))),
    # conditionalPanel(
    #   condition = "input.method == 'Protein level'",
      # numericInput("n_pept",
      #              label = "Select proteins with at least",
      #              value = 1,
      #              min = 1,
      #              max = 20),
    #   radioButtons("pept_types",
    #                      label = "peptides belonging to:",
    #                      choices = levels(sample1[["type"]])),
    #   checkboxInput("plot_filter",
    #                 "Apply density plot threshold",
    #                 FALSE)),
    actionButton('filter',
                 "Filter"),
    actionButton("reset",
                 "Reset"),
    includeMarkdown("filter_reset.md")),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Data",
               includeMarkdown("filtering.md"),
               dataTableOutput("table")),
      tabPanel("Proteins", 
               includeMarkdown("proteins.md"),
               dataTableOutput("proteins")),
      tabPanel("Heatmap",
               includeMarkdown("heatmap.md"),
               plotOutput("heatmap"))      
    )
  )
))