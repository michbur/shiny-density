ui <-  (fluidPage(
  
  sidebarPanel(
    radioButtons("method",
                 label = "Method of filtering:",
                 choices = c("Density plot", "Minimum/maximum values"), 
                 selected = character(0)),
    selectInput("variable", 
                label = "Select a variable:", 
                choices = c("ben_WT", "ben_mal", "mal_WT")), 
    br(),
    textOutput("text_density_plot"),
    plotOutput("density_plot", 
               click = "plot_click"),
    textOutput("test"),
    br(),
    
    conditionalPanel(
      condition = "input.method == 'Density plot'",
      radioButtons("thresh_button", 
                   "Select values lower or higher than threshold?",
                   choices = list("Lower" = 1, "Higher" = 2))),
    
    conditionalPanel(
      condition = "input.method == 'Minimum/maximum values'",  
      radioButtons("type",
                   "Select the minimum or maximum values?",
                   choices = list("Minimum" = 1, "Maximum" = 2)),
      numericInput("x",
                   label = "Number of sequences to show:",
                   value = 0,
                   min = 0,
                   max = 100)),
    actionButton('filter',
                 "Filter"),
    actionButton("reset",
                 "Reset")),
  
  mainPanel(
    dataTableOutput("table"))
))