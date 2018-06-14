ui <-  (navbarPage("Data Table", 

      ### density plot panel ###             
                   
tabPanel(("Density plot"),
    sidebarPanel(
      selectInput("variable", 
                  label="Select a variable:", 
                  choices = c("ben_WT", "ben_mal", "mal_WT"), 
                  selected = "ben_WT"), 
      br(),
      textOutput("text_density_plot"),
      plotOutput("density_plot")),
  
    mainPanel(
      dataTableOutput("table"))),

      ### sequences with min/max values panel ##

tabPanel(("Min/Max values"),
    sidebarPanel(   
      numericInput("x",
                   label="Select a number of sequences to show:",
                   value=0,
                   min=0,
                   max=100),      
      radioButtons("column",
                  label ="Variable:",
                  choices = list("ben_WT", "ben_mal", "mal_WT")),
      radioButtons("type",
                  label="Values:",
                  choices = list("Minimum"=1, "Maximum"=2)),
      actionButton("do", "Show")),
    
    mainPanel(dataTableOutput("filtered_sequences")
             ))
))