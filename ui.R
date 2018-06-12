ui <- fluidPage(
  titlePanel("Data table"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("variable", 
                  label="Select a variable:", 
                  choices = c("ben_WT", "ben_mal", "mal_WT"), 
                  selected = "ben_WT"), 
      br(),
      textOutput("text_density_plot"),
      
      plotOutput("density_plot")),
    mainPanel(
              dataTableOutput("table")
              )
             )
)