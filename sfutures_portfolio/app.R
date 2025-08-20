library(shiny)
library(shinydashboard)
library(ggplot2)

# Sample data for plots
data <- mtcars

ui <- dashboardPage(
  dashboardHeader(title = "Simple Dashboard"),
  dashboardSidebar(),
  dashboardBody(
    fluidRow(
      box(plotOutput("plot1"), width = 6),
      box(plotOutput("plot2"), width = 6)
    ),
    fluidRow(
      box(actionButton("btn", "Click Me"), width = 12)
    )
  )
)

server <- function(input, output, session) {
  
  output$plot1 <- renderPlot({
    ggplot(data, aes(x = wt, y = mpg)) +
      geom_point() +
      ggtitle("Plot 1: Weight vs MPG")
  })
  
  output$plot2 <- renderPlot({
    ggplot(data, aes(x = factor(cyl), fill = factor(cyl))) +
      geom_bar() +
      ggtitle("Plot 2: Cylinder Count")
  })
  
  observeEvent(input$btn, {
    showModal(modalDialog(
      title = "Button Clicked",
      "You clicked the button!"
    ))
  })
}

shinyApp(ui, server)
