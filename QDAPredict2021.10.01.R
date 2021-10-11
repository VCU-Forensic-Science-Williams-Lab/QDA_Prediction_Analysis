# ---
# title: "QDA Analysis with 10 fold cross-validation"
# output: html_document
# runtime: shiny
# author: m. creighton
# ---

# Libraries
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinyjs)
library(rsconnect)
library(ggplot2)
library(tidyverse)
library(forcats)
library(rio)
library(rpart)
library(rpart.plot)
library(stringr)
library(norm)
library(condMVNorm)
library(MASS)


# UI Setup
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Dashboard"),
  
  dashboardSidebar(
    width = 450,
    
    # Necessary Instuctions
    # * All miRs must have a numerical value in order to run the QDA Analysis.
    # * Differential expression or delta cycle threshold (ΔCq) values are calculated by 
    # subtracting the average Cq of let-7g and let-7i from the Cq value of the target miRNA.
    # * Formula: ΔCq = Cq(miRNA target) – Cq(avg let-7g & let-7i)
    
    h3("Instructions"),
    h4(
      "* All miRs must have a numerical value in order to run the QDA Analysis."
    ),
    h4(
      "* Differential expression or delta cycle threshold (ΔCq) values are calculated by
       subtracting the average Cq of let-7g and let-7i from the Cq value of the target miRNA."
    ),
    h4("* Formula: ΔCq = Cq(miRNA target) – Cq(avg let-7g & let-7i)"),
    
    # Sidebar with numeric inputs and action button
    numericInput("miR200b", "miR200b", value = 0),
    numericInput("miR320c", "miR320c", value = 0),
    numericInput("miR10b", "miR10b", value = 0),
    numericInput("miR891a", "miR891a", value = 0),

    actionButton("submit", "Submit")
  ),
  
  dashboardBody(
    
    # Show plot
    plotOutput("results", width = "100%"))
)

server <- function(input, output, session) {
  
  # Proprietary QDA model
  final.qda.model <- readRDS('final.qda.model.rds')
  
  result.plot <- eventReactive(input$submit, {
    
    # Pull data from user input
    data = data.frame(input$miR200b, input$miR320c, input$miR10b, input$miR891a)
    
    # Apply prediction model to data input
    pred2.qda <-
      predict(final.qda.model, newdata = data, method = "predictive")
    
    # Redefine axes input
    posterior = data.frame(t(pred2.qda$posterior))
    post_data = tibble::rownames_to_column(posterior, "class")
    
    # Graph
    ggplot(post_data, aes(x = class, y = t.pred2.qda.posterior.)) +
      geom_col() +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.3,
          hjust = .25,
          size = 13
        ),
        axis.text.y = element_text(vjust = 0.5)
      ) +
      xlab("Body Fluid") +
      ylab("Probability")
    
  })
  
  output$results <- renderPlot(result.plot())
  
  # options = list(height = 600)
  
}

# Run shiny app
shinyApp(ui, server)
