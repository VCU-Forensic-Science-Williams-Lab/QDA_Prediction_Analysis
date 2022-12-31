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
    numericInput("miR141", "miR141 (DNA Only)", value = 0),
    numericInput("miR412", "miR412 (DNA Only)", value = 0),
    numericInput("miR205", "miR205 (DNA Only)", value = 0),
    
    actionButton("submitRNA", "Submit RNA"),
    actionButton("submitDNA", "Submit DNA")
  ),
  
  dashboardBody(# Show plot
    plotOutput("results", width = "100%"))
  
)


server <- function(input, output, session) {
  # Proprietary QDA model
  finalRNA.qda.model <- readRDS('final.qda.model.rds')
  finalDNA.qda.model <- readRDS('finalDNA.qda.model.rds')
  
  # If RNA Submit button selected
  observeEvent(input$submitRNA, {
    data <- data.frame(input$miR200b,
                       input$miR320c,
                       input$miR10b,
                       input$miR891a)
    
    # Apply prediction model to data input
    pred2.qda <-
      predict(finalRNA.qda.model,
              newdata = data,
              method = "predictive")
    
    # Redefine axes input
    posterior = data.frame(t(pred2.qda$posterior))
    post_data = tibble::rownames_to_column(posterior, "class")
    
    # Graph
    output$results <- renderPlot(
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
    )
  })
  
  # If DNA Submit button selected
  observeEvent(input$submitDNA, {
    data <- data.frame(
      input$miR200b,
      input$miR320c,
      input$miR10b,
      input$miR891a,
      input$miR141,
      input$miR412,
      input$miR205
    )
    
    # Apply prediction model to data input
    pred2.qda <-
      predict(finalDNA.qda.model,
              newdata = data,
              method = "predictive")
    
    # Redefine axes input
    posterior = data.frame(t(pred2.qda$posterior))
    post_data = tibble::rownames_to_column(posterior, "class")
    
    # Graph
    output$results <- renderPlot(
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
    )
  })
  
}


# Run shiny app
shinyApp(ui, server)
