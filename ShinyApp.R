#################################################-
## USER INTERFACE
#################################################-
## Preliminaries ----
#################################################-
library(shiny)
library(shinyMatrix)
library(scales)
library(tidyverse)
library(medicaldata)
library(mlbench)
library(plotly)
library(ggplot2)
library(ggpubr)
library(rmarkdown)
library(knitr)
library(pander)
library(ggforce)

#################################################-
## Define UI ----
#################################################-
#setwd("/Users/lisabuche/Documents/Projects/Spaak/shiny")
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel(
    h1("How to compute niche and fitness differences with intraspecific interactions", h2("An inclusive approach to coexistence theory"))
  ),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h4("Scenario of two species interacting."), 
      h4("With intrinsic growth rate of species j, constant and equal to 1."),
        h4("Only change the interaction coefficients specific to species i."),
      hr(),
      h4("Scenarios:"),
      actionButton("alphamat_to_NFD",
                   "Select interactions value"),
      actionButton("NFD_to_alphamat_with_intra",
                   "Select ND and FD, with α(i,i) being facilitative"),
      actionButton("NFD_to_alphamat_no_intra",
                   "Select ND and FD, with α(i,i) being competitive") ,
      hr(),
      h4("Time of time series:"),
      sliderInput("time",
                  label =withMathJax(""),
                  min = 5,
                  max = 50,
                  value = 5,
                  step=1),
      hr(),
      h4("Species parameters:"),
      sliderInput("mu_i",
                  label =withMathJax("Intrinsic growth rate \\(\\mu_{i}\\)"),
                  min = -1,
                  max = 1,
                  value = 1,
                  step=2),
      hr(),
      h4("Species interactions:"),
      shinyMatrix::matrixInput(inputId = "alphamat",
                               value = matrix(c(1, 0.4, 
                                                0.2,1),
                                              nrow = 2,
                                              dimnames = list(c("α(i,_)", "α(j,_)"), c("α(_,i)", "α(_,j)")),
                                              byrow = TRUE),
                               class = "numeric",
                               rows = list(names = TRUE,
                                           editableNames = FALSE),
                               cols = list(names = TRUE,
                                           editableNames = FALSE)),
      hr(),
      h4("Niche and fitness difference of species i:"),
      shinyMatrix::matrixInput(inputId = "NFD",
                               value = matrix(c(0.7171573,
                                                0.293),
                                              nrow = 1,
                                              dimnames = list(c(""), c("ND", "FD")),
                                              byrow = TRUE),
                               class = "numeric",
                               rows = list(names = FALSE,
                                           editableNames = FALSE),
                               cols = list(names = TRUE,
                                           editableNames = FALSE))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(fluidRow(
      splitLayout(cellWidths = c("40%", "55%"), 
                  tags$b("a. Species interactions"),
                  tags$b("b. Abundance over time when species i invades species j \nand for species i in monoculturure")),
    ),
    fluidRow(
      splitLayout(cellWidths = c("40%", "55%"), 
                  plotOutput("interaction.plot"), 
                  plotOutput("abundance.plot"))
      ),
    fluidRow(
      splitLayout(cellWidths = c("55%", "40%"), 
                tags$b('c. Growth rate landscape of species i'),
                tags$b("d. Graph of niche and fitness difference MAP"))
      ),
    fluidRow(
        splitLayout(cellWidths = c("55%", "40%"), 
                    plotlyOutput("ThreeDplot"), 
                    plotOutput("TwoDplot"))
        )
      )
    )
  )


#################################################-
# Shiny server logic ----
#################################################-

#source("shiny/Shinytoolbox.R")
# This section takes inputs from the ui and displays the desired
# information using the code from Saavedra et al. above

server <- function(input, output, session) {
  
  observeEvent(input$alphamat_to_NFD, {
    mu <- c(input$mu_i,1)
    A <- input$alphamat
    pars <- list(mu = mu, A = A)
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    pars$sign_A11 <- sign(pars$A[1,1])
    pars <- compute_NFD(pars)
 
    updateMatrixInput(session = session,
                      inputId = "alphamat",
                      value = matrix(c(pars$A[1,1], pars$A[1,2], 
                                       pars$A[2,1],pars$A[2,2]),
                                     nrow = 2,
                                     dimnames = list(c("α(i,_)", "α(j,_)"), c("α(_,i)", "α(_,j)")),
                                     byrow = TRUE))
    updateMatrixInput(session = session,
                      inputId = "NFD",
                      value = matrix(c(pars$ND, pars$FD),
                                     nrow = 1,
                                     dimnames = list(c(""), c("ND", "FD")),
                                     byrow = TRUE))
  })
  
  observeEvent(input$NFD_to_alphamat_with_intra,{
    mu <- c(input$mu_i,1)
    ND <- input$NFD[1]
    FD <- input$NFD[2]
    mu <- c(input$mu_i,1)
    A <- input$alphamat
    pars <- list(mu = mu, A = A)
    pars$sign_A11 <- -1
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    pars <- compute_LV_from_NFD(ND,FD,pars)
    updateMatrixInput(session = session,
                      inputId = "alphamat",
                      value = matrix(c(pars$A[1,1], pars$A[1,2], 
                                       pars$A[2,1],pars$A[2,2]),
                                     nrow = 2,
                                     dimnames = list(c("α(i,_)", "α(j,_)"), c("α(_,i)", "α(_,j)")),
                                     byrow = TRUE))
  })
  
  observeEvent(input$NFD_to_alphamat_no_intra,{
    mu <- c(input$mu_i,1)
    ND <- input$NFD[1]
    FD <- input$NFD[2]
    mu <- c(input$mu_i,1)
    A <- input$alphamat
    pars <- list(mu = mu, A = A)
    pars$sign_A11 <- 1
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    pars <- compute_LV_from_NFD(ND,FD,pars)
    updateMatrixInput(session = session,
                      inputId = "alphamat",
                      value = matrix(c(pars$A[1,1], pars$A[1,2], 
                                       pars$A[2,1],pars$A[2,2]),
                                     nrow = 2,
                                     dimnames = list(c("α(i,_)", "α(j,_)"), c("α(_,i)", "α(_,j)")),
                                     byrow = TRUE))
  })
  
  
  output$interaction.plot <- renderPlot({
    mu <- c(input$mu_i,1)
    A <- input$alphamat
    A[2,1] <- 0.2
    A[2,2] <- 1
    pars <- list(mu = mu, A = A)
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    pars$sign_A11 <- sign(pars$A[1,1])
    
    mk_graph_2sp(A,mu)
    
  })
  output$abundance.plot <- renderPlot({
    mu <- c(input$mu_i,1)
    A <- input$alphamat
    pars <- list(mu = mu, A = A)
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    pars$sign_A11 <- sign(pars$A[1,1])
    
    pars <- compute_NFD(pars)
    
    graph_density_lv(pars,input$time)
    
  })
  
  output$ThreeDplot <- renderPlotly({
    mu <- c(input$mu_i,1)
    A <- input$alphamat
    pars <- list(mu = mu, A = A)
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    pars$sign_A11 <- sign(pars$A[1,1])
    
    pars <- compute_NFD(pars)
 
    fig <- graph_NFD_plan(pars)
    
  })
  
  output$TwoDplot <- renderPlot({
    ND <- input$NFD[1]
    FD <- input$NFD[2]
    A <- input$alphamat
    pars <- list(ND = ND,  FD =  FD,A = A)
    pars$A[2,1] <- 0.2
    pars$A[2,2] <- 1
    
    graph_NFD_2D(pars)
    
  })
  

  
}

#################################################-
## Run the application  ----
#################################################-
shinyApp(ui = ui, server = server)
