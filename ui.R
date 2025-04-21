library(shiny)
library(tidyverse)
library(ggplot2)
library(shinythemes)
options(scipen = 999)

fc <- read.csv('20250418_NL_BioID_FoldChanges.csv')

ui <- fluidPage(theme=shinytheme('spacelab'),
  titlePanel("Neuroligin BioID data"),

    sidebarLayout(
              sidebarPanel(selectInput("Gene", 
                             label = "Choose a gene",
                             choices = unique(fc$GeneSymbol),
                             selected='Nlgn2'),
                 
                     "Count Distributions",
                        plotOutput("Hist"),
                        plotOutput("HistN")
                     ),

    
    mainPanel(
      tabsetPanel(
        type='tabs',
        tabPanel("Enrichment",
        "Fold Change Relative to Experimental Control",
              h6('Enrichment compared to cytosolic Turbo BirA.
              P-values in each bar'),
              plotOutput("FC", height='500px', width='700px'),
              "Dotted line indicates 1.5 fold enrichment compared to control. \n Values above bar are p-values [bold and red indicates p-value < 0.05"),
        tabPanel("Reference",
                 "These data were generated using in vivo BioID in P21 mouse cortex. \n Please see publication at:")
    ))))

