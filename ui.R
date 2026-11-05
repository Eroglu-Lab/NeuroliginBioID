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
              "Dotted line indicates 1.5 fold enrichment compared to control. \n Values above bar are p-values [bold and red indicates p-value < 0.05]"),
                tabPanel("Reference",
                 "These data were generated using in vivo BioID in P21 mouse cortex.",
                         br(),
                         br(),
                         "Please see publication at:",
                         br(),
                         br(),
                  "Sakers K, Ramirez JJ, Elazar N, Nagendren L, Soderblom E, Eroglu C. Neuroligin-2 is ubiquitinated by Nedd4l to control developmental astrocyte morphogenesis.
                  J Cell Biol. 2026 Jul 6;225(7):e202512111. doi: 10.1083/jcb.202512111. Epub 2026 Jun 1. PMID: 42223381; PMCID: PMC13276441.",
                         br(),
                         br(),
                 tags$p(
                   tags$a("Click here",
                          href = "https://rupress.org/jcb/article/225/7/e202512111/282648/Neuroligin-2-is-ubiquitinated-by-Nedd4l-to-control",
                          target = "_blank")
                 )
        )   # closes tabPanel("Reference")
      )     # closes tabsetPanel
    )       # closes mainPanel
  )         # closes sidebarLayout
)           # closes fluidPage

