server <- function(input, output) {
  
  counts <- read.csv('20240910_normalizedCounts.csv')
  fc <- read.csv('20240910_NL_BioID_FoldChanges.csv')
  pct <- read.csv('20240910_percentiles.csv')
  pct$conditionLong <- paste(pct$cell, pct$BioID, sep=' ')
  
##############
# Prep  Data #
##############
  
  pct2 <- pct %>% 
    dplyr::mutate(conditionLong = factor(conditionLong, 
                                      levels = c("Astrocyte Turbo1","Astrocyte NL1", "Astrocyte Turbo2","Astrocyte NL2",
                                                 "Astrocyte Turbo3","Astrocyte NL3", "Neuron Turbo", "Neuron NL2" )))
  
  counts <- counts %>% 
    dplyr::mutate(condition = factor(condition, 
                                     levels = c("Turbo1_A", "NL1_A","Turbo2_A","NL2_A", "Turbo3_A", "NL3_A",
                                                "Turbo_N", "NL2_N" )))
  
  fc$pval_rnd <- formatC(fc$pval, format='e', digits=2)
  
  fc <- fc %>% 
    dplyr::mutate(BioID = factor(BioID, 
                                     levels = c("NL1", "NL2","NL3")),
                  cell = factor(cell, 
                                levels = c('Astrocyte', 'Neuron')))
  
  
  

  counts_mean <- aggregate(counts$counts, by=list(counts$gs_long, counts$condition), mean)
  counts_mean_A <- counts_mean[grep("_A", counts_mean$Group.2), ]
  counts_mean_N <- counts_mean[grep("_N", counts_mean$Group.2), ]

  tab_pct <- reactive({pct2 %>%
      filter(gs_long2 == input$Gene)})
  
  tab_counts <- reactive({ counts %>%
      filter(gs_long == input$Gene)})
  
  tab_meanC_A <- reactive({ counts_mean_A %>% 
      filter(Group.1 == input$Gene)})
  
  tab_meanC_N <- reactive({ counts_mean_N %>% 
      filter(Group.1 == input$Gene)})
  
  observeEvent(input$Gene, {
    output$Hist <- renderPlot({
    
    hist(log(counts_mean_A$x), breaks='FD', border='lightgray', ylim=c(0,800), xlim=c(0,25),
         main="Astrocyte Count Distribution", xlab='log(Normalized Counts)')
    abline(v=log(tab_meanC_A()[,3]), 
           col=c("palegreen3", "palegreen4","hotpink3", "hotpink4", "darkorchid1", "darkorchid4"), lwd=4)
    legend(x=0, y=700, legend = c('Turbo1', 'NL1', 'Turbo2', 'NL2', 'Turbo3', 'NL3'),
           col=c("palegreen3", "palegreen4","hotpink3", "hotpink4", "darkorchid1", "darkorchid4"),
           lwd=4, bty='n')
  })

})
  
  observeEvent(input$Gene, {
    output$HistN <- renderPlot({
      
      hist(log(counts_mean_N$x), breaks='FD', border='lightgray', ylim=c(0,400), xlim=c(0,25), main="Neuron Count Distribution", xlab='log(Normalized Counts)')
      abline(v=log(tab_meanC_N()[,3]), 
             col=c("firebrick1", "firebrick4"), lwd=4)
      legend(x=0, y=350, legend = c('Turbo', 'NL2'),
             col=c("firebrick1", "firebrick4"),
             lwd=4, bty='n')
    })
    
  })

  
tab_fc <- reactive({
    fc %>%
      filter(GeneSymbol == input$Gene)})
  
observeEvent(input$Gene, {
  output$FC <- renderPlot({ 
    
    ggplot(tab_fc(), aes(x=BioID, y=FC, fill=cell)) + 
      geom_bar(stat='identity', position = 'dodge') + scale_x_discrete(drop=F) + 
      geom_hline(linetype=2, col='black', yintercept=0, lwd=1) + 
      geom_text(aes(label=pval_rnd), position=position_dodge(width=0.9), vjust=-0.5, alpha=1, size=5, color= 'black') +
      theme_bw(base_size=20)+ theme(panel.border=element_blank(), axis.line=element_line(), 
                                    axis.text.x=element_text(size=15, color='black'), 
                                    axis.text.y=element_text(size=15, color='black'),
                                    axis.title.x = element_text(size=14),
                                    axis.title.y = element_text(size=14)) +
      ylab('Fold change \n relative to control') + xlab('Bait Protein') +
      scale_fill_manual(values = c("palegreen4","hotpink4"))
    
  })
})

}
