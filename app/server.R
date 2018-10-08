library(shiny)
library(openCyto)
library(flowCore)
library(ggcyto)

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

  output$distPlot <- renderPlot({
    d <- read.FCS('~/cytoui/data/tcell-cd4-cd8.fcs')
    chnl <- c("ciCD4", "ciCD8")
    trans <- transformList(chnl, logTransform())
    d2 <- transform(d, trans)
    g1 <- openCyto:::.mindensity(d2, channels = chnl[1])
    g2 <- openCyto:::.mindensity(d2, channels = chnl[2])
    autoplot(d2, chnl[1], chnl[2]) + geom_gate(g1) + geom_gate(g2)
  })

})