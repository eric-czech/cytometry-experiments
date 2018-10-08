library(openCyto)
library(flowCore)
library(ggcyto)
library(data.table)

# Example: http://opencyto.org/gating_using_opencyto.html

d <- read.FCS('~/analysis/data/tcell-cd4-cd8.fcs')
chnl <- c("ciCD4", "ciCD8", "ciDAPI")

#trans <- transformList(chnl, logTransform())
trans <- transformList(chnl, logicleTransform())
d2 <- transform(d, trans)

g1 <- openCyto:::.mindensity(d2, channels = chnl[1])
g2 <- openCyto:::.mindensity(d2, channels = chnl[2])
autoplot(d2, chnl[1], chnl[2]) + geom_gate(g1) + geom_gate(g2)
getPopStats(d)


dgs <- d[,c('ciCD4', 'ciCD8', 'ciDAPI', 'celldiameter', 'cellsolidity')]
gt <- gatingTemplate('~/analysis/R/scripts/tcell-example/gating.csv')
transformer <- transformerList(c('ciCD4', 'ciCD8'), logicle_trans())
gs <- transform(GatingSet(flowSet(dgs)), transformer)

gating(gt, gs)
plotGate(gs[[1]], default.y="ciDAPI", xbin=0)

outFile <- tempfile(fileext = ".wsp")
outFile <- '~/analysis/test.wsp'
GatingSet2flowJo(gs, outFile)
