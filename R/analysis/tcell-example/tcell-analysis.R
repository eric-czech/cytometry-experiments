library(openCyto)
library(flowCore)
library(ggcyto)
library(data.table)

# Example: http://opencyto.org/gating_using_opencyto.html
# Vignette: https://bioconductor.org/packages/release/bioc/vignettes/openCyto/inst/doc/openCytoVignette.html
# Workshop (incremental gating): https://www.bioconductor.org/help/course-materials/2015/BioC2015/OpenCytoWorkshop.html
# Manual Gating: http://rpubs.com/wjiang2/395270

d <- read.FCS('~/analysis/data/tcell-cd4-cd8.fcs')

chnl <- c("ciCD4", "ciCD8")
d <- transform(d, transformList(chnl, logTransform()))

ggcyto(d, aes(x='celldiameter')) + geom_density()

ggcyto(d, aes(x='ciCD4', y='ciCD8')) + geom_hex() 

ggcyto(d, aes(x="ciCD4", y="ciCD8")) +
  geom_hex() +
  geom_gate(openCyto:::quadGate.tmix(d, channels=c("ciCD4", "ciCD8"), K=3)) +
  geom_stats()

ggcyto(d, aes(x="ciCD4")) + geom_histogram(bins=30) +
  geom_gate(openCyto:::.mindensity(d, channels="ciCD4")) +
  geom_stats()

# Building a template
gs <- GatingSet(flowSet(d))
template = add_pop(
  gs, alias="nonLargeCells", pop="-", parent="root", dims="celldiameter", 
  gating_method="tailgate", gating_args="tol=0.02,side='right'"
)
plotGate(gs, 'nonLargeCells', default.y='ciDAPI')
getPopStats(gs)

template = rbind(template, add_pop(
  gs, alias="cells", pop="+", parent="nonLargeCells", dims="celldiameter", 
  gating_method="tailgate", gating_args="tol=0.02,side='left'"
))

plotGate(gs, 'cells', default.y='ciDAPI')
getPopStats(gs)

template = rbind(template, add_pop(
  gs, alias="*", pop="+/-+/-", parent="cells", dims="ciCD4,ciCD8", 
  gating_method="mindensity", gating_args=""
))

plotGate(gs, default.y='ciDAPI')

plotGate(gs[[1]], default.y='ciDAPI')




dgs <- d[,c('ciCD4', 'ciCD8', 'ciDAPI', 'celldiameter', 'cellsolidity')]
gt <- gatingTemplate('~/analysis/R/scripts/tcell-example/gating.csv')
transformer <- transformerList(c('ciCD4', 'ciCD8'), logicle_trans())
gs <- transform(GatingSet(flowSet(dgs)), transformer)



gating(gt, gs)
plotGate(gs[[1]], default.y="ciDAPI", xbin=0)
outFile <- tempfile(fileext = ".wsp")
outFile <- '~/analysis/test.wsp'
GatingSet2flowJo(gs, outFile)
