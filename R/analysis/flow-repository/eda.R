# install.packages("BiocManager")
# BiocManager::install("FlowRepositoryR", update=F)
# BiocManager::install("flowDensity", update=F)

library(FlowRepositoryR)
library(flowCore)
library(openCyto)
library(ggcyto)
library(dplyr)
library(purrr)


# This is crucial as quadGate.tmix uses filter symbol from parent env (a 
# call that is wrapped in a silent try/catch leading to an obscure error)
filter <- flowCore::filter
data_dir <- '/tmp/fr'
dsids <- flowRep.ls()

dsids <- flowRep.search("OMIP-021")

ds <- flowRep.get(dsids[[1]])
summary(ds)

ds <- download(ds, dirpath=data_dir)
summary(ds)

ds_files <- list.files(path=data_dir, pattern="Donor*", full.names=T)

### Programmatic Definitions

# path <- localpath(fcs.files(ds)[[1]])
path <- ds_files[1]
d <- read.FCS(path)

cm <- with(parameters(d)@data, setNames(name, coalesce(desc, name)))

#d <- transform(d, transformList(c('525/50Violet-A'), logicleTransform()))
d <- transform(d, transformList(cm, logicleTransform()))

# d_comp <- compensate(d_raw, keyword(d_raw)$SPILL)

# Supplementary info doc: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fcyto.a.22475&file=cytoa22475-sup-0005-suppinfo01.doc
# States in 4.2 "Compensation Description":
# > Manual compensation using BD DIVA software and single stained cells
# > Compensation was performed immediately prior to cell acquisition
# > Compensation controls uploaded to flowrepository.org (see item 4.1 above)

# OMIP-021
# Innate-like T-cell Panel
# Donor1.fcs
# https://sci-hub.tw/10.1002/cyto.a.22475
# Figure 1


apply_gates <- function(d, gates) purrr::reduce(gates, flowCore::Subset, .init=d)

# ggcyto(d, aes(x='FSC-A', y='L/D Aqua - Viability')) + geom_hex(bins=100) 

# g <- gate_tail(d, cm[['L/D Aqua - Viability']], num_peaks=3, ref_peak=3)
# g1 <- polygonGate(.gate=tribble(
#   ~"FSC-A", ~"525/50Violet-A",
#   0, 0,
#   4.6, 0,
#   4.6, 2.2,
#   4.0, 2.0,
#   3.5, 1.4,
#   0, 1
# ))
#g1 <- flowClust.2d(d, cm[['FSC-A']], cm[['L/D Aqua - Viability']], target=c(3.8, 2.8), K=6)
gates <- list()
gates[[1]] <- openCyto:::.singletGate(d, channels=c(cm[['FSC-A']], cm[['L/D Aqua - Viability']]))

# ggcyto(d, aes(x='L/D Aqua - Viability')) + geom_histogram(bins=100) 
# add_pop(
#   gs, alias="LiveCells", pop="LiveCells", parent="root", dims="FSC-A,525/50Violet-A", 
#   gating_method="mindensity", gating_args=""
# )
ggcyto(d, aes(x='FSC-A', y='L/D Aqua - Viability')) + geom_hex(bins=100) + geom_gate(gates[[1]]) + geom_stats() +
  xlim(2.5, 4.5) + ylim(0, 3.5)

gates[[2]] <- gate_flowClust_2d(apply_gates(d, gates[1]), 'FSC-A', 'SSC-A', K=3, quantile=.95)
gates[[2]] <- flowClust.2d(apply_gates(d, gates[1]), 'FSC-A', 'SSC-A', K=3, quantile=.95)
ggcyto(apply_gates(d, gates[1]), aes(x='FSC-A', y='SSC-A')) + 
  geom_hex(bins=100) + geom_gate(gates[[2]]) + geom_stats() +
  xlim(3.5, 4.5) + ylim(3.0, 4.5)

# gates[[3]] <- polygonGate(.gate=tribble(
#   ~"FSC-A", ~"FSC-H",
#   2.8, 0,
#   4.5, 4.2,
#   4.5, 4.5,
#   0, 1
# ))
gates[[3]] <- openCyto:::.singletGate(apply_gates(d, gates[1:2]), channels=c(cm[['FSC-A']], cm[['FSC-H']]), maxit=1000, wider_gate=T, prediction_level=.9999999)
ggcyto(apply_gates(d, gates[1:2]), aes(x='FSC-A', y='FSC-H')) + 
  geom_hex(bins=100) + geom_gate(gates[[3]]) + geom_stats()


ggcyto(d, aes(x='BV785 - CD3', y='BV605 - CD161')) + geom_hex(bins=100) #+ geom_gate(g3) + geom_stats()


# ggcyto(d, aes(x='BV785 - CD3', y='FITC - gdTCR')) + geom_hex(bins=100) + scale_y_logicle() + scale_x_logicle()

##### Template Version

.polyGate <- function(fr, pp_res, channels, filterId="polygate", ...){ 
  args <- list(...)
  g <- data.frame(x=args$x, y=args$y)
  colnames(g) <- channels
  flowCore::polygonGate(.gate=g, filterId=filterId)
}
registerPlugins(fun=.polyGate, methodName='polyGate',dep=NA)

path <- ds_files[1]
d <- read.FCS(path)
colnames(d) <- make.names(colnames(d))
colnames(d@description$SPILL) <- make.names(colnames(d@description$SPILL))
cm <- with(parameters(d)@data, setNames(name, coalesce(desc, name)))
transformer <- transformerList(colnames(d@description$SPILL), logicle_trans())
gs <- transform(GatingSet(flowSet(d)), transformer)
gt <- gatingTemplate('~/analysis/R/scripts/flow-repository/gating.csv')
gating(gt, gs)

plot_density2d <- function(fr, x, y, point_alpha=.3, point_size=.2){
  ggcyto(fr, aes_string(x=x, y=y)) + 
    geom_point(alpha=point_alpha, size=point_size, color='blue') +
    stat_density2d(aes(alpha=..level.., fill=..level..), geom="polygon") +
    scale_fill_gradientn(colours = rev(rainbow(20))) +
    theme_bw() + theme(
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_blank()
    )
}

autoplot(gs[[1]], strip.text = "gate", bins=100, merge=F, axis_inverse_trans=F) +
  theme_bw() + theme(
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_blank()
  )

autoplot(gs, gate='Live', bins=100)
autoplot(gs, gate='Lymphocytes', bins=100)
autoplot(gs, gate='Singlets', bins=100)


# xvar <- 'BV785 - CD3'
# yvar <- 'BV421 - TETaGC'
# xvar <- 'BV785 - CD3'
# yvar <- 'FITC - gdTCR'
xvar <- 'PE - TCR Va7'
yvar <- 'BV605 - CD161'

# g <- openCyto::flowClust.2d(
#   getData(gs[[1]], 'abTCells'), 
#   cm[[xvar]], cm[[yvar]], 
#   K=7, target=c(3.5, 3.5)
# )

g <- openCyto:::.quadGate.seq(
  getData(gs[[1]], 'abTCells'), 
  channels=c("X582.15Yellow.A", "X540.30Violet.A"),
  gFun='tailgate', tol = 0.01
)
ggcyto(getData(gs[[1]], 'abTCells'), aes_string(x=paste('`', xvar, '`', sep=""), y=paste('`', yvar, '`', sep=""))) + 
  geom_hex(bins=100) + geom_gate(g) + geom_stats()

ggcyto(getData(gs[[1]], 'abTCells'), aes_string(x=cm[[xvar]], y=cm[[yvar]])) + 
  geom_hex(bins=100) + geom_gate(g) + geom_stats()

ggcyto(getData(gs[[1]], 'Live'), aes(x='BV785 - CD3', y='BV421 - TETaGC')) + 
  geom_hex(bins=100)


plot(gs)
# Remove a gate: Rm('maiTCells', gs)
.polyGate <- function(fr, pp_res, channels, filterId="polygate", ...){ 
  args <- list(...)
  g <- data.frame(x=args$x, y=args$y)
  colnames(g) <- channels
  print(g)
  flowCore::polygonGate(.gate=g, filterId=filterId)
}
registerPlugins(fun=.polyGate, methodName='polyGate',dep=NA)
add_pop(
  gs, alias="maiTCells", pop="+", parent="abTCells", dims="X582.15Yellow.A,X540.30Violet.A", 
  gating_method="polyGate", gating_args="x=c(2.7,5,5,2.7),y=c(2.5,2.5,5,5)"
)
ggcyto(getData(gs[[1]], 'abTCells'), aes(x='X582.15Yellow.A', y='X540.30Violet.A')) + 
  geom_hex(bins=100) + geom_gate(getGate(gs[[1]], 'maiTCells')) + geom_stats()
plot_density2d(getData(gs[[1]], 'maiTCells'), 'X582.15Yellow.A', 'X540.30Violet.A')

### Manual Gating
library(flowDensity)
flowDensity::plotDens(getData(gs[[1]], 'abTCells'), c('X582.15Yellow.A', 'X540.30Violet.A'))
coords <- locator(n = 6, type = "p", lwd = 2, pch = 16, col = "red")
gate <- data.frame(coords)
colnames(gate) <- c('X582.15Yellow.A', 'X540.30Violet.A')
g <- polygonGate(.gate=gate)

plot_density2d(getData(gs[[1]], 'abTCells'), 'X582.15Yellow.A', 'X540.30Violet.A') +
  geom_gate(g) + geom_stats()
# No way to add these to gs with openCyto 1.18.0

