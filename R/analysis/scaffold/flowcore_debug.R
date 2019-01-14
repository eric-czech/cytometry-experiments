library(flowCore)
library(ggcyto)
library(openCyto)

set.seed(1)
n <- 10000
get_fr <- function(m) flowFrame(as.matrix(data.frame(f1=rnorm(n, mean=m), f2=rnorm(n, mean=m))))
fs <- as(list(sample1=get_fr(0), sample2=get_fr(3.5)), 'flowSet')
gs <- GatingSet(fs)

add_pop(gs, alias='A,B,C,D', pop='*', parent='root', dims='f1,f2', gating_method = 'quadGate.seq', gating_args = "gFunc='mindensity2'")

#recompute(gs)
#add_pop(gs, alias='A,B,C,D', pop='*', parent='root', dims='f1,f2', gating_method = 'quadGate.seq', gating_args = "gFunc='tailgate'")
autoplot(gs, c('/A', '/B', '/C', '/D')) + xlim(-5, 8) + ylim(-5, 8)

getGate(gs, '/A')
