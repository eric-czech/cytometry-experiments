library(flowWorkspace)
library(flowCore)
library(ggcyto)

ws <- openWorkspace(
  '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/human-pbmc-cytof/data/attachments/081116-Mike-HIMC_controls.wsp')
gs <- parseWorkspace(
  ws, path='/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/human-pbmc-cytof/data/wsp',
  #subset=c('081216-Mike-HIMC ctrls-001_01_1.fcs'),
  sampNloc='sampleNode'
)
getData(gs)[[1]]
getNodes(gs)
autoplot(gs[[1]], '/Cells/Intact cells/Intact singlets')