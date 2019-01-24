library(tidyverse)
library(flowWorkspace)
library(flowCore)
library(ggcyto)

ws <- openWorkspace(
  '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/human-pbmc-cytof/data/attachments/081116-Mike-HIMC_controls.wsp')
gs <- parseWorkspace(
  ws, path='/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/human-pbmc-cytof/data/fcs',
  #subset=c('081216-Mike-HIMC ctrls-001_01_1.fcs'),
  sampNloc='sampleNode',
  #name='PBMC control cells-MATLAB normalized'
  name='PBMC control cells-Fluidigm-ver2 normalized'
)
# getData(gs)[[1]]
getNodes(gs)

cols_linear <- c(
  'Time', 
  'Event_length', 
  # Environmental contaminants
  'Sn120Di', 'I127Di', 'Xe131Di', 'Ba138Di', 'Pb208Di',
  'In115Di', # 115In_Dead
  'Ce140Di', # 140Ce_Bead 
  'Ir191Di', # 191Ir_DNA1
  'Ir193Di', # 191Ir_DNA2
  'Center', 
  'Offset',
  'Width',
  'Residual'
)
cols_marker <- colnames(gs) %>% setdiff(cols_linear)

# Transform on gatingset is "fasinh" but inversion does little to move min-max range
# - min/max in parameters df do not match actual data somehow so they are rescaled here
#trans <- transformerList(cols_trans, flowCore_asinht_trans(a=0, b=1/5, c=0))
# gs_inv_trans <- getTransformations(gs[[1]], inverse=T, only.function=F)
# gs_inv_trans <- transformerList(names(gs_inv_trans), gs_inv_trans)
# gsi <- transform(gs, gs_inv_trans)

# fs <- getData(gs)
# getSampleGroups(ws)

# getData(gs, '/Cells/Intact cells/Intact singlets/Live intact singlets/Nonbasophils')

getNodes(gs)
getDescendants(gs[[1]], '/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+')
autoplot(gs, gate= "/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+",y='Ir191Di', bins=64, axis_inverse_trans=F, merge=F) + 
  ggcyto_par_set(limits='data')


# getGate(gs[[1]], '/Cells/Intact cells/Intact singlets/Live intact singlets/Nonbasophils')
# plotGate(gs, "/Cells/Intact cells", xlim='data', ylim='data')
# autoplot(gs, gate= "/Cells/Intact cells/Intact singlets/Live intact singlets/CD14+ CD33+/CD14+ CD16+", bins=64, axis_inverse_trans=F, merge=F) + ggcyto_par_set(limits='data')
# ggcyto(gs, aes(x='140Ce_Bead', y='191Ir_DNA1')) + geom_gate('/Beads') + ggcyto_par_set(limits='data')
# ggcyto(getData(gs), aes(x='Lu175Di', y='Yb173Di')) + geom_hex(bins=64)

# Gate mappings
# CD123+ HLADR- = Basophils


cols_cluster <- cols_marker %>% 
  discard(~. %in% c('BCKG190Di')) %>% 
  map(~getChannelMarker(getData(gs)[[1]], .)) %>% 
  bind_rows %>% pull('desc')


# FlowSOM

library(FlowSOM)

set.seed(1)

cols_cd4 <- c(
  # Lineage + differentation
  '160Gd_CCR7', '162Dy_CD45RA', '176Yb_CD25', '165Ho_CD127',
  
  # Th vs Tfh
  '158Gd_CXCR5', '159Tb_CXCR3', '155Gd_CCR6',
  
  # Stem cells
  '152Sm_CD27', '167Er_CD28',
  
  # Th1 v Th2
  '156Gd_CD94',
  
  # Treg subset
  '170Er_CD161', '169Tm_ICOS',
  
  # Activation and function
  '175Lu_HLADR', '151Eu_CD38',
  
  # Chronic + aging
  '113In_CD57', '147Sm_CD85j', '172Yb_PD-1'
)
cols_cd4_stem <- c('152Sm_CD27', '167Er_CD28', '160Gd_CCR7', '162Dy_CD45RA')

#fr_som <- as(getData(gs[1:3]), 'flowFrame')
#fr_som <- getData(gs)[[1]]
#fr_som <- getData(gs, '/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+/CD4+ CD8-')[['081216-Mike-HIMC ctrls-001_01_1.fcs_250000']]
fr_som <- getData(gs, '/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+/CD4+ CD8-')[['081216-Mike-HIMC ctrls-Pat03_01_1.fcs_250000']]
get_marker_indexes <- function(x) which(fr_som@parameters@data$desc %in% x)
prep_query <- function(q) names(q) %>% modify(~getChannelMarker(fr_som, .)) %>% bind_rows %>% pull('name') %>% set_names(q, .)
#mi_som <- get_marker_indexes(cols_cluster)
mi_som <- get_marker_indexes(cols_cd4)

# Prepare the flowFrame by scaling intensities (-mean/sd)
fi_som <- ReadInput(fr_som, scale=T)

# Use a 64 node (8x8) SOM instead of the default 100 node (10x10) SOM 
# * 100 is too high for visualization of ~30 latent clusters
bi_som <- BuildSOM(fi_som, colsToUse=mi_som, xdim=14, ydim=14)

# Build the MST connecting the SOM nodes
mst_som <- BuildMST(bi_som, tSNE=T)


#markers <- get_marker_indexes(c('143Nd_CD4', '144Nd_CD8', '150Nd_CD3', '162Dy_CD45RA'))
# markers <- get_marker_indexes(c(
#   '150Nd_CD3', '143Nd_CD4', '144Nd_CD8', '162Dy_CD45RA', '176Yb_CD25', 
#   '172Yb_PD-1', '142Nd_CD19', '158Gd_CXCR5', '155Gd_CCR6'
# #   ))

# # Tfh
# markers <- get_marker_indexes(c(
#   '162Dy_CD45RA', '176Yb_CD25', 
#   '158Gd_CXCR5', '155Gd_CCR6',
#   '159Tb_CXCR3', '165Ho_CD127'
# ))

#markers <- get_marker_indexes(cols_cd4_stem)
markers <- get_marker_indexes(cols_cd4)

# mst_som <- UpdateNodeSize(mst_som) # Size by count
mst_som <- UpdateNodeSize(mst_som, reset=TRUE) # Reset sizes (to 1?)
mst_som$MST$size <- mst_som$MST$size/1.5
PlotStars(mst_som, view="MST", markers=markers, colorPalette=rainbow, legend=TRUE) 
PlotStars(mst_som, view="tSNE", markers=markers, colorPalette=rainbow, legend=TRUE) 

#query <- prep_query(c("162Dy_CD45RA"='low', "176Yb_CD25"='low', "158Gd_CXCR5"='low', '165Ho_CD127'='low'))
#query <- prep_query(c("162Dy_CD45RA"='low', "160Gd_CCR7"='high'))
#query <- prep_query(c("162Dy_CD45RA"='low', "160Gd_CCR7"='low'))
query <- prep_query(c("162Dy_CD45RA"='high', "160Gd_CCR7"='low'))
target <- 'CM'
query_res <- QueryStarPlot(mst_som, query, plot=F)
bg_values <- factor(rep("Other", mst_som$map$nNodes), levels=c("Other", target))
bg_values[query_res$selected] <- target
bg_colors <- c("#FFFFFF00", "#0000FF33")
PlotStars(mst_som, view="tSNE", markers=markers, colorPalette=rainbow, legend=TRUE, 
          backgroundValues = bg_values, backgroundColor = bg_colors) 


queries <- list(
  EM=prep_query(c("162Dy_CD45RA"='low', "160Gd_CCR7"='low')),
  EMRA=prep_query(c("162Dy_CD45RA"='high', "160Gd_CCR7"='low')),
  CM=prep_query(c("162Dy_CD45RA"='low', "160Gd_CCR7"='high')),
  NAIVE=prep_query(c("162Dy_CD45RA"='high', "160Gd_CCR7"='high'))
)
bg_values <<- factor(rep("Other", mst_som$map$nNodes), levels=c("Other", names(queries)))
for (qn in names(queries)){
  query_res <- QueryStarPlot(UpdateNodeSize(mst_som, reset=TRUE), queries[[qn]], plot=F)
  bg_values[query_res$selected] <- qn
}
bg_colors <- c("#FFFFFF00", "#0000FF33", "#00FF0033", "#FF000033", "#FFFF0033")
PlotStars(mst_som, view="tSNE", markers=markers, colorPalette=rainbow, legend=TRUE, 
          backgroundValues = bg_values, backgroundColor = bg_colors) 


#PlotNumbers(mst_som)
#PlotMarker(mst_som, get_marker_indexes("144Nd_CD8"))


#nodes <- getDescendants(gs[[1]], '/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+', path='auto')
nodes <- getDescendants(gs[[1]], '/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+/CD4+ CD8-', path='auto')

fr_gates <- nodes %>% map(~getIndices(gs[[1]], .)) %>% bind_cols %>% 
  set_names(nodes) %>% apply(., 1, function(x){
    if (sum(x) < 1) NA_character_
    else {
      types <- nodes[x]
      types[rev(order(str_count(types, '/')))[1]]
    }
  })
#table(fr_gates)

n_clusters <- length(unique(fr_gates))
meta_clustering <- metaClustering_consensus(mst_som$map$codes, k=n_clusters)
dev.new(width=12, height=12)
PlotPies(mst_som, fr_gates, view="MST", backgroundValues = as.factor(meta_clustering))
PlotPies(mst_som, fr_gates, view="tSNE", backgroundValues = as.factor(meta_clustering))



# Scaffold Analysis

# export_dir <- '/tmp/human-pbmc-cytof'
# fs <- getData(gs[1:3])
# write.flowSet(fs, outdir=file.path(export_dir, 'fcs'))
# 
# nodes <- getDescendants(gs[[1]], '/Cells/Intact cells/Intact singlets/Live intact singlets/CD14- CD33-/CD3+', path='auto')
# # fr <- as(flowWorkspace::getData(gs[1:3], nodes[1]), 'flowFrame')
# to_filename <- function(n) n %>% str_replace_all(' ', '') %>% str_replace_all('/', ':') %>% paste0('landmark_', ., '.fcs')
# for (n in nodes){
#   fr <- as(flowWorkspace::getData(gs[1:3], n), 'flowFrame')
#   write.FCS(fr, file.path(export_dir, 'landmarks', to_filename(n)))
# }

# grappolo::cluster_fcs_files_groups(
#   list(all=list.files(file.path(export_dir, 'fcs'), full.names = T)),
#   num.cores = 1,
#   col.names = cols_cluster,
#   num.clusters = 250,
#   asinh.cofactor = NULL,
#   downsample.to = 30000,
#   output.dir = file.path(export_dir, 'clustered')
# )
# 
# landmarks_data <- vite::load_landmarks(
#   list.files(file.path(export_dir, 'landmarks'), full.names = T, pattern = ".*\\.fcs$"),
#   asinh.cofactor = NULL, transform.data = F
# )
# 
# cluster_files <- list.files(file.path(export_dir, 'clustered'), pattern='.*\\.clustered\\.txt', full.names = T)
# set.seed(1)
# vite::run_scaffold_analysis(
#   cluster_files,
#   ref.file=cluster_files[[1]],
#   landmarks_data,
#   cols_cluster,
#   downsample.to = 10000,
#   out.dir=file.path(export_dir, 'scaffold')
#   #inter.cluster.col.names=
#   #inter.cluster.weight.factor=.7,
#   #min.similarity=.9,
#   #overlap.method = "repel"
# )
# 
# panorama::panorama()

