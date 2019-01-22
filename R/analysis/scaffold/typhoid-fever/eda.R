

##### Marker Subsets

library(flowCore)
library(tidyverse)
fr <- read.FCS('~/Downloads/SurfacePanel_1-D0')
#fr <- read.FCS('~/Downloads/SurfacePanel_2-D0')
#fr <- read.FCS('~/Downloads/IntracellularCytokinePanel_1-D0')

# See: https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22788
omip30 <- c('CCR10', 'CD127', 'CD161', 'CXCR3', 'CXCR5', 'CCR4', 'CCR6', 'CCR7', 'CD20', 'CD25', 'CD3', 
            'CD38', 'CD4', 'CD45RA', 'CD8', 'Ki67')

panel <- fr %>% markernames %>% sort

# set1 <- c('alphaGalcer', 'BCL2', 'CCR2', 'CCR4', 'CCR5', 'CCR6', 'CCR7', 'CCR9', 'CD103', 'CD127', 'CD14', 'CD14', 'CD16', 'CD161', 'CD19', 
#           'CD25', 'CD27', 'CD3', 'CD38', 'CD4', 'CD45', 'CD45RA', 'CD45RO', 'CD49a', 'CD49d', 'CD56', 'CD57', 'CD8', 
#           'Cisplatin', 'CLA', 'CX3CR1', 'CXCR3', 'CXCR5', 'DNA', 'DNA', 'Foxp3', 'gdTCR', 'Granzyme_B', 'HLA-DR', 
#           'ICOS', 'Intb7', 'Ki67', 'PD-1', 'Va7_2', 'Vd1', 'Vd2')


cat(setdiff(panel, omip30))
# for SurfacePanel:
# alphaGalcer BC BCL2 CCR2 CCR5 CCR9 CD103 CD14 CD16 CD19 CD27 CD45 CD45RO CD49a CD49d 
# CD56 CD57 Cisplatin CLA CX3CR1 DNA Foxp3 gdTCR Granzyme_B HLA-DR ICOS Intb7 PD-1 Va7_2 Vd1 Vd2
cat(setdiff(panel, omip30))
# for IntracellularCytokinePanel:
# alphaGalcer BC CCR2 CCR9 CD103 CD107a CD14 CD16 CD19 CD40L CD45 CD56 CD57 Cisplatin CTLA-4 CX3CR1 DNA 
# gdTCR GM-CSF HLA-DR IFNg IL-10 IL-17A IL-2 IL-22 IL-4 IL-8 Intb7 MIP-1b PD-1 TNFa Va7.2 Vd1 Vd2
cat(setdiff(omip30, panel))
# OMIP30 - SurfacePanel: CCR10 CD20
# OMIP30 - IntracellularCytokinePanel: CCR10 CD127 CXCR3 CD20 Ki67

##### CyTOF v Flow Comparison

library(tidyverse)
library(fs)
library(flowCore)
library(ggcyto)
library(openCyto)
library(FlowRepositoryR)

get_omip30_data <- function(){
  data_dir <- dir_create('~/repos/hammer/t-cell-data/data/flow_files/OMIP-030')
  zip_url <- 'https://storage.googleapis.com/t-cell-data/flow_files/OMIP-030.zip'
  zip_file <- fs::path(data_dir, 'OMIP-030.zip')
  if (!file_exists(zip_file)) {
    curl_download(zip_url, zip_file) %>% unzip(exdir=data_dir)
  }
  ws <- openWorkspace(fs::path(data_dir, 'OMIP-030.wsp'))
  parseWorkspace(ws, name='TCells', path=data_dir, isNcdf=FALSE)
}
gs <- get_omip30_data()

# fr_flow <- getData(gs, 'T Cells')[[1]]
#frt <- read.FCS('/Users/eczech/repos/cytometry-experiments/data/omip/omip-030/downloads/PBMC_control1-4_cct.fcs')
# frt <- read.FCS('/Users/eczech/repos/hammer/t-cell-data/data/flow_files/OMIP-030/PBMC_control1-4_cct_clean.fcs')

# Extract T Cell subset and invert scales back to original
#gh_flow <- transform(getData(gs, 'T Cells')[[1]], transformList(names(gs_inv_trans), gs_inv_trans))
#getNodes(gs[[1]], '/T Cell')
node_paths <- getDescendants(gs, 'T Cells', path='auto') %>% keep(~str_detect(., '^CD[4|8]\\+'))

fs_flow <- map(node_paths, ~ getData(gs, .)[[1]]) %>% set_names(node_paths) %>% as('flowSet')
gs_inv_trans <- getTransformations(gs[[1]], inverse=TRUE)
fs_flow <- transform(fs_flow, transformList(names(gs_inv_trans), gs_inv_trans))
fs_flow <- transform(fs_flow, transformList(colnames(fs_flow), arcsinhTransform(a=0, b=1/150, c=0)))

cols_flow <- fs_flow[[1]]@parameters@data %>% 
  mutate(desc=as.character(desc)) %>%
  dplyr::filter(!is.na(desc)) %>%
  dplyr::filter(desc != 'L+D') %>%
  (function(df) set_names(df$desc, df$name)) 

fs_flow <- fs_flow[,names(cols_flow)]
markernames(fs_flow) <- cols_flow 
colnames(fs_flow) <- cols_flow


get_cytof_data <- function(dataset_id, pattern='.*'){
  # Note that download will choke without a full path (no ~)
  data_dir <- dir_create("/Users/eczech/repos/hammer/t-cell-data/data/cytof_files")
  dataset_path <- fs::path(data_dir, dataset_id)
  dataset <- flowRep.get(dataset_id)
  if (!dir_exists(dataset_path))
    dataset <- download(dataset, dirpath=dataset_path, only.files=pattern, show.progress=F)
  read.flowSet(list.files(dataset_path, pattern=pattern, full.names = T))
}
fs <- get_cytof_data('FR-FCM-ZYW2', pattern='SurfacePanel_.*-D0')
#fs <- get_cytof_data('FR-FCM-ZYW2', pattern='SurfacePanel_.*-D28')
#fs <- get_cytof_data('FR-FCM-ZYW2', pattern='SurfacePanel_.*-D4')

#trans <- transformList(colnames(fs), arcsinhTransform(a=0, b=1/5, c=0))

# ggcyto(fs, aes(x='CCR7')) + geom_histogram(bins=64) + scale_x_logicle()
# ggcyto(transform(fs, trans), aes(x='CD38')) + geom_histogram(bins=64) 
# ggcyto(fs, aes(x='Event_length')) + geom_histogram(bins=64) 
# 
# ggcyto(getData(gs[[1]]), aes(x='CD197')) + geom_histogram(bins=64) + scale_x_logicle()
# ggcyto(getData(gs[[1]]), aes(x='CD45RA')) + geom_histogram(bins=64) + scale_x_logicle()



fr_cytof <- as(fs, 'flowFrame')

cols_cytof <- c(
  CCR2='CD192',
  CCR4='CD194',
  CCR5='CD195',
  CCR6='CD196',
  CCR7='CD197',
  CCR9='CDw199',
  CXCR3='CD183',
  CXCR5='CD185-bio',
  `Original Frame`='Donor'
) 
cols_cytof <- fr_cytof@parameters@data %>% 
  mutate(marker=map_if(desc, desc %in% names(cols_cytof), ~ unname(cols_cytof[.]))) %>%
  dplyr::filter(!is.na(marker)) %>%
  dplyr::filter(!duplicated(marker)) %>%
  dplyr::filter(!(marker %in% c('DNA', 'BC'))) %>%
  (function(df) set_names(df$marker, df$name)) %>% unlist

fr_cytof <- fr_cytof[,names(cols_cytof)]
markernames(fr_cytof) <- cols_cytof 
colnames(fr_cytof) <- cols_cytof

cytofAsinh <- function(transformationId = "defaultCytofAsinh", cofactor = 5, sd=.01){
  t <- new("transform", .Data = function(value) {
    value <- value-1
    loID <- which(value < 0)
    if(length(loID) > 0)
      value[loID] <- rnorm(length(loID), mean = 0, sd = sd)
    value <- value / cofactor
    value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
    return(value)
  })
  t@transformationId <- transformationId
  t
}
#fr_cytof <- transform(fr_cytof, transformList(discard(colnames(fr_cytof), ~.=='Donor'), arcsinhTransform(a=0, b=1/5, c=0)))
fr_cytof <- transform(fr_cytof, transformList(discard(colnames(fr_cytof), ~.=='Donor'), cytofAsinh(cofactor=1, sd=.5)))

# ggcyto(fr_cytof, aes(x='CD3')) + geom_histogram(bins=64) + xlim(.01, 6) 
# #ggcyto(as(fs_flow, 'flowFrame'), aes(x='Ki67')) + geom_histogram(bins=64) + xlim(-2, 6)
# ggcyto(transform(fs_flow[['CD4+',]], transformList(colnames(fs_flow), truncateTransform(a=0))), aes(x='CD3')) + 
#   geom_histogram(bins=64) + xlim(.01, 6)

cols_shared <- base::intersect(colnames(fs_flow), colnames(fr_cytof))
# fs_all <- list(
#   flow=fs_flow[['CD4+',cols_shared]],
#   cytof=fr_cytof[,cols_shared]
# ) %>% as('flowSet')
fs_all <- as(fs_flow[,cols_shared], 'list')
fs_all$cytof <- fr_cytof[,cols_shared]
fs_all <- as(fs_all, 'flowSet')

# flowViz::densityplot(~., fs_all)
#flowViz::histogram(~., fs_all)
fs_norm <- flowStats::warpSet(fs_all, cols_shared, target='cytof')
# flowViz::densityplot(~., fs_norm)
flowViz::densityplot(~., fs_norm[c('CD4+ FH Th1-like', 'CD4+ Naive Treg', 'CD8+ EM', 'CD4+ Naive', 'CD4+ Th22', 'cytof'),])

#warp_funs <- warpSet(fs_all, cols_shared, warpFuns = T)

# Use normalized flow data as-is
fs_flow_norm <- fs_norm[sampleNames(fs_norm) %>% discard(~.=='cytof'),]

# Combine normalized cytof data with extra markers
fr_cytof_norm <- cbind(fr_cytof[,!colnames(fr_cytof) %in% cols_shared], exprs(fs_norm[['cytof']]))
markernames(fr_cytof_norm) <- set_names(colnames(fr_cytof_norm), colnames(fr_cytof_norm))
colnames(fr_cytof_norm) <- markernames(fr_cytof_norm)

# Export
write.FCS(fr_cytof_norm, '/tmp/cd4/cytof.fcs')
write.flowSet(fs_flow_norm, outdir='/tmp/cd4/landmarks', filename=sampleNames(fs_flow_norm) %>% str_replace_all('[\\+| |\\-]+', ''))

# write.FCS(fr_cytof, '/tmp/cd4/cytof.fcs')
# write.flowSet(fs_flow, outdir='/tmp/cd4/landmarks', filename=sampleNames(fs_flow) %>% str_replace_all('[\\+| |\\-]+', ''))


# Scaffold Analysis

col_names_cluster <- colnames(fr_cytof_norm) %>% discard(~.%in%c('Cisplatin', 'Donor', 'CD3', 'CD8', 'CD4', 'CD45'))
grappolo::cluster_fcs_files_groups(
  list(all='/tmp/cd4/cytof.fcs'), 
  num.cores = 1, 
  col.names = col_names_cluster, 
  num.clusters = 100, 
  asinh.cofactor = NULL, 
  downsample.to = 100000, 
  output.dir = "/tmp/cd4/output/clustered"
)

landmarks_data <- vite::load_landmarks(
  #list.files('/tmp/cd4/landmarks', full.names = T, pattern = ".*fcs") %>% keep(~str_detect(., 'CD4(FH|Th|CM|EM|Naive|Eff|Mem)')), 
  list.files('/tmp/cd4/landmarks', full.names = T, pattern = ".*fcs"),
  asinh.cofactor = NULL, transform.data = F
)

cluster_files <- list.files('/tmp/cd4/output/clustered', pattern='.*\\.clustered\\.txt', full.names = T)
set.seed(1)
vite::run_scaffold_analysis(
  cluster_files, 
  ref.file=cluster_files[[1]], 
  landmarks_data, 
  cols_shared, 
  downsample.to = 10000,
  out.dir="/tmp/cd4/output/scaffold",
  inter.cluster.col.names=cols_shared,#col_names_cluster,
  inter.cluster.weight.factor=.1,
  min.similarity=.9,
  overlap.method = "repel"
)

panorama::panorama()


