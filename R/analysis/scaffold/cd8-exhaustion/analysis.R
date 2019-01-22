library(flowCore)
library(tidyverse)
library(vite)
library(grappolo)

source('cd8-exhaustion/common.R')
data_dir <- file.path(data_dir, 'fcs')
#fcs_files <- list.files(data_dir, '.*\\.fcs')[c(1:3, 7:9, 25:27)]
fcs_files <- c(
  # Healthy
  'CD8_HC_001.fcs', 'CD8_HC_002.fcs', 'CD8_HC_003.fcs', 'CD8_HC_004.fcs', 'CD8_HC_006.fcs',
  
  # Severe untreated HIV
  'CD8_HIV_002.fcs', 'CD8_HIV_003.fcs', 'CD8_HIV_006.fcs', 'CD8_HIV_015.fcs', 'CD8_HIV_027.fcs',
  
  sprintf('CD8_LUCA_00%s_PBMC.fcs', 1:5),
  sprintf('CD8_LUCA_00%s_TIL.fcs', 1:5),
  sprintf('CD8_LUCA_00%s_Lung.fcs', 1:5)
)
col_names <- col_names_core

#### Independent Clustering

cluster_fcs_files(
  file.path(data_dir, fcs_files), 
  asinh.cofactor=5,
  num.cores=1,
  col.names=col_names,
  num.clusters=50,
  output.dir = "output/bengsch_2018/clustered_single_samples"
)

cluster_files <- dplyr::tibble(
  filename=list.files(path = "output/bengsch_2018/clustered_single_samples", pattern = ".*\\.clustered\\.txt$", full.names = TRUE),
  group=filename %>% basename %>% str_split('_') %>% purrr::map(~.[[2]]) %>% as.character
)


# Unsupervised
G <- vite::get_unsupervised_graph_from_files(
  cluster_files$filename,
  col.names = col_names,
  metadata.tab = as.data.frame(cluster_files),
  metadata.filename.col = 'filename',
  filtering.threshold = 15,
  process.clusters.data=T,
  downsample.to=1000,
  clusters.data.out.dir = "output/bengsch_2018/clustered_single_samples/unsupervised_graph"
)
# trace(vite::write_clusters_data, browser)
# untrace(vite::write_clusters_data)
vite::write_graph(G, "output/bengsch_2018/clustered_single_samples/unsupervised_graph/graph.graphml")
panorama::panorama()



# Pooled clustering

fcs_file_groups <- list(all.pooled=file.path(data_dir, fcs_files))
cluster_fcs_files_groups(
  fcs_file_groups, 
  num.cores = 1, 
  col.names = col_names, 
  num.clusters = 50, 
  asinh.cofactor = 5, 
  downsample.to = 5000, 
  output.dir = "output/bengsch_2018/clustered_all_samples"
)


cluster_files <- list.files(path = "output/bengsch_2018/clustered_all_samples", pattern = ".*\\.clustered\\.txt$", full.names = TRUE)
G <- vite::get_unsupervised_graph_from_files(
  cluster_files, 
  col.names = col_names,
  filtering.threshold = 15,
  downsample.to=1000,
  clusters.data.out.dir = "output/bengsch_2018/clustered_all_samples/unsupervised_graph"
)
vite::write_graph(G, "output/bengsch_2018/clustered_all_samples/unsupervised_graph/graph.graphml")


### Scaffold Analysis

fcs_file_groups <- list(
  healthy_pbmc=fcs_files %>% purrr::keep(~str_detect(., 'HC')),
  hiv_pbmc=fcs_files %>% purrr::keep(~str_detect(., 'HIV')),
  luca_pbmc=fcs_files %>% purrr::keep(~str_detect(., 'LUCA_.*_PBMC')),
  luca_lung=fcs_files %>% purrr::keep(~str_detect(., 'LUCA_.*_Lung')),
  luca_til=fcs_files %>% purrr::keep(~str_detect(., 'LUCA_.*_TIL'))
) %>% purrr::map(~file.path(data_dir, .))

# Spitzer data has 80 files and is clustered by 4 tissue types with 10k per file
# giving about 20 * 10,000 = 200,000 events to be separated into 200 clusters
# i.e. aim for num events = 1000 * num clusters
grappolo::cluster_fcs_files_groups(
  fcs_file_groups, 
  num.cores = 1, 
  col.names = col_names, 
  num.clusters = 75, 
  asinh.cofactor = 5, 
  downsample.to = 15000, 
  output.dir = "output/bengsch_2018/clustered_by_group"
)

cluster_files <- list.files(path = "output/bengsch_2018/clustered_by_group", pattern = ".*\\.clustered\\.txt$", full.names = TRUE)
landmarks_data <- vite::load_landmarks_from_dir(file.path(data_dir, 'landmarks'), asinh.cofactor = 5, transform.data = T)
vite::run_scaffold_analysis(
  cluster_files, 
  ref.file=cluster_files[[1]], 
  landmarks_data, col_names, 
  downsample.to = 10000,
  out.dir="output/bengsch_2018/clustered_by_group/scaffold",
  inter.cluster.weight.factor=.7
)
panorama::panorama()
# CTLA4
# /Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/output/bengsch_2018/clustered_by_group/scaffold/healthy_pbmc.graphml

# fr <- read.FCS('/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/spitzer_2015/landmarks/BM2_cct_normalized_01_BM2_cct_normalized_concat.fcs_B cells.fcs')


