library(flowCore)
library(tidyverse)
library(vite)
library(grappolo)

data_dir <- '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/bengsch_2018/fcs'
fcs_files <- list.files(data_dir, '.*\\.fcs')[c(1:3, 7:9, 25:27)]
col_names <- c(
  "2B4", "CCR7", "Cd112Di", "Cd114Di", "CD127", "CD16", "CD160", "CD19", "CD26", "CD27", "CD28", "CD3", 
  "CD38", "CD39", "CD4", "CD45", "CD45RA", "CD45RO", "CD57", "CD7", "CD73", "CD8", "Cisplatin", "CTLA4", "CXCR5", 
  "Eomes", "GranzymeB", "GranzymeK", "Helios", "Ki67", "LAG3", 
  "PD-1", "Tbet", "TIGIT", "Tim3", "Time", "TOX"
)

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
  healthy=fcs_files %>% purrr::keep(~str_detect(., 'HC')),
  hiv=fcs_files %>% purrr::keep(~str_detect(., 'HIV')),
  luca=fcs_files %>% purrr::keep(~str_detect(., 'LUCA'))
) %>% purrr::map(~file.path(data_dir, .))

grappolo::cluster_fcs_files_groups(
  fcs_file_groups, 
  num.cores = 1, 
  col.names = col_names, 
  num.clusters = 100, 
  asinh.cofactor = 5, 
  downsample.to = 1000, 
  output.dir = "output/bengsch_2018/clustered_by_group"
)

cluster_files <- list.files(path = "output/bengsch_2018/clustered_by_group", pattern = ".*\\.clustered\\.txt$", full.names = TRUE)
landmarks_data <- vite::load_landmarks_from_dir(file.path(data_dir, 'landmarks'), asinh.cofactor = NULL, transform.data = F)
vite::run_scaffold_analysis(
  cluster_files, 
  ref.file=cluster_files[[1]], 
  landmarks_data, col_names, 
  downsample.to = 10000,
  out.dir="output/bengsch_2018/clustered_by_group/scaffold"
)
panorama::panorama()

# fr <- read.FCS('/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/spitzer_2015/landmarks/BM2_cct_normalized_01_BM2_cct_normalized_concat.fcs_B cells.fcs')


