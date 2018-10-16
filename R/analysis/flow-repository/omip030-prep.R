# install.packages("BiocManager")
# BiocManager::install('flowDensity')
library(magrittr)
library(tidyverse)
library(fs)
library(FlowRepositoryR)
library(flowCore)
source('lib.R')

# Dataset: https://flowrepository.org/id/FR-FCM-ZZWU
data_dir <- dir_create("/tmp/flow/")
dataset_id <- "FR-FCM-ZZWU"
dataset_path <- fs::path(data_dir, dataset_id)
dataset <- flowRep.get(dataset_id)

if (!dir_exists(dataset_path)) 
  dataset %<>% download(dirpath=dataset_path, only.files=".*fcs", show.progress=FALSE)
fr_raw <- fs::path(dataset_path, 'PBMC_control 1-4_cct.fcs') %>% read.FCS %>% clean_flow_frame
fr_raw@description$SPILL <- fr_raw@description$`$SPILLOVER`

path <- '/home/eczech/analysis/data/omip/omip-030/PBMC_control1-4_cct_clean.fcs'
write.FCS(fr_raw, path)

# Read and rewrite resulting file so that file paths stored represent original source
fr_clean <- read.FCS(path)
write.FCS(fr_clean, path)
