library(flowCore)
library(ggcyto)
library(tidyverse)

source('cd8-exhaustion/common.R')
data_dir <- file.path(data_dir, 'fcs')

# Read in data for all healthy donors
fs <- list.files(data_dir, full.names = T) %>% purrr::keep(~str_detect(., 'HC_')) %>% read.flowSet

# Apply transforms to most fields
#fs <- transform(fs, transformList(col_names_dist, arcsinh_trans))
fst <- transform(fs, transformList(col_names_dist, cytofAsinh(sd=1)))

gate <- quadGate("CCR7"=1, "2B4"=1)


fsp <- fst
#fsp <- Subset(fst, rectangleGate('CD45RA'=c(1, Inf))) # Check split in CD45RA+ cells only
ggcyto(fsp, aes(x='CCR7', y='2B4')) + 
  geom_hex(bins=64) + geom_gate(gate) + geom_stats() +
  xlim(-1, 5) + ylim(-1, 4) +
  labs(x='CCR7', y='CD244 (aka 2B4)') + 
  ggtitle('Mutual Exclusivity of CD244 and CCR7 in 6 Healthy Humans\n/Live/Singlet/CD45+/CD3+/CD8+')