library(flowCore)
library(ggcyto)
library(tidyverse)

source('cd8-exhaustion/common.R')
data_dir <- file.path(data_dir, 'fcs')

fcs_files <- list.files(data_dir, '.*HC_.*\\.fcs', full.names = TRUE)
#fcs_files <- list.files(data_dir, '.*\\.fcs', full.names = TRUE)[c(1:3, 7:12, 25:30)]
fs <- read.flowSet(fcs_files)
col_names <- col_names_dist

trans <- transformList(col_names, arcsinhTransform(a=0, b=1/5, c=0))

# Create transformed flowSet for visualization only
fst <- transform(fs, trans)
ggcyto(fst, aes(x=`PD-1`)) + geom_histogram(bins=32)
ggcyto(fst, aes(x='CD45RA', y='CCR7')) + geom_hex(bins=32) + xlim(0, 5) + ylim(0, 5) 
#ggcyto(fst, aes(x='CD45RO', y='CCR7')) + geom_hex(bins=32) + xlim(0, 5) + ylim(0, 5) 

g1 <- rectangleGate(filterId='PD-1-', `PD-1`=c(-Inf, 1)) 
g2 <- rectangleGate(filterId='PD-1+', `PD-1`=c(1, Inf))
g3 <- quadGate(filterId="Type", .gate=c("CD45RA"=1, "CCR7"=1)) 
gate <- g2 | (g1 & g3)


# Show PD-1 gate
ggcyto(fst, aes(x=`PD-1`)) + geom_histogram(bins=64) + xlim(0, 5) + geom_gate(g2) + geom_stats()

# Show PD-1- subsets
ggcyto(Subset(fst, g1), aes(x='CD45RA', y='CCR7')) + geom_hex(bins=32) + xlim(0, 5) + ylim(0, 5) + geom_gate(g3) + geom_stats()

# Split into 5 groups
fs_pop <- split(fs, (g2 | (g1 & g3)) %on% trans)

names(fs_pop) <- c(
  `CD45RA+CCR7+`='TN', 
  `CD45RA-CCR7+`='TCM', 
  `CD45RA+CCR7-`='TEMRA', 
  `CD45RA-CCR7-`='TEM', 
  `PD-1+`='PD1'
)[names(fs_pop)]
stopifnot(all(!is.na(names(fs_pop))))

for (grp_name in names(fs_pop)){
  path <- file.path(data_dir, 'landmarks', sprintf('landmarks_%s.fcs', grp_name))
  if (!dir.exists(dirname(path)))
    dir.create(dirname(path))
  # Squash flowSet into single flowFrame
  fr <- as(fs_pop[[grp_name]], "flowFrame")
  message('Writing fcs file for population ', grp_name, ' with size ', nrow(fr), ' to path ', path)
  write.FCS(fr, path)
}
