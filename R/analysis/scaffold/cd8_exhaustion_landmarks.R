library(flowCore)
library(ggcyto)
library(tidyverse)

data_dir <- '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/bengsch_2018/fcs'
fcs_files <- list.files(data_dir, '.*HC_.*\\.fcs', full.names = TRUE)
#fcs_files <- list.files(data_dir, '.*\\.fcs', full.names = TRUE)[c(1:3, 7:12, 25:30)]
fs <- read.flowSet(fcs_files)

col_names <- c(
  "2B4", "CCR7", "Cd112Di", "Cd114Di", "CD127", "CD16", "CD160", "CD19", "CD26", "CD27", "CD28", "CD3", 
  "CD38", "CD39", "CD4", "CD45", "CD45RA", "CD45RO", "CD57", "CD7", "CD73", "CD8", "Cisplatin", "CTLA4", "CXCR5", 
  "Eomes", "GranzymeB", "GranzymeK", "Helios", "Iridium191", "Iridium193", "Ki67", "LAG3", "LD", 
  "PD-1", "Tbet", "TIGIT", "Tim3", "Time", "TOX"
)
trans <- transformList(col_names, arcsinhTransform(a=0, b=1/5, c=0))
fst <- transform(fs, trans)

gate <- quadGate(filterId="subset", .gate=c("CD45RO"=1, "CCR7"=1))
#ggcyto(fst, aes(x='CD45RO')) + geom_histogram(bins=64) + xlim(0, 3)
#ggcyto(fst, aes(x='CCR7')) + geom_histogram(bins=64) + xlim(0, 3)
ggcyto(fst, aes(x='CD45RO', y='CCR7')) + geom_hex(bins=32) + xlim(0, 3) + ylim(0, 4) + geom_gate(gate) + geom_stats()

fs_pop <- split(fst, gate)
names(fs_pop) <- c(`CD45RO+CCR7+`='TCM', `CD45RO-CCR7+`='TN', `CD45RO+CCR7-`='TEFF', `CD45RO-CCR7-`='TEMRA')[names(fs_pop)]
for (grp_name in names(fs_pop)){
  path <- file.path(data_dir, 'landmarks', sprintf('landmarks_%s.fcs', grp_name))
  if (!dir.exists(dirname(path)))
    dir.create(dirname(path))
  # Squash flowSet into single flowFrame
  fr <- as(fs_pop[[grp_name]], "flowFrame")
  message('Writing fcs file for population ', grp_name, ' to path ', path)
  write.FCS(fr, path)
}
