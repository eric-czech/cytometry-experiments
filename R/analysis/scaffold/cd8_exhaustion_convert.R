library(flowCore)
library(stringr)
library(purrr)
library(dplyr)
library(readr)
data_dir <- '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/bengsch_2018'


files <- tibble(
    filename=list.files(data_dir, 'CD8_.*\\.fcs\\.txt$'), 
    group=str_split(filename, '_') %>% map(~.[[2]]) %>% as.character,
    filepath=file.path(data_dir, filename)
  ) %>% 
  dplyr::filter(!str_detect(filename, 'batch'))

# to_df <- function(path) read_delim(path, delim='\t', skip=1) %>% select(-one_of('Event #'))
# to_fr <- function(df) flowFrame(as.matrix(df))
# fs <- files$filepath %>% map(to_df) %>% map(to_fr) %>% flowSet
# cts <- fs %>% map(colnames) %>% unlist %>% table %>% sort
# cts
# tet_CMV   tet_FLU_A2       CD200R        Gzm_M       Ptger2       HLA-DR         TCF1  tet_HIV_ILK  tet_HIV_SL9         ICOS 
# 7            8           12           12           12           18           18           18           21           23 
# KLRG1         CD36        Foxp3     PERFORIN        CD103          2B4        beads         CCR7      Cd112Di      Cd114Di 
# 26           32           32           32           37           44           44           44           44           44 
# CD127         CD16        CD160         CD19         CD26         CD27         CD28          CD3         CD38         CD39 
# 44           44           44           44           44           44           44           44           44           44 
# CD4         CD45       CD45RA       CD45RO         CD57          CD7         CD73          CD8    Cisplatin        CTLA4 
# 44           44           44           44           44           44           44           44           44           44 
# CXCR5        Eomes Event_length    GranzymeB    GranzymeK       Helios   Iridium191   Iridium193         Ki67         LAG3 
# 44           44           44           44           44           44           44           44           44           44 
# LD    PCA_1_ESG    PCA_2_ESG         PD-1  Pheno_ESG_1  Pheno_ESG_2         Tbet        TIGIT         Tim3         Time 
# 44           44           44           44           44           44           44           44           44           44 
# TOX 
# 44
# cts[cts == max(cts)] %>% names %>% paste(collapse='", "') %>% cat

# Set list of fields present in all panels to be used to standardize all FCS files
common_fields <- c(
  "2B4", "beads", "CCR7", "Cd112Di", "Cd114Di", "CD127", "CD16", "CD160", "CD19", "CD26", "CD27", "CD28", "CD3", 
  "CD38", "CD39", "CD4", "CD45", "CD45RA", "CD45RO", "CD57", "CD7", "CD73", "CD8", "Cisplatin", "CTLA4", "CXCR5", 
  "Eomes", "Event_length", "GranzymeB", "GranzymeK", "Helios", "Iridium191", "Iridium193", "Ki67", "LAG3", "LD", 
  "PCA_1_ESG", "PCA_2_ESG", "PD-1", "Pheno_ESG_1", "Pheno_ESG_2", "Tbet", "TIGIT", "Tim3", "Time", "TOX"
)

# Define functions for reading Cytobank csv files as well as convertion to flowFrame
to_df <- function(path) read_delim(path, delim='\t', skip=1) %>% select(one_of(common_fields))
to_fr <- function(df) flowFrame(as.matrix(df))

# Convert each csv to a flowFrame and group all flowFrames into a flowSet
fs <- files$filepath %>% map(to_df) %>% map(to_fr) %>% flowSet

# Define asinh transform for all expression fields
# trans <- transformList(common_fields %>% discard(~. %in% c('beads', 'Event_length', 'Time')), arcsinhTransform())
# fst <- transform(fs, trans)

# Save results as individual FCS files with common expression fields
write.flowSet(fs, outdir=file.path(data_dir, 'fcs'), str_replace(files$filename, '.txt', ''))

# fr <- read.FCS(file.path(data_dir, 'fcs', 'CD8_HC_001.fcs'))
