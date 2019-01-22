
data_dir <- '/Users/eczech/repos/cytometry-experiments/R/analysis/scaffold/data/bengsch_2018'
col_names_shared <- c(
  "2B4", "beads", "CCR7", "Cd112Di", "Cd114Di", "CD127", "CD16", "CD160", "CD19", "CD26", "CD27", "CD28", "CD3", 
  "CD38", "CD39", "CD4", "CD45", "CD45RA", "CD45RO", "CD57", "CD7", "CD73", "CD8", "Cisplatin", "CTLA4", "CXCR5", 
  "Eomes", "Event_length", "GranzymeB", "GranzymeK", "Helios", "Iridium191", "Iridium193", "Ki67", "LAG3", "LD", 
  "PCA_1_ESG", "PCA_2_ESG", "PD-1", "Pheno_ESG_1", "Pheno_ESG_2", "Tbet", "TIGIT", "Tim3", "Time", "TOX"
)
col_names_dist <- c(
  "2B4", "CCR7", "Cd112Di", "Cd114Di", "CD127", "CD16", "CD160", "CD19", "CD26", "CD27", "CD28", "CD3", 
  "CD38", "CD39", "CD45RA", "CD45RO", "CD57", "CD7", "CD73", "Cisplatin", "CTLA4", "CXCR5", 
  "Eomes", "GranzymeB", "GranzymeK", "Helios", "Ki67", "LAG3", "PD-1", "Tbet", "TIGIT", "Tim3", "TOX"
)

# tSNE/PhenoGraph params from Table S4 in 
# https://www.cell.com/cms/10.1016/j.immuni.2018.04.026/attachment/eff6a6bf-f812-4cf9-b3f1-735588df23b7/mmc1.pdf
col_names_core <- c(
  "CTLA4", "CD7", "CD73", "CD127", "CD39", "GranzymeK", "Helios", "PD-1", "Eomes", "CD38", "TOX",
  "TIGIT", "CXCR5", "2B4", "LAG3", "CCR7", 
  "CD45RA" # CD45RA not in doc but included here for consistency
)

# This is parameterized to be equivalent to (default) transformation used in ParkerICI/vite
arcsinh_trans <- arcsinhTransform(a=0, b=1/5, c=0)

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
