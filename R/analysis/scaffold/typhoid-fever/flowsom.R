

library(tidyverse)
library(fs)
library(flowCore)
library(ggcyto)
library(openCyto)
library(FlowRepositoryR)


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


library(FlowSOM)

fr_som <- fr_cytof
cols_cluster <- colnames(fr_som) %>% discard(~.%in%c('Cisplatin', 'Donor', 'CD3', 'CD8', 'CD4', 'CD45'))
get_marker_indexes <- function(x) which(fr_som@parameters@data$desc %in% x)
prep_query <- function(q) names(q) %>% modify(~getChannelMarker(fr_som, .)) %>% bind_rows %>% pull('name') %>% set_names(q, .)
mi_som <- get_marker_indexes(cols_cluster)

# Prepare the flowFrame by scaling intensities (-mean/sd)
set.seed(1)
fi_som <- ReadInput(fr_som, scale=T)

# Use a 64 node (8x8) SOM instead of the default 100 node (10x10) SOM 
# * 100 is too high for visualization of ~30 latent clusters
bi_som <- BuildSOM(fi_som, colsToUse=mi_som, xdim=12, ydim=12)

# Build the MST connecting the SOM nodes
mst_som <- BuildMST(bi_som, tSNE=T)

markers <- get_marker_indexes(c(
  'Vd2', 'HLA-DR', 'Va7_2', 'CLA', 'gdTCR', 'CD57', 'CD27', 'Foxp3', 'CD192',
  'CD56', 'CD16', 'CDw199', 'CD185-bio', 'PD-1', 'CD103', 'Ki67', 'CD49a',
  'ICOS', 'Intb7', 'BCL2', 'CD127',
  'CD194', 'CD183'
))
# markers <- get_marker_indexes(c(
#   'CD183', 'CD185-bio', 'CD183', 'CD196', 'CD194', 'CD25', 'CD127',
#   'CD45RO', 'CD45RA', 'CD197'
# ))
queries <- list(
  EM=prep_query(c("CD45RA"='low', "CD197"='low')),
  EMRA=prep_query(c("CD45RA"='high', "CD197"='low')),
  CM=prep_query(c("CD45RA"='low', "CD197"='high')),
  NAIVE=prep_query(c("CD45RA"='high', "CD197"='high'))
)
bg_values <<- factor(rep("Other", mst_som$map$nNodes), levels=c("Other", names(queries)))
for (qn in names(queries)){
  query_res <- QueryStarPlot(UpdateNodeSize(mst_som, reset=TRUE), queries[[qn]], plot=F)
  bg_values[query_res$selected] <- qn
}
bg_colors <- c("#FFFFFF00", "#0000FF33", "#00FF0033", "#FF000033", "#FFFF0033")

mst_som <- UpdateNodeSize(mst_som)
mst_som$MST$size <- mst_som$MST$size/2
PlotStars(mst_som, view="MST", markers=markers, colorPalette=rainbow, legend=TRUE, 
          backgroundValues = bg_values, backgroundColor = bg_colors) 

# mst_som <- UpdateNodeSize(mst_som, reset=TRUE) # Reset sizes (to 1?)
# mst_som$MST$size <- mst_som$MST$size/1.5
PlotMarker(mst_som, get_marker_indexes('CD25'), view='MST' )
# mst_som <- UpdateNodeSize(mst_som)

par(1)
#par(mfrow=c(2,2))
marker_color <- c(
  'Vd2', 'HLA-DR', 'Va7_2', 'CLA', 'gdTCR', 'CD57', 'CD27', 'Foxp3', 'CD192', 
  'CD56', 'CD16', 'CDw199', 'CD185-bio', 'PD-1', 'CD103', 'Ki67', 'CD49a',
  'ICOS', 'Intb7', 'BCL2', 'CD127'
)
for (m in marker_color){
  print(paste0('/tmp/som_', m, '.png'))
  png(paste0('/tmp/som_', m, '.png'))
  PlotMarker(mst_som, get_marker_indexes(m), view='MST', main=m)
  dev.off()
}
dev.off()

