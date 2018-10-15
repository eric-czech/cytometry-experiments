# install.packages("BiocManager")
# BiocManager::install('flowDensity')
library(magrittr)
library(tidyverse)
library(fs)
library(FlowRepositoryR)
library(flowCore)
library(openCyto)
library(ggcyto)

map <- purrr::map
add <- flowCore::add
filter <- dplyr::filter
source('lib.R')

# Dataset: https://flowrepository.org/id/FR-FCM-ZZWU
data_dir <- dir_create("/tmp/flow/")
dataset_id <- "FR-FCM-ZZWU"
dataset_path <- fs::path(data_dir, dataset_id)
dataset <- flowRep.get(dataset_id)

if (!dir_exists(dataset_path)) 
  dataset %<>% download(dirpath=dataset_path, only.files=".*fcs", show.progress=FALSE)
fr_raw <- fs::path(dataset_path, 'PBMC_control 1-4_cct.fcs') %>% read.FCS %>% clean_flow_frame

# Find field names in first non-null compensation matrix
# fluoro_channels <- fr_raw %>% spillover %>% discard(is_null) %>% extract2(1) %>% colnames
# frame_transform <- transformList(fluoro_channels, biexponentialTransform())
# frame_transform <- transformList(fluoro_channels, logicleTransform())
# fr <- transform(fr_raw, frame_transform)

# Keyword searching 
# fr %>% keyword %>% names %>% str_extract('.*P.*G$') %>% discard(is.na) %>% 
#  magrittr::extract(fr@description, .) %>% as_vector

fr_comp <- compensate(fr, fr@description$`$SPILLOVER`)
ggcyto(fr, aes(x='FSC_A', y='FSC_W')) + geom_hex(bins=100)
autoplot(fr, 'CD4', bins=100, axis_inverse_trans=T) + xlim(2, 7) + ylim(2, 8)

gate_live <- draw_gate(fr, c('FSC_A', 'FSC_W'))
#g <- flowClust.2d(fr, 'FSC_A', 'FSC_W', quantile = .95)
#g <- quadGate('FSC_A'=12, 'FSC_W'=12)
ggcyto(fr, aes(x='FSC_A', y='FSC_W')) + geom_hex(bins=100) + geom_gate(gate_live) + geom_stats()

Subset(fr, gate_live)

get_gate_coord_string(gate_live)
gate_live@boundaries %>% as.data.frame %>% c


### TSNE

library(Rtsne)
fr_samp <- fr@exprs[,fluoro_channels] %>% as_tibble %>% sample_n(2500)
tsne <- Rtsne(as.matrix(fr_samp))
plot(tsne$Y)

### FlowJo

?flowWorkspace::openWorkspace
?flowWorkspace::parseWorkspace

#ws <- flowWorkspace::openWorkspace('~/analysis/data/omip/omip-030/OMIP_PBMC_example_template.wsp')
#flowWorkspace::closeWorkspace(ws)
ws <- flowWorkspace::openWorkspace('~/analysis/data/omip/omip-030/clean/OMIP_PBMC_clean.wsp')
gs <- flowWorkspace::parseWorkspace(
  ws, name='TCells', path='/home/eczech/analysis/data/omip/omip-030/clean',
  min.limit=NULL, isNcdf=F, alter.names=T
)
gh <- gs[[1]]

# Show plot from t cells on
plot(gs, 'T Cells', width=6, height=1, fontsize=18, shape='rectangle')
# plot(gs, width=6, height=1, fontsize=18, shape='rectangle')

# Show CD8+ tree
# plot(gs, 'CD8+', shape='rectangle')
# autoplot(gs[[1]], getDescendants(gs, 'CD8+', path=1))

all_nodes <- getNodes(gs, path=1)
cell_nodes <- getNodes(gs, path=1) %>% str_match(., 'CD[4|8]\\+ .*') %>% discard(is.na)
cell_class <- cell_nodes %>% str_detect('CD4+') %>% ifelse('CD4+', 'CD8+') %>% set_names(cell_nodes)
#cd4_nodes <- cell_nodes %>% str_match(., 'CD4\\+ .*') %>% discard(is.na)
#cd8_nodes <- cell_nodes %>% str_match(., 'CD8\\+ .*') %>% discard(is.na)

# Plot terminal node separation
autoplot(gs[[1]], cell_nodes, bins=100, strip.text='gate', axis_inverse_trans=F) + 
  ggcyto_par_set(facet=facet_wrap(~name, scales='free'), limits = list(x=c(0, 225), y=c(0, 225)))

# Plot type distribution
cell_count <- getTotal(gh, 'T Cells')
getPopStats(gs) %>% filter(Population %in% cell_nodes) %>% select(Population, Count) %>%
  mutate(Percent=Count/cell_count, Class=cell_class[Population]) %>% arrange(desc(Percent)) %>% 
  mutate(Population=fct_reorder(Population, Percent)) %>%
  ggplot(aes(x=Population, y=Percent, fill=Class)) +
  geom_bar(stat='identity') +
  scale_fill_brewer(palette='Set1') +
  scale_y_continuous(labels = scales::percent) +
  ggtitle('T Cell Type Distribution') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
# Choose cell type for each event

# Show assignment distribution
cell_nodes %>% map(~getIndices(gh, .)) %>% bind_cols %>% 
  apply(., 1, sum) %>% table

# For now, choose cell type arbitrarily when multiple are available
cell_types <- cell_nodes %>% map(~getIndices(gh, .)) %>% bind_cols %>% 
  apply(., 1, function(x) cell_nodes[x][1]) 
  #base::apply(., 1, function(x) paste(cell_nodes[x], collapse=',')) %>% table %>% sort

# Create data frame with measurments and cell type
fr <- getData(gh)
name_map <- parameters(fr)@data %>% filter(str_detect(name, 'Comp')) 
name_map <- set_names(name_map$name, name_map$desc)
df <- fr %>% exprs %>% as_tibble %>% 
  select(starts_with('Comp')) %>% 
  rename(!!!name_map) %>%
  mutate(type=cell_types)



