
OMIP 030
========

``` r
library(magrittr)
library(tidyverse)
library(fs)
library(janitor)
library(FlowRepositoryR)
library(flowCore)
library(openCyto)
library(ggcyto)
library(gplots)
library(Rtsne)

map <- purrr::map
add <- flowCore::add
filter <- dplyr::filter
```

Load a flow workspace, which will contain all gaiting and cell expression information:

``` r
ws <- openWorkspace('~/analysis/data/omip/omip-030/clean/OMIP_PBMC_clean.wsp')
gs <- parseWorkspace(ws, name='TCells', 
                     path='/home/eczech/analysis/data/omip/omip-030/clean', isNcdf=FALSE)
# Extract single GatingHierarchy from GatingSet
gh <- gs[[1]]
```

``` r
print(getData(gh))
```

    ## flowFrame object 'Set 1_B1-4_cct.fcs_1664607'
    ## with 1664607 cells and 22 observables:
    ##                             name      desc       range minRange
    ## $P1                        FSC_A      <NA> 262144.0000  0.00000
    ## $P2                        FSC_W      <NA> 262144.0000  0.00000
    ## $P3                        SSC_A      <NA> 262144.0000  0.00000
    ## $P4                        SSC_W      <NA> 262144.0000  0.00000
    ## $P5        Comp-LIVE_DEAD_BLUE_A       L+D    206.2050 50.70977
    ## $P6  Comp-BRILLIANT_VIOLET_421_A      CD38    206.2049 50.70982
    ## $P7                  Comp-V500_A       CD4    206.2050 50.70978
    ## $P8                Comp-BV_570_A    CD45RA    206.2050 50.70978
    ## $P9              Comp-QDOT_605_A       CD3    206.2050 50.70976
    ## $P10             Comp-QDOT_655_A       CD8    206.2049 50.70979
    ## $P11               Comp-BV_711_A      CD25    206.2049 50.70984
    ## $P12               Comp-BV_785_A     CD196    206.2050 50.70976
    ## $P13                 Comp-FITC_A     CD183    206.2050 50.70978
    ## $P14   Comp-PER_CP_E_FLUOR_710_A     CD161    206.2050 50.70976
    ## $P15                   Comp-PE_A     CD194    206.2050 50.70977
    ## $P16         Comp-PE_TEXAS_RED_A     CD197    206.2050 50.70977
    ## $P17             Comp-PE_CY5_5_A      CD20    206.2050 50.70976
    ## $P18               Comp-PE_CY7_A CD185-bio    206.2050 50.70976
    ## $P19                  Comp-APC_A     CCR10    206.2050 50.70978
    ## $P20      Comp-ALEXA_FLUOR_700_A      Ki67    206.2050 50.70977
    ## $P21              Comp-APC_CY7_A     CD127    206.2049 50.70979
    ## $P22                        TIME      <NA>  32768.0000  0.00000
    ##         maxRange
    ## $P1  262144.0000
    ## $P2  262144.0000
    ## $P3  262144.0000
    ## $P4  262144.0000
    ## $P5     256.9147
    ## $P6     256.9147
    ## $P7     256.9147
    ## $P8     256.9147
    ## $P9     256.9147
    ## $P10    256.9147
    ## $P11    256.9147
    ## $P12    256.9147
    ## $P13    256.9147
    ## $P14    256.9147
    ## $P15    256.9147
    ## $P16    256.9147
    ## $P17    256.9147
    ## $P18    256.9147
    ## $P19    256.9147
    ## $P20    256.9147
    ## $P21    256.9147
    ## $P22  32768.0000
    ## 298 keywords are stored in the 'description' slot

Plot the gating workflow starting from where T Cells are first identified:

``` r
plot(gs, 'T Cells', width=6, height=1, fontsize=18, shape='rectangle')
```

![](omip030_files/figure-markdown_github/show_workflow-1.png)

Identify terminal nodes in the workflow which were, by convention, named with either `CD4+` or `CD8+` prefixed to their population names:

``` r
all_nodes <- getNodes(gs, path=1)
cell_nodes <- getNodes(gs, path=1) %>% str_match(., 'CD[4|8]\\+ .*') %>% discard(is.na)
cell_type_to_class <- function(x) str_detect(x, 'CD4+') %>% ifelse('CD4+', 'CD8+') 
cell_class <- cell_type_to_class(cell_nodes) %>% set_names(cell_nodes)
cell_nodes_th2 <- c('CD4+ Th2g1', 'CD4+ Th2g2')
cell_nodes
```

Show the gating for terminal nodes, which correspond to specific cell types:

``` r
autoplot(gs[[1]], cell_nodes, bins=100, strip.text='gate', axis_inverse_trans=FALSE) + 
  ggcyto_par_set(facet=facet_wrap(~name, scales='free'), limits = list(x=c(0, 225), y=c(0, 225)))
```

![](omip030_files/figure-markdown_github/leaf_gates-1.png)

Plot the distribution of different cell types noting that many cells are assigned to more than one cell type, and here any duplicate assignments are resolved by counting *all* cells assigned to any one type:

``` r
cell_count <- getTotal(gh, 'T Cells')

cell_pop <- getPopStats(gs) %>% 
  filter(Population %in% cell_nodes) %>% 
  select(Population, Count) %>% 
  arrange(Count)

cell_pop %>%
  mutate(Percent=Count/cell_count, Class=cell_type_to_class(Population)) %>% 
  arrange(desc(Percent)) %>% 
  mutate(Population=fct_reorder(Population, Percent)) %>%
  ggplot(aes(x=Population, y=Percent, fill=Class)) +
  geom_bar(stat='identity') +
  scale_fill_brewer(palette='Set1') +
  scale_y_continuous(labels = scales::percent) +
  ggtitle('T Cell Type Distribution') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
```

![](omip030_files/figure-markdown_github/cell_type_distribution-1.png)

To move forward though, determining a single cell type for each cell is done by using the label with the lowest frequency to maximize the size of populations relating to more rare cells. Before doing that, here is the frequency of all cell type assignments where multiple assignments are comma separated:

``` r
cell_nodes %>% map(~getIndices(gh, .)) %>% bind_cols %>% 
  apply(., 1, function(x) paste(cell_nodes[x], collapse=', ')) %>% 
  tabyl(var1='assignment') %>%
  arrange(n) %>% adorn_pct_formatting() %>%
  knitr::kable()
```

| .                                            |       n| percent |
|:---------------------------------------------|-------:|:--------|
| CD4+ CM, CD4+ FH Th17-like, CD4+ Naive Treg  |       1| 0.0%    |
| CD4+ CM, CD4+ Memory Treg                    |       1| 0.0%    |
| CD4+ CM, CD4+ Th17, CD4+ Naive Treg          |       1| 0.0%    |
| CD4+ CM, CD4+ Th2g2, CD4+ Naive Treg         |       1| 0.0%    |
| CD4+ EM, CD4+ Naive Treg                     |       1| 0.0%    |
| CD4+ EM, CD4+ Th22, CD4+ Naive Treg          |       1| 0.0%    |
| CD4+ EM, CD4+ Th2g2, CD4+ Memory Treg        |       1| 0.0%    |
| CD4+ FH Th17-like                            |       1| 0.0%    |
| CD4+ Th17                                    |       1| 0.0%    |
| CD4+ Th2g1                                   |       1| 0.0%    |
| CD4+ CM, CD4+ FH Th17-like, CD4+ Memory Treg |       2| 0.0%    |
| CD4+ Th22                                    |       3| 0.0%    |
| CD4+ Th2g2                                   |       3| 0.0%    |
| CD4+ EM, CD4+ Th2g1, CD4+ Memory Treg        |       4| 0.0%    |
| CD4+ FH Th2-like                             |       4| 0.0%    |
| CD4+ Th9                                     |       4| 0.0%    |
| CD4+ CM, CD4+ Th9, CD4+ Memory Treg          |       5| 0.0%    |
| CD4+ EM, CD4+ Th17, CD4+ Memory Treg         |       5| 0.0%    |
| CD4+ CM, CD4+ ThG , CD4+ Memory Treg         |       6| 0.0%    |
| CD4+ Effector, CD4+ Naive Treg               |       6| 0.0%    |
| CD4+ Naive, CD4+ Naive Treg                  |       7| 0.0%    |
| CD4+ CM, CD4+ Th22, CD4+ Memory Treg         |      12| 0.0%    |
| CD4+ CM, CD4+ Th2g1, CD4+ Memory Treg        |      12| 0.0%    |
| CD4+ EM, CD4+ Th9, CD4+ Memory Treg          |      14| 0.0%    |
| CD4+ CM, CD4+ Th17, CD4+ Memory Treg         |      15| 0.0%    |
| CD8+ mNKT-MAIT, CD8+ Effector                |      18| 0.0%    |
| CD4+ EM, CD4+ Th22, CD4+ Memory Treg         |      25| 0.0%    |
| CD8+ mNKT-MAIT, CD8+ CM                      |     653| 0.0%    |
| CD8+ mNKT-MAIT, CD8+ EM                      |     826| 0.0%    |
| CD4+ EM, CD4+ FH Th1-like                    |    1335| 0.1%    |
| CD4+ EM, CD4+ FH Th17-like                   |    1693| 0.1%    |
| CD4+ EM, CD4+ FH Th2-like                    |    2453| 0.1%    |
| CD4+ EM, CD4+ ThG                            |    2628| 0.2%    |
| CD4+ EM, CD4+ Th2g1                          |    5223| 0.3%    |
| CD4+ EM, CD4+ Th2g2                          |    5256| 0.3%    |
| CD4+ EM, CD4+ Th17                           |    6185| 0.4%    |
| CD4+ CM, CD4+ ThG                            |    6873| 0.4%    |
| CD4+ CM, CD4+ FH Th1-like                    |    7019| 0.4%    |
| CD4+ CM, CD4+ Th22                           |    7381| 0.4%    |
| CD4+ EM, CD4+ Th22                           |    7564| 0.5%    |
| CD8+ mNKT-MAIT                               |   11612| 0.7%    |
| CD4+ CM, CD4+ FH Th2-like                    |   11868| 0.7%    |
| CD4+ Effector                                |   12677| 0.8%    |
| CD4+ CM, CD4+ Th17                           |   15732| 0.9%    |
| CD4+ CM, CD4+ FH Th17-like                   |   15891| 1.0%    |
| CD4+ Naive Treg                              |   18770| 1.1%    |
| CD4+ EM                                      |   18937| 1.1%    |
| CD8+ CM                                      |   19226| 1.2%    |
| CD4+ Memory Treg                             |   19398| 1.2%    |
| CD4+ EM, CD4+ Th9                            |   21043| 1.3%    |
| CD4+ CM, CD4+ Th2g1                          |   22107| 1.3%    |
| CD4+ CM, CD4+ Th9                            |   23019| 1.4%    |
| CD4+ CM, CD4+ Th2g2                          |   23151| 1.4%    |
| CD4+ CM                                      |   23465| 1.4%    |
| CD8+ Effector                                |   39677| 2.4%    |
| CD8+ EM                                      |   43847| 2.6%    |
| CD8+ Naive                                   |  119979| 7.2%    |
| CD4+ Naive                                   |  191831| 11.5%   |
|                                              |  957133| 57.5%   |

Choose one cell type for each cell:

``` r
# For now, favor less frequent cell types when assigning single type
get_cell_type <- function(x){
  ctyps <- cell_nodes[x]
  if (is_empty(ctyps)) NA
  # Select rarest cell type if there are multiple
  else cell_pop$Population[which.max(cell_pop$Population %in% ctyps)]
}
cell_types <- cell_nodes %>% map(~getIndices(gh, .)) %>% bind_cols %>% 
  apply(., 1, get_cell_type) 
cell_types <- case_when(cell_types %in% cell_nodes_th2 ~ 'CD4+ Th2', TRUE ~ cell_types)

# Create data frame with measurments and cell type
fr <- getData(gh)
name_map <- parameters(fr)@data %>% filter(str_detect(name, 'Comp'))
name_map <- set_names(name_map$name, name_map$desc)
df <- fr %>% exprs %>% as_tibble %>% 
  select(starts_with('Comp')) %>% 
  rename(!!name_map) %>%
  mutate(type=cell_types) %>%
  filter(!is.na(type)) %>%
  mutate(class=cell_type_to_class(type))
```

Calculate the Kolmogrov Smirnov statistic for each cell type and marker, where the reference sample for any one combination is all *other* T Cells (for the same marker):

``` r
dfg <- df %>% group_by(type) %>% do({
  d <- .
  #dp <- df %>% filter(type != d$type[1] & class == d$class[1]) %>% select_if(is.numeric)
  dp <- df %>% filter(type != d$type[1]) %>% select_if(is.numeric)
  dt <- d %>% select_if(is.numeric)
  assertthat::are_equal(colnames(dp), colnames(dt))
  stats <- lapply(colnames(dp), function(col) {
    x <- pull(dt, col)
    y <- pull(dp, col)
    stat <- ks.test(x, y)$statistic[[1]]
    sign(median(x) - median(y)) * stat
  })
  names(stats) <- colnames(dp)
  data.frame(stats)
}) %>% ungroup

dfm <- dfg %>% select(-type) %>% as.matrix
row.names(dfm) <- dfg$type
heatmap.2(
  t(dfm), Rowv=TRUE, Colv=TRUE, 
  col = colorRampPalette(RColorBrewer::brewer.pal(11, 'RdBu'))(16), cexCol=.8, cexRow=.8, 
  margins = c(8, 8), 
  dendrogram = 'col', trace='none', scale='none',
  colsep=1:nrow(dfm), rowsep=1:ncol(dfm),
  density.info='none', sepwidth=c(.05,.05), key.title='KS', key.xlab='', 
  keysize=.5, key.par =list(cex=.3)
)
```

![](omip030_files/figure-markdown_github/ks-1.png)

Plot TSNE decomposition:

``` r
set.seed(1)
df_samp <- df %>% sample_n(10000)
tsne <- df_samp %>% select(-type, -class) %>% as.matrix %>% Rtsne
df_tsne <- tsne$Y %>% as_tibble %>%
  mutate(type=df_samp$type, class=df_samp$class)

df_tsne %>% mutate(type=factor(type, levels=sample(unique(type)))) %>%
  ggplot(aes(x=V1, y=V2, color=type, shape=class)) + geom_point() +
  flow_theme
```

![](omip030_files/figure-markdown_github/tsne-1.png)
