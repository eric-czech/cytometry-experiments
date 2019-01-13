##### EDA

df <- read_delim('~/Downloads/CD8_HC_001.fcs_raw_events.txt', delim='\t', skip=1)
#df <- read_delim('~/Downloads/CD8_HC_003.fcs_raw_events.txt', delim='\t', skip=1)
#df <- read_delim('~/Downloads/CD8_HIV_005.fcs_raw_events.txt', delim='\t', skip=1)
dft <- df %>% mutate_at(vars(-Time, -Event_length), asinh) 
library(ggplot2)
library(dplyr)



# Density 2D gist
devtools::source_gist('bf63a4833cf991a1595e3fc503856c4f')

biplot <- function(df, x, y){
  df %>%
    mutate(density=get_density(!!x, !!y, nbins=64)) %>%
    #dplyr::sample_n(100000) %>%
    ggplot(aes(x=!!x, y=!!y, color=density)) + geom_point(size=1, alpha=.5) +
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    xlim(0, 5) +
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
}

biplot(dft, as.name('CD45RA'), as.name('CCR7'))
biplot(dft, as.name('CD45RO'), as.name('CCR7'))


dff <- dft %>% dplyr::filter(CD45RO < .5) 

# dff %>% dplyr::select(starts_with('CD'), CCR7) %>%
#   sample_n(1000) %>%
#   melt(id.vars='CCR7') %>% 
#   group_by(variable) %>% dplyr::mutate(density=get_density(CCR7, value)) %>% ungroup %>%
#   mutate(is_naive=CCR7 > 1) %>%
#   ggplot(aes(x=CCR7, y=value, color=density)) +
#   geom_point(size=1, alpha=.5) +
#   scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
#   facet_wrap(~variable)

#dff %>% dplyr::select(starts_with('CD'), CCR7) %>%
dft %>% dplyr::filter(CD45RO < .5)  %>% 
  dplyr::select(CCR7, Tbet, CD38, TIGIT, Eomes, TOX, Helios, `2B4`) %>%
  melt(id.vars='CCR7') %>% 
  mutate(is_naive=CCR7 > 1) %>%
  ggplot(aes(x=is_naive, y=value)) +
  geom_violin() + 
  facet_wrap(~variable, scales='free')


## Batch Analysis

dfa <- bind_rows(lapply(c(1, 2, 3, 4, 6, 7), function(i){
  filename <- sprintf('CD8_HC_00%s.fcs_raw_events.txt', i)
  path <- file.path('~/Downloads', filename)
  df <- read_delim(path, delim='\t', skip=1)  
  df <- df %>% mutate_at(vars(-Time, -Event_length), asinh) 
  df$Sample <- filename
  if (nrow(df) > 50000)
    df %>% sample_n(50000)
  else
    df
}))


x <- as.name('CD45RO')
y <- as.name('CCR7')
# x <- as.name('CD4')
# y <- as.name('CD8')
dfa %>%
  group_by(Sample) %>% mutate(density=get_density(!!x, !!y, nbins=64)) %>% ungroup %>%
  ggplot(aes(x=!!x, y=!!y, color=density)) + geom_point(size=1, alpha=.5) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
  facet_wrap(~Sample) +
  xlim(0, 5) + ylim(0, 7) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

# CD244 = 2B4
dfa %>% dplyr::filter(CD45RO < .5)  %>% 
  dplyr::select(Sample, CCR7, Tbet, CD38, TIGIT, Eomes, TOX, Helios, `2B4`, `PD-1`) %>%
  melt(id.vars=c('Sample', 'CCR7')) %>% 
  mutate(is_naive=CCR7 > 1) %>%
  ggplot(aes(x=is_naive, y=value, fill=Sample)) +
  geom_violin(scale='width', draw_quantiles=c(0.25, .5, 0.75)) + 
  scale_fill_brewer(palette = 'Set1') +
  ggtitle('CD8+, CD45RO- Healthy Human T Cell Expression Distributions\n*Lines denote 25, 50, and 75th percentiles') +
  labs(x='CCR7+ (i.e. naive)') +
  facet_wrap(~variable, scales='free', ncol = 2)

  