
library(snakecase)
library(tidyverse)
library(flowDensity)

clean_flow_frame <- function(fr){
  spill_fields <- fr %>% spillover %>% discard(is_null) %>% names
  if (!is_empty(spill_fields))
    for (field in spill_fields)
      colnames(fr@description[[field]]) %<>% clean_names
  colnames(fr) %<>% clean_names
  fr
}

clean_names <- function(names, case='all_caps') {
  # Lift from janitor::clean_names
  names %>% gsub("'", "", .) %>% gsub("\"", "", .) %>% 
    gsub("%", ".percent_", .) %>% gsub("#", ".number_", .) %>% 
    gsub("^[[:space:][:punct:]]+", "", .) %>% make.names(.) %>% 
    snakecase::to_any_case(case = case, sep_in = "\\.", transliterations = c("Latin-ASCII"), parsing_option = 4)
}

get_gate_coord_string <- function(gate, digits=4){
  # Convert gate boundary matrix (which is Nx2) to list of column vectors
  gate@boundaries %>% round(digits=digits) %>% as.data.frame %>% c %>% set_names(c('x', 'y')) %>%
  
    # Convert each column vector to csv string enclosed by vector constructor
    map(~paste(., collapse=',')) %>%
    map(~str_glue('c({.})')) %>%
      
    # Add key/value declaration
    imap(~paste(.y, .x, sep='=')) %>%
      
    # Collapse to single line like "x=c(1,2,3),y=c(4,5,6)"
    paste(collapse=',')
}


draw_gate <- function(fr, channels){
  flowDensity::plotDens(fr, channels)
  coords <- locator(n = 4, type = "p", lwd = 2, pch = 16, col = "red")
  gate <- coords %>% data.frame %>% set_names(channels)
  polygonGate(.gate=gate)
}
