test
================

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
