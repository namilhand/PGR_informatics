library(tidyverse)
# h2aw6/7 show DCO distance comparable to WT

# dco cotable
dir_cotable="/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/dco_cotable/combine_by_genotype"
col_dco.file="col_dco.cotable.txt"
hta6_dco.file="hta6_dco.cotable.txt"
hta7_dco.file="hta7_dco.cotable.txt"
hta12_dco.file="hta12_dco.cotable.txt"
hta67_dco.file="hta67_dco.cotable.txt"
h2aw_dco.file="h2aw_dco.cotable.txt"
#mmH_dco.cotable.txt
#mmS_dco.cotable.txt
#mmHS_dco.cotable.txt

col_dco <- read_tsv(file.path(dir_cotable, col_dco.file))
hta6_dco <- read_tsv(file.path(dir_cotable, hta6_dco.file))
hta67_dco <- read_tsv(file.path(dir_cotable, hta67_dco.file))

calcDistance <- function(dat){
    dat <- col_dco
    dat <- arrange(dat, lib, chr, cos)
    dco_start <- filter(dat, row_number()%%2 ==1)
    dco_end <- filter(dat, row_number()%%2 ==0)
}
