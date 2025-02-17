library(tidyverse)

allvar="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/04_result/hcr7_allvar.se.tsv"
dirout="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/04_result"
chromosome="Chr3"
upperb=17000000
lowerb=25000000

prefix <- basename(allvar) %>% str_replace(".se.tsv", "")

cand <- read_tsv(allvar, col_names=T) %>%
    filter(chr == chromosome) %>%
    filter(pos > upperb, pos < lowerb) %>%
    filter(ratio > 0.35)

write_tsv(cand, file.path(dirout, paste0(prefix, ".candidate.tsv")))
