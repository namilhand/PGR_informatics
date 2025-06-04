library(tidyverse)

args <- commandArgs(trailingOnly=T)

dco_block.file <- args[1]
output <- args[2]

#dco_block.file <- "/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/smooth_dco_block/20201028_GBS_Col_dco.smooth.txt"

dco_block <- read_tsv(dco_block.file)

# trim libname
dco_block2 <- dco_block %>%
    mutate(lib = str_replace(lib, "results/06_tiger/lib", "")) %>%
    mutate(lib = str_replace(lib, "_MappedOn_tair10", "")) %>%
    mutate(lib = as.numeric(lib)) %>%
    arrange(lib, chr, start)

# dco_block to cotable
cotable <- NULL
for(i in unique(dco_block2$lib)){
    bylib <- filter(dco_block2, lib == i)

    for(j in 1:5){
        bylibchr <- filter(bylib, chr == j)

        if(nrow(bylibchr) == 0){
            next;
        }

       cotable.start <- bylibchr$end 
       cotable.start <- cotable.start[-length(cotable.start)]
       cotable.end <- bylibchr$start
       cotable.end <- cotable.end[-1]
       cotable.width <- cotable.end - cotable.start
       cotable.cos <- cotable.start + round(cotable.width/2)
       cotable.lib <- rep(i, length(cotable.cos))
       cotable.chr <- rep(j, length(cotable.cos))

       cotable.bylibchr <- cbind(lib=cotable.lib, chr=cotable.chr, start=cotable.start, end=cotable.end, cos=cotable.cos, width=cotable.width)

       #cotable.bylib <- rbind(cotable.bylib, cotable.bylibchr)
       cotable <- rbind(cotable, cotable.bylibchr)
    }
}

cotable.tb <- as_tibble(cotable)
write_tsv(cotable.tb, file=output)
