library(tidyverse)

# read smooth genotype fragments
args <- commandArgs(trailingOnly=T)
dir_tiger <- args[1]
libname <- args[2]
dirout <- args[3]
#dir_tiger="/datasets/data_4/nison/GBS/GBS_marker_v2/20220602_mmH2AW_set2/results/06_tiger"

smooth_list.file <- list.files(dir_tiger, pattern="smooth", full.names=T)
smooth_list.df <- lapply(smooth_list.file, read_tsv, col_names=F)

# find DCO block in a library
##  input = a tiger output (smooth)
findDCO <- function(dat){
#    dat <- smooth_list.df[[85]]
    colnames(dat) <- c("lib", "chr", "start", "end", "gt")

    # correct genotype error of tiger pipeline
    dat2 <- dat %>%
        mutate(gt = case_when(chr != 5 & gt == "LL" ~ "CL",
                              chr != 5 & gt == "CL" ~ "LL",
                              TRUE ~ gt))

    homhethom.list <- list()

    # find HOM-HET-HOM genotype blocks in each chromosome
    for(i in 1:5){
#        i <- 5
        dat2.chr <- filter(dat2, chr == i)

        # If the number of block is smaller than 3, pass the chromosome
        if(nrow(dat2.chr) < 3){
            next;
        }

        # If the number of block is greater than 3,
        # find the start position of HOM-HET-HOM genotype block
        dat2.chr$gt2 <- c(dat2.chr$gt[2:length(dat2.chr$gt)], "empty")
        dat2.chr$gt3 <- c(dat2.chr$gt[3:length(dat2.chr$gt)], "empty", "empty")

        dat2.chr.find_dco <- mutate(dat2.chr, dco = case_when(gt == "LL" & gt2 == "CL" & gt3 == "LL" ~ "dco",
                                                     gt == "CC" & gt2 == "CL" & gt3 == "CC" ~ "dco",
                                                     TRUE ~ "not_dco"))
        dco_index <- which(dat2.chr.find_dco$dco == "dco")

        # If there is no DCO, pass the chromosome
        if(length(dco_index) == 0){
            next;
        }

        dco_index_series <- sort(unique(c(dco_index, dco_index + 1, dco_index + 2)))
        # If there is a DCO, retain only HOM-HET-HOM genotype block from the chromosome
        dat2.chr.homhethom <- dat2.chr[dco_index_series,]

        # HOM-HET-HOM block of each chromosome is stored in the list
        homhethom.list[[i]] <- dat2.chr.homhethom
    }

    # now all # HOM-HET-HOM block of an F2 is saved in list
    # combine them and return
    homhethom <- bind_rows(homhethom.list)

    return(homhethom)
}

# apply findDCO to all smooth genotype fragment files
all_smooth.dco.list <- lapply(smooth_list.df, findDCO)


all_smooth.dco <- bind_rows(all_smooth.dco.list) %>%
    dplyr::select(1:5)

write_tsv(all_smooth.dco, file=file.path(dirout, paste0(libname, "_dco.smooth.txt")))
