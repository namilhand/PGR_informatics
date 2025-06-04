library(tidyverse)

dirout <- "/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/v2/01_dco_table"
libname <- "recq4"
dir_tiger <- "/datasets/data_4/nison/GBS/GBS_marker_v2/20210322_recq4/results/06_tiger"
smooth_list.file <- list.files(dir_tiger, pattern="smooth", full.names=T)
smooth_list.df <- lapply(smooth_list.file, read_tsv, col_names=F)

# find DCO block in a library and return double crossover table for a library
##  input = a tiger output (smooth)
findDCO <- function(dat){
    #dat <- filter(smooth_list.df[[12]])
    colnames(dat) <- c("lib", "chr", "start", "end", "gt")

    no_lib <- str_replace(dat$lib[1], "results/06_tiger/lib", "") %>%
                str_replace("_MappedOn_tair10", "") %>%
                as.numeric()

    # correct genotype error of tiger pipeline
    dat2 <- dat %>%
        mutate(gt = case_when(chr != 5 & gt == "LL" ~ "CL",
                              chr != 5 & gt == "CL" ~ "LL",
                              TRUE ~ gt))

    lib_dco <- tibble()

    # find HOM-HET-HOM genotype blocks in each chromosome
    for(i in 1:5){
#        i <- 1
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

        # double crossover table for each chromsome
        libchr_dco <- tibble()
        for(i_dco in dco_index){
            #i_dco <- 2
            dco.start <- dat2.chr.find_dco$end[i_dco:(i_dco+1)]
            dco.end <- dat2.chr.find_dco$start[(i_dco +1):(i_dco+2)]
            dco_temp <- tibble(lib=no_lib, chr=i, start=dco.start, end=dco.end, cos = (dco.start + dco.end)/2, width = dco.end - dco.start)
            libchr_dco <- bind_rows(libchr_dco, dco_temp)
        }
        lib_dco <- bind_rows(lib_dco, libchr_dco) 
    }

    # now all # HOM-HET-HOM block of an F2 is saved in list
    # combine them and return
    # homhethom <- bind_rows(homhethom.list)

    return(lib_dco)
}

# apply findDCO to all smooth genotype fragment files
all.dco_table.list <- lapply(smooth_list.df, findDCO)


all.dco_table <- bind_rows(all.dco_table.list)
    dplyr::select(1:5)

write_tsv(all.dco_table, file=file.path(dirout, paste0(libname, "_dco.smooth.txt")))
