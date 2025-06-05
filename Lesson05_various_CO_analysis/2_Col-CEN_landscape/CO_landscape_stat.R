# Description 
library(tidyverse)
library(parallel)
library(GenomicRanges)

# GENOME INFO SETTING ====================================
#' 0. chromosome info 
wd <- "/home/namilhand/01_Projects/PGR_informatics/Lesson05_various_CO_analysis/2_Col-CEN_landscape"
col_cen <- read_tsv(file.path(wd, "data/colcen/Col-CEN_v1.2_chr15.fasta.fai"), col_names=c("chr", "length", "offset", "linelength", "linelength2")) %>%
    dplyr::select(c(chr, length))
chrs <- col_cen$chr
chr.ends <- col_cen$length

tha.cum <- c(0, cumsum(chr.ends))
tha.tot <- tha.cum[length(tha.cum)]

nrz_left <- c(14100001, 3100001, 13500001, 3000001, 11000001)
nrz_right <- c(17700001, 6400001, 16500001, 7100001, 14800001)

lrz_left <- c(13200001, 2400001, 11800001, 1500001, 10600001)
lrz_right <- c(18200001, 800001, 17600001, 8100001, 16200001)

cen_left <- c(14841110, 3823792, 13597188, 4203902, 11784131)
cen_right <- c(17559778, 6045243, 15733925, 6977949, 14551809)

nrz_mid <- (nrz_left + nrz_right)/2

col_cen.coord <- col_cen %>%
    add_column(nrz_left, nrz_right) %>%
    add_column(nrz_left.cum = nrz_left + tha.cum[1:5]) %>%
    add_column(nrz_right.cum = nrz_right + tha.cum[1:5]) %>%
    add_column(lrz_left.cum = lrz_left + tha.cum[1:5]) %>%
    add_column(lrz_right.cum = lrz_right + tha.cum[1:5]) %>%
    add_column(cen_left.cum = cen_left + tha.cum[1:5]) %>%
    add_column(cen_right.cum = cen_right + tha.cum[1:5])

#================ STATISTICS ===============

# make 1 Mb window stepping by 100 kb

makeWindow <- function(winsize, stepsize){
        wg.bin <- tibble(chrs=character(), bin.start=numeric(), bin.end=numeric(), cum.start=numeric(), cum.end=numeric())
        for(i in 1:5){
                coord.start <- seq(1, chr.ends[i], by=stepsize)
                coord.end <- coord.start + winsize -1
                coord.end <- c(coord.end[-length(coord.end)], chr.ends[i])
                cum.start <- coord.start+ tha.cum[i]
                cum.end <- coord.end + tha.cum[i]
                nchr <- rep(i, length(coord.start))
                bins <- tibble(chrs=paste0("Chr",nchr), bin.start=coord.start, bin.end=coord.end, cum.start=cum.start, cum.end=cum.end)
                wg.bin <- add_row(wg.bin, bins)
        }
        return(wg.bin)
}

window_1mb_100kb <- makeWindow(10^6, 100 * 10^3)

binCotable <- function(dat, genomeBin){
    # dat <- cotable.gbs.wt_male 
    # genomeBin <- window_1mb_100kb
    # genomeBin <- window_1mb_50kb
    winSize=genomeBin$bin.end[1] - genomeBin$bin.start[1] + 1

    total_co <- nrow(dat)
    co_by_chr <- group_by(dat, chrs) %>%
        summarise(coInChr = n())

    binned <- NULL
    collect.bin <- NULL
    for(i in 1:5){
        chr.dat <- filter(dat, chrs == paste0("Chr", i))
        chr.bin <- filter(genomeBin, chrs == paste0("Chr", i))
        for(j in 1:nrow(chr.bin)){
            dat.in.bin <- sum(chr.dat[["cos"]] >= chr.bin$bin.start[j] & chr.dat[["cos"]] <= chr.bin$bin.end[j])
            collect.bin <- c(collect.bin, dat.in.bin)
        }
    }
    coBin <- tibble(genomeBin, coInWindow=collect.bin, libSize=length(unique(dat$lib)))

    resultList <- mclapply(1:5, function(x){
        # x <- 1
        chrCoBin <- filter(coBin, chrs == paste0("Chr", x))
        if(with(chrCoBin[nrow(chrCoBin),], bin.end - bin.start) +1 < winSize){
            chrCoBin[nrow(chrCoBin),]$coInWindow <- chrCoBin[nrow(chrCoBin)-1,]$coInWindow
        }
        return(chrCoBin)
    }, mc.cores=5)

    result <- bind_rows(resultList) %>%
        left_join(co_by_chr, by = "chrs") %>%
        add_column(total_co=total_co)

    return(result)
}

# Read cotables
dir_cotable <- "data/crossovers/cotable"

## cotable
wt_i <- read_csv(file.path(dir_cotable, "Selected_inhouse_WT_Col-CEN_cotable.txt"))
wt_r <- read_csv(file.path(dir_cotable, "Selected_Rowan_WT_Col-CEN_cotable.txt")) %>%
    mutate(lib = lib + 240)
wt_cotable <- bind_rows(wt_i, wt_r)

mj3_1 <- read_csv(file.path(dir_cotable, "pJ3-mJ3_set1_Col-CEN_cotable.txt"))
mj3_2 <- read_csv(file.path(dir_cotable, "pJ3-mJ3_set2_Col-CEN_cotable.txt")) %>%
    mutate(lib = lib + 96)
mj3_3 <- read_csv(file.path(dir_cotable, "pJ3-mJ3_set3_Col-CEN_cotable.txt")) %>%
    mutate(lib = lib + 96 + 96)
mj3_cotable <- bind_rows(mj3_1, mj3_2, mj3_3)

h2aw_1 <- read_csv(file.path(dir_cotable, "h2aw-null_set1_Col-CEN_cotable.txt"))
h2aw_2 <- read_csv(file.path(dir_cotable, "h2aw-null_set2_Col-CEN_cotable.txt")) %>%
    mutate(lib = lib + 96)
h2aw_3 <- read_csv(file.path(dir_cotable, "h2aw-null_set3_Col-CEN_cotable.txt")) %>%
    mutate(lib = lib + 96 + 96)
h2aw_cotable <- bind_rows(h2aw_1, h2aw_2, h2aw_3)

cotable_list <- list(wt=wt_cotable, mj3=mj3_cotable, h2aw=h2aw_cotable)

cotable_list.binned_1mb100kb <- lapply(cotable_list, binCotable, window_1mb_100kb)

exacttest <- function(dat1, dat2, cutoff){
    dat1 <- dat1 %>%
        mutate(coOutWindow=total_co - coInWindow) %>%
        relocate(coOutWindow, .after = coInWindow) %>%
        dplyr::select(-libSize) %>%
        # dplyr::select(-total_co) %>%
        mutate(relCo=coInWindow/total_co)

    dat2 <- dat2 %>%
        mutate(coOutWindow=total_co - coInWindow) %>%
        relocate(coOutWindow, .after = coInWindow) %>%
        dplyr::select(-libSize) %>%
        # dplyr::select(-total_co) %>%
        mutate(relCo=coInWindow/total_co)

    dat12 <- left_join(dat1, dat2, by=c("chrs", "bin.start", "bin.end", "cum.start", "cum.end"))
    colnames(dat12) <- c("chrs", "bin.start", "bin.end", "cum.start", "cum.end", "dat1.inside", "dat1.outside",
                        "dat1.chr", "dat1.total", "dat1.relCo", "dat2.inside", "dat2.outside", "dat2.chr", "dat2.total", "dat2.relCo")
    
    dat12 <- mutate(dat12, diff.relCo = dat2.relCo - dat1.relCo)

    

    # xsq_bins <- sort(dat1$cum.start)
    
    pvalues.fisher <- c()
    # pvalues.prop <- c()
    # pvalues.barnard <- c()
    # pvalues.midp <- c()
    for(bin in 1:nrow(dat12)){
        xsq_dat <- matrix(c(dat12[bin,]$dat1.inside, dat12[bin,]$dat1.outside, dat12[bin,]$dat2.inside, dat12[bin,]$dat2.outside), nrow=2)

        # fisher_p <- fisher.test(xsq_dat, simulate.p.value=T, B=1e5)$p.value
        fisher_p <- fisher.test(xsq_dat)$p.value
        # prop_p <- prop.test(x=xsq_dat)$p.value
        # barnard_p <- exact.test(data=xsq_dat)
        # exact2x2_p <- exact2x2(xsq_dat, tsmethod="blaker")

        pvalues.fisher <- c(pvalues.fisher, fisher_p)
        # pvalues.prop <- c(pvalues.prop, prop_p)
    }
    pvalues.adj <- p.adjust(pvalues.fisher, method="BH") 
    result <- dplyr::select(dat12, chrs, cum.start, cum.end, dat1.inside, dat1.outside, dat1.chr, dat1.relCo, dat2.inside, dat2.outside, dat2.chr, dat2.relCo, diff.relCo) %>%
        add_column(pval=pvalues.fisher, fdr=pvalues.adj) %>%
        mutate(xsq.result=case_when(fdr < cutoff ~ "diff",
                                    TRUE ~ "ns"))

    return(result)
}

# exploring proper window and step size


# (which(exacttest(cotable_list.binned_1mb100kb[["gbs_wt_male"]], cotable_list.binned_1mb100kb[["comapper_wt_pollen"]], 0.05)$fdr < 0.05) * 0.1 * 10^6 - 0.1 * 10^6 + 2*10^6)/10^6
# # 2

# (which(exacttest(cotable_list.binned_1mb100kb[["gbs_recq4_male"]], cotable_list.binned_1mb100kb[["comapper_recq4_pollen"]], 0.05)$fdr < 0.05) * 0.1 * 10^6 - 0.1 * 10^6 + 2*10^6)/10^6
# # 5

# (which(exacttest(cotable_list.binned_1mb100kb[["gbs_wt"]], cotable_list.binned_1mb100kb[["comapper_wt_f2"]], 0.05)$fdr < 0.05) * 0.1 * 10^6 - 0.1 * 10^6 + 2*10^6)/10^6
# # 0

# (which(exacttest(cotable_list.binned_1mb100kb[["gbs_recq4"]], cotable_list.binned_1mb100kb[["comapper_recq4_f2"]], 0.05)$fdr < 0.05) * 0.1 * 10^6 - 0.1 * 10^6 + 2*10^6)/10^6
# # 6



######
# result:
# 1 mb window 100 kb step yields best result
######

# wt vs mj3
exacttest.wt_mj3 <- exacttest(cotable_list.binned_1mb100kb[["wt"]], cotable_list.binned_1mb100kb[["mj3"]], 0.05)
sum(exacttest.wt_mj3$xsq.result == "diff")
# result: 302

# wt vs h2aw
exacttest.wt_h2aw <- exacttest(cotable_list.binned_1mb100kb[["wt"]], cotable_list.binned_1mb100kb[["h2aw"]], 0.05)
filter(exacttest.wt_h2aw, xsq.result == "diff")
sum(exacttest.wt_h2aw$xsq.result == "diff")
# 35

dir.create(file.path("stat"), recursive=T)

write_csv(file="stat/wt-vs-mj3_landscape.csv", exacttest.wt_mj3)
write_csv(file="stat/wt-vs-h2aw_landscape.csv", exacttest.wt_h2aw)

# Combining overlapping regions

reduceRegions <- function(dat, fdr_cutoff){
    dat <- filter(dat, fdr < fdr_cutoff)
    dat.gr <- makeGRangesFromDataFrame(dat, seqnames.field = "chrs", start.field="cum.start", end.field="cum.end")
    dat.gr_reduced <- reduce(dat.gr)
    dat.reduced <- as_tibble(dat.gr_reduced) %>% dplyr::select(-strand)
    colnames(dat.reduced) <- c("chrs", "cum.start", "cum.end", "width")
    return(dat.reduced)
}

exacttest.wt_mj3.reduced.pval05 <- reduceRegions(exacttest.wt_mj3, 0.05)
sum(exacttest.wt_mj3.reduced.pval05$width)
# 61,800,000
sum(exacttest.wt_mj3.reduced.pval05$width)/tha.tot
# 0.47

exacttest.wt_h2aw.reduced.pval05 <- reduceRegions(exacttest.wt_h2aw, 0.05)
sum(exacttest.wt_h2aw.reduced.pval05$width)
# 8,000,000
sum(exacttest.wt_h2aw.reduced.pval05$width)/tha.tot
# 0.06