library(tidyverse)
library(foreach)
library(doParallel)

#====== for parallel processing ======
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

#====== input ======
dirtsv <- "/home/namilhand/01_Projects/PGR_informatics/Lesson02_FilteringSNPs/result/04_filtering_variants/variant_tables"

tsv_files <- list.files(dirtsv, pattern=".tsv")
tsv_list <- lapply(file.path(dirtsv, tsv_files), read_tsv)
## total 334795 variants

markers <- bind_rows(tsv_list) %>% arrange(CHROM, POS) %>%
    filter(NCALLED != 0)
## variants with NCALLED = 0 result in errors during fisher's exact test

colnames(markers) <- c("chr", "pos", "ref", "alt", "type", "dp", "af", "nsamples", "homref", "het", "homvar", "ncalled", "nocall", "var")

# NOTES: how are these markers obtained
# 1. allele frequency higher than 0.1
# 2. only snp

#====== output ======
dirout="result/05_final_snp"
dir.create(dirout, recursive=T)

#====== mendelian segregation test =======

mendelTest <- function(dat){
    # dat <- markers[1,]
    observed <- c(RR=dat$homref, RA=dat$het, AA=dat$homvar)
    expected <- c(RR=1/4, RA=1/2, AA=1/4) * sum(observed)
    chisq.table <- rbind(observed, expected)
    p <- fisher.test(chisq.table)$p.value

    # why fisher's exact test not chi-square test?
    return(p)
}

# mendel_test_result <- foreach(i=1:nrow(markers), .combine=c) %dopar%{
#     mendelTest(markers[i,])
# }

mendel_test_result <- foreach(i=1:nrow(markers), .combine=c) %dopar%{
    mendelTest(markers[i,])
}
markers$pval <- mendel_test_result

#======= exploratory analysis on SNPs ======
# 334,795 variants

my_theme <- theme_classic() +
    theme(text = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.text = element_text(size=6, family="Arial", colour="black")) +
    theme(plot.title = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.title = element_text(size=6, family="Arial", colour="black")) +
    theme(legend.margin = margin(0.5,0.5,0.5,0.5)) +
    theme(strip.text = element_text(size=6, family="Arial", colour="black")) +
    theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))
theme_set(my_theme)

# 1. NCALLED

mean.ncalled <- mean(markers$ncalled) # 39.05

p.hist_ncalled <- ggplot() +
    geom_histogram(data=markers, aes(x=ncalled), binwidth=1) +
    geom_vline(xintercept=20, colour="red", linewidth=0.4)

png(file=file.path(dirout, "hist_ncalled.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.hist_ncalled)
dev.off()

## ncalled > 10 : 332601
## ncalled > 15 : 326175
## ncalled > 20 : 311,608

# 2. allele frequency

mean.af <- mean(markers$af) # 0.47

p.hist_af <- ggplot() +
    geom_histogram(data=markers, aes(x=af), binwidth=0.01) +
    geom_vline(xintercept=0.4, colour="red", linewidth=0.4) +
    geom_vline(xintercept=0.6, colour="red", linewidth=0.4)

png(file=file.path(dirout, "hist_allele_frequency.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.hist_af)
dev.off()

## AF 0.4 - 0.6 : 235,004 variants

# 3. fisher's exact test p-value

p.pval <- ggplot() +
    geom_histogram(data=markers, aes(x=pval), binwidth=0.01) +
    geom_vline(xintercept=0.1, colour="red", linewidth=0.4)

png(file=file.path(dirout, "hist_exact-test.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.pval)
dev.off()

## variants with pvalue < 0.05 : 34,036
## variants with pvalue < 0.10 : 52,228

# 4. pvalue by allele frequency

p.pval_af <- ggplot() +
    geom_point(data=markers, aes(x=af, y=pval), shape=21, colour="NA", fill="black", alpha=0.3, size=0.1) +
    geom_vline(xintercept=c(0.4, 0.6), colour="red", linewidth=0.4) +
    geom_hline(yintercept=0.1, colour="red", linewidth=0.4) +
    scale_y_continuous(limits=c(0, 1)) +
    scale_x_continuous(limits=c(0, 1))

png(file=file.path(dirout, "plot_pvalue-by-af.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.pval_af)
dev.off()

#====== Final SNP sets ======
# 1. NCALLED > 20
# 2. AF 0.4 - 0.6
# 3. chi-square p-value > 0.1

markers.fin <- filter(markers, ncalled > 20, af > 0.4 & af < 0.6, pval > 0.1)
# 223,870 snps


markers.fin.tsv <- dplyr::select(markers.fin, c(chr, pos, ref, alt))

write_tsv(markers.fin.tsv, file="result/05_final_snp/coller_gbs_markers.tsv", col_names=F)