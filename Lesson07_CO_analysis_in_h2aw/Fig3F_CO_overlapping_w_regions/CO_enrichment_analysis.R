library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(ggsci)
library(rtracklayer)
library(extrafont)
font_import(pattern="Arial", prompt=F)
loadfonts(device="pdf")

#================================
# cotables
#================================

# read cotables
dir_cotable="/home/namilhand/01_Projects/H2AW/h2aw_server/kb_scale_co/data/cotables"
co.wt_ir <- read_csv(file.path(dir_cotable, "wt_all_cov1.5_cotable.txt"), col_names=T)
co.h2aw=read_csv(file.path(dir_cotable, "h2aw_all_cotable.txt"), col_names=T) %>%
    mutate(lib = case_when(genotype == "mmH" ~ lib,
                            genotype == "hta6" ~ lib + 192,
                            genotype == "hta7" ~ lib + 192 + 192,
                            genotype == "hta67" ~ lib + 192 + 192 + 192,
                            genotype == "h2aw" ~ lib + 192 + 192 + 192 + 192))


co.wt.gr <- makeGRangesFromDataFrame(co.wt_ir, seqnames.field="chrs", start.field="start", end.field="stop", starts.in.df.are.0based=F, keep.extra.columns=T)
co.h2aw.gr <- makeGRangesFromDataFrame(co.h2aw, seqnames.field="chrs", start.field="start", end.field="stop", starts.in.df.are.0based=F, keep.extra.columns=T)
### Note: h2aw includes 2 F2 library with no CO site: 3, 136
### h2aw null lib #3 --> 3 + 192*4 = 771;
### h2aw null lib #136 --> 136 + 192*4 = 904;

#================================
# various genomic features to count overlap with CO
#================================
# readBed: read BED file and convert into GRanges
readBed <- function(dat, start_open=T){
    dat <- read_tsv(dat, col_names=F)
    
    colnames(dat) <- c("seqnames", "start", "end", "name", "score", "strand")
    if(start_open==T){
        result <- makeGRangesFromDataFrame(dat, starts.in.df.are.0based=T)
    } else {
        result <- makeGRangesFromDataFrame(dat, starts.in.df.are.0based=F)
    }
    #result <- dat
    return(result)
}


dir_features="/home/namilhand/01_Projects/H2AW/h2aw_server/kb_scale_co/data/features"

peri.gene <- readBed(dat=file.path(dir_features, "Araport11_genes_pericentromere.bed"))
# n=4552
peri.te <- readBed(dat=file.path(dir_features, "Araport11_TEs_pericentromere.bed"))
# n=17756
peri <- readBed(dat=file.path("/home/namilhand/01_Projects/H2AW/h2aw_server/data/tair10_pericen_coord/ath_pericen_coord.bed")) 

# short te ( < 1000; n=13436)
peri.teshort <- peri.te[width(peri.te) < 1000]
# long te ( > 2000; n=2188)
peri.telong <- peri.te[width(peri.te) > 2000]

# upstream/downstream 2kb of genes and TEs
peri.gene.upstr <- promoters(peri.gene, upstream=2000, downstream=0)
peri.gene.downstr <- peri.gene %>% flank_downstream(width=2000)

# peri.te.upstr <- promoters(peri.te, upstream=2000, downstream=0)
# peri.te.downstr <- peri.te %>% flank_downstream(width=2000)

peri.teshort.upstr <- promoters(peri.teshort, upstream=2000, downstream=0)
peri.teshort.downstr <- peri.teshort %>% flank_downstream(width=2000)

peri.telong.upstr <- promoters(peri.telong, upstream=2000, downstream=0)
peri.telong.downstr <- peri.telong %>% flank_downstream(width=2000)

# export short/long TEs
export(peri.teshort, file.path(dir_features, "TEs_pericen_lt-1kb.bed"))
export(peri.telong, file.path(dir_features, "TEs_pericen_gt-2kb.bed"))


#================================
# masking features within wide SNP interval
#================================

# regions where SNP interval is long (> 4kb)
marker_v2 <- read_tsv(file.path(dir_features, "coller_gbs_marker_v2.tsv"), col_names=F)
chr_ends <- c(30427671, 19698289, 23459830, 18585056, 26975502)

snp_interval <- list()
for (i in 1:5) {
    target_chr=paste0("Chr", i)
    marker_chr <- filter(marker_v2, X1 == target_chr) %>% dplyr::select(X1, X2)
    marker_chr$end <- c(marker_chr$X2[2:nrow(marker_chr)], chr_ends[i])
    colnames(marker_chr) <- c("chr", "start", "end")
    snp_interval[[i]] <- marker_chr
}

snp_interval <- bind_rows(snp_interval) %>% mutate(width = end - start)
snp_interval.wide <- filter(snp_interval, width > 4000)
snp_interval.blacklist <- makeGRangesFromDataFrame(snp_interval.wide, starts.in.df.are.0based=F)

# export blacklist
# export(snp_interval.blacklist, file.path(dir_features, "pericen_snp_interval_gt-4000.blacklist.bed"))

#
peri.masked <- setdiff(peri, snp_interval.blacklist, ignore.strand=T)

peri.gene.masked <- setdiff(peri.gene, snp_interval.blacklist, ignore.strand=T)
peri.gene.upstr.masked <- setdiff(peri.gene.upstr, snp_interval.blacklist, ignore.strand=T)
peri.gene.downstr.masked <- setdiff(peri.gene.downstr, snp_interval.blacklist, ignore.strand=T)
# peri.te.masked <- setdiff(peri.te, snp_interval.blacklist, ignore.strand=T)
# peri.te.upstr.masked <- setdiff(peri.te.upstr, snp_interval.blacklist, ignore.strand=T)
# peri.te.downstr.masked <- setdiff(peri.te.downstr, snp_interval.blacklist, ignore.strand=T)

peri.teshort.masked <- setdiff(peri.teshort, snp_interval.blacklist, ignore.strand=T)
peri.teshort.upstr.masked <- setdiff(peri.teshort.upstr, snp_interval.blacklist, ignore.strand=T)
peri.teshort.downstr.masked <- setdiff(peri.teshort.downstr, snp_interval.blacklist, ignore.strand=T)

peri.telong.masked <- setdiff(peri.telong, snp_interval.blacklist, ignore.strand=T)
peri.telong.upstr.masked <- setdiff(peri.telong.upstr, snp_interval.blacklist, ignore.strand=T)
peri.telong.downstr.masked <- setdiff(peri.telong.downstr, snp_interval.blacklist, ignore.strand=T)




# how many CO locate in masked pericentromere?
# WT
sum(countOverlaps(co.wt.gr, peri))
# 727
sum(countOverlaps(co.wt.gr, peri.masked))
# 668 (91.88%)

# h2aw
sum(countOverlaps(co.h2aw.gr, peri))
# 2435
sum(countOverlaps(co.h2aw.gr, peri.masked))
# 2224 (91.33%)

#=================================
# count the number of crossovers overlapping regions per F2
#=================================

countFeature <- function(dat, feature){
    # dat <- co.wt.gr
    # feature <- peri.gene.upstr.masked
    count.overlap <- c()
    for (i in unique(dat$lib)){
        # i <- 407
        dat.lib <- dat[dat$lib == i]
        count.overlap <- c(count.overlap, length(unique(queryHits(findOverlaps(dat.lib, feature)))))
        # count.overlap <- c(count.overlap, sum(countOverlaps(dat.lib, feature)))
    }
    return(count.overlap)
}

## compensate one library with no crossovers
wt.peri.gene.upstr <- c(countFeature(co.wt.gr, peri.gene.upstr.masked), 0)
wt.peri.gene <- c(countFeature(co.wt.gr, peri.gene.masked), 0)
wt.peri.gene.downstr <- c(countFeature(co.wt.gr, peri.gene.downstr.masked), 0)

wt.peri.teshort.upstr <- c(countFeature(co.wt.gr, peri.teshort.upstr.masked), 0)
wt.peri.teshort <- c(countFeature(co.wt.gr, peri.teshort.masked), 0)
wt.peri.teshort.downstr <- c(countFeature(co.wt.gr, peri.teshort.downstr.masked), 0)

wt.peri.telong.upstr <- c(countFeature(co.wt.gr, peri.telong.upstr.masked), 0)
wt.peri.telong <- c(countFeature(co.wt.gr, peri.telong.masked), 0)
wt.peri.telong.downstr <- c(countFeature(co.wt.gr, peri.telong.downstr.masked), 0)

## compensate two libraries with no crossovers (h2aw-null #3 and #136)
h2aw.peri.gene.upstr <- c(countFeature(co.h2aw.gr, peri.gene.upstr.masked), 0, 0)
h2aw.peri.gene <- c(countFeature(co.h2aw.gr, peri.gene.masked), 0, 0)
h2aw.peri.gene.downstr <- c(countFeature(co.h2aw.gr, peri.gene.downstr.masked), 0, 0)

h2aw.peri.teshort.upstr <- c(countFeature(co.h2aw.gr, peri.teshort.upstr.masked), 0, 0)
h2aw.peri.teshort <- c(countFeature(co.h2aw.gr, peri.teshort.masked), 0, 0)
h2aw.peri.teshort.downstr <- c(countFeature(co.h2aw.gr, peri.teshort.downstr.masked), 0, 0)

h2aw.peri.telong.upstr <- c(countFeature(co.h2aw.gr, peri.telong.upstr.masked), 0, 0)
h2aw.peri.telong <- c(countFeature(co.h2aw.gr, peri.telong.masked), 0, 0)
h2aw.peri.telong.downstr <- c(countFeature(co.h2aw.gr, peri.telong.downstr.masked), 0, 0)

#==================================
# wilcoxon test
#==================================

n_wt <- length(unique(co.wt.gr$lib)) + 1
n_h2aw <- length(unique(co.h2aw.gr$lib)) + 2

wt.counting_result <- tibble(genotype = rep("wt", times=n_wt), peri.gene.body=wt.peri.gene, peri.gene.upstr=wt.peri.gene.upstr, peri.gene.downstr=wt.peri.gene.downstr, peri.teshort.body=wt.peri.teshort, peri.teshort.upstr=wt.peri.teshort.upstr, peri.teshort.downstr=wt.peri.teshort.downstr, peri.telong.body =wt.peri.telong, peri.telong.upstr = wt.peri.telong.upstr, peri.telong.downstr = wt.peri.telong.downstr)

h2aw.counting_result <- tibble(genotype = rep("h2aw", times=n_h2aw), peri.gene.body=h2aw.peri.gene, peri.gene.upstr=h2aw.peri.gene.upstr, peri.gene.downstr=h2aw.peri.gene.downstr, peri.teshort.body=h2aw.peri.teshort, peri.teshort.upstr=h2aw.peri.teshort.upstr, peri.teshort.downstr=h2aw.peri.teshort.downstr, peri.telong.body =h2aw.peri.telong, peri.telong.upstr = h2aw.peri.telong.upstr, peri.telong.downstr = h2aw.peri.telong.downstr)

counting_result <- bind_rows(wt.counting_result, h2aw.counting_result) %>%
    pivot_longer(cols=2:ncol(.), names_to="region.feature1.feature2", values_to="count") %>%
    separate(col=region.feature1.feature2, sep="\\.", into=c("region", "feature1", "feature2"))

counting_result.summary <- group_by(counting_result, genotype, region, feature1, feature2) %>%
    summarise(mean = mean(count), se = sd(count)/sqrt(dplyr::n())) %>%
    mutate(feature1 = factor(feature1, levels=c("gene", "teshort", "telong"))) %>%
    mutate(feature2 = factor(feature2, levels=c("upstr", "body", "downstr"))) %>%
    # mutate(feature = factor(feature, levels=c("gene_upstr", "gene", "gene_downstr", "te_upstr", "te", "te_downstr", "teshort_upstr", "teshort", "teshort_downstr", "telong_upstr", "telong", "telong_downstr"))) %>%
    mutate(genotype = factor(genotype, levels=c("wt", "h2aw"))) %>%
    arrange(genotype, region, feature1, feature2)

write_tsv(counting_result.summary, file="regions_overlapping_co-interval_summary_wt-cov-1.5.tsv")

calcWilcox <- function(dat, region.sel, feature1.sel, feature2.sel){
    # dat <- counting_result
    # region.sel <- "peri"
    # feature.sel <- "gene_upstr"
    wilcox_dat <- filter(dat, region == region.sel & feature1 == feature1.sel & feature2 == feature2.sel) %>% dplyr::select(c(genotype, count))
    pval <- wilcox.test(count ~ genotype, data=wilcox_dat)$p.value
    pval.sci <- formatC(pval, format = "e", digits=2)
    return(pval.sci)
}

pval.gene <- calcWilcox(counting_result, "peri", "gene", "body")
# 1.27e-02
pval.gene.upstr <- calcWilcox(counting_result, "peri", "gene", "upstr")
# 2.84e-05
pval.gene.downstr <- calcWilcox(counting_result, "peri", "gene", "downstr")
# 1.32e-03

pval.teshort <- calcWilcox(counting_result, "peri", "teshort", "body")
# 3.15e-03
pval.teshort.upstr <- calcWilcox(counting_result, "peri", "teshort", "upstr")
# 4.12e-04
pval.teshort.downstr <- calcWilcox(counting_result, "peri", "teshort", "downstr")
# 2.23e-06

pval.telong <- calcWilcox(counting_result, "peri", "telong", "body")
# 7.89e-01
pval.telong.upstr <- calcWilcox(counting_result, "peri", "telong", "upstr")
# 2.25e-01
pval.telong.downstr <- calcWilcox(counting_result, "peri", "telong", "downstr")
# 3.21e-01

sink(file="CO_enrichment_test_wilcox_pvalue.txt")
print("gene body")
print(pval.gene)
print("gene upstream 2kb")
print(pval.gene.upstr)
print("gene downstream 2kb")
print(pval.gene.downstr)

print("teshort body")
print(pval.teshort)
print("teshort upstream 2kb")
print(pval.teshort.upstr)
print("teshort downstream 2kb")
print(pval.teshort.downstr)

print("telong body")
print(pval.telong)
print("telong upstream 2kb")
print(pval.telong.upstr)
print("telong downstream 2kb")
print(pval.telong.downstr)
sink()



#==========================
# plotting
#==========================
my_theme <- theme_classic() +
    theme(text = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.text = element_text(size=6, family="Arial", colour="black")) +
    theme(plot.title = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.title = element_text(size=6, family="Arial", colour="black")) +
    theme(legend.position = "none") +
    # theme(legend.margin = margin(0.5,0.5,0.5,0.5)) +
    theme(strip.text = element_text(size=6, family="Arial", colour="black")) +
    theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))
    # theme(legend.spacing.y = unit(1, "cm"))
theme_set(my_theme)

plot.co_over_feature <- ggplot(counting_result.summary, aes(x=feature2, y=mean, fill=genotype)) +
    facet_grid(cols=vars(feature1)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(x=feature2, ymin=mean - se, ymax=mean + se), position=position_dodge(width=0.9), width=0.2) +
    scale_fill_jco()

fig3g <- ggplot(filter(counting_result.summary, feature1 != "te"), aes(x=feature2, y=mean, fill=genotype)) +
    facet_grid(cols=vars(feature1)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(x=feature2, ymin=mean - se, ymax=mean + se), position=position_dodge(width=0.9), linewidth=0.2, width=0.2) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(size=6, family="Arial", colour="black", angle=30, hjust=1, vjust=1)) +
    scale_fill_jco()

dir.create("results/plots", recursive=T)
pdf(file="results/plots/bar_co_over_pericentromeric_feature.pdf", width=5, height=2)
print(plot.co_over_feature)
dev.off()

pdf(file="results/plots/fig3f.pdf", width=2.8, height=1.6)
print(fig3g)
dev.off()
