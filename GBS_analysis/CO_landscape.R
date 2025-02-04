# Description 
library(tidyverse)
library(khroma)
library(scales)
library(extrafont)
library(ggsci)
library(RColorBrewer)
# install.packages("COMPoissonReg")
font_import(pattern="Arial", prompt=FALSE)
loadfonts(device="pdf")

# GENOME INFO SETTING ====================================
#' 0. chromosome info 
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)

tha.cum <- c(0, cumsum(chr.ends))
tha.tot <- tha.cum[length(tha.cum)]
centromeres.cum <- centromeres + tha.cum[1:5]
# pericentromeric + centromeric region
# (Ref.: Underwood, 2018, Genome Research)
north.start <- c(11420001, 910001, 10390001, 1070001, 8890001)
south.end <- c(18270000, 7320000, 16730000, 6630000, 15550000)

coord_subpericen <- read_csv("/home/namilhand/01_Projects/H2AW/h2aw_server/data/tair10_coord/tair10_cumulative_coord.csv")

# INPUT DATA =============================================================
wd <- "/home/namilhand/01_Projects/H2AW/h2aw_server/genome_profile/CO_profile/gbs_marker_v2/hta"
setwd(wd)
dir.create(file.path(wd, "plots"), recursive=T)

dir_gbs <- "/home/namilhand/01_Projects/H2AW/h2aw_server/data/GBS/gbs_marker_v2/tsv"
dir_chip <- "/home/namilhand/01_Projects/H2AW/h2aw_server/data/ChIP"
dir_dname <- "/home/namilhand/01_Projects/H2AW/h2aw_server/data/dname"

## GBS landscape
wt_i <- read_tsv(file.path(dir_gbs, "Selected_inhouse_WT_v2_genomeBin100kb.tsv"))
wt_r <- read_tsv(file.path(dir_gbs, "Selected_rowan_WT_v2_genomeBin100kb.tsv"))
hta67_1 <- read_tsv(file.path(dir_gbs, "hta67_set1_genomeBin100kb.tsv"))
hta67_2 <- read_tsv(file.path(dir_gbs, "hta67_set2_genomeBin100kb.tsv"))
hta6_1 <- read_tsv(file.path(dir_gbs, "hta6_set1_genomeBin100kb.tsv"))
hta6_2 <- read_tsv(file.path(dir_gbs, "hta6_set2_genomeBin100kb.tsv"))
hta7_1 <- read_tsv(file.path(dir_gbs, "hta7_set1_genomeBin100kb.tsv"))
hta7_2 <- read_tsv(file.path(dir_gbs, "hta7_set2_genomeBin100kb.tsv"))
hta12_1 <- read_tsv(file.path(dir_gbs, "hta12_set1_genomeBin100kb.tsv"))
hta12_2 <- read_tsv(file.path(dir_gbs, "hta12_set2_genomeBin100kb.tsv"))
h2aw_1 <- read_tsv(file.path(dir_gbs, "h2aw-null_set1_genomeBin100kb.tsv"))
h2aw_2 <- read_tsv(file.path(dir_gbs, "h2aw-null_set2_genomeBin100kb.tsv"))
h2aw_3 <- read_tsv(file.path(dir_gbs, "h2aw-null_set3_genomeBin100kb.tsv"))

## ChIPseq landscape
mnase_chip <- read_tsv(file.path(dir_chip, "MNase_WT-bud_log2ChIP_binSize100kb.tsv"))
k9me2_chip <- read_tsv(file.path(dir_chip, "K9me2_WT-bud_log2ChIP_binSize100kb.tsv"))
k9me2_chip <- add_column(k9me2_chip, feature = "H3K9me2")

h2aw6_chip <- read_tsv(file.path(dir_chip, "H2AW6_WT_seedling_log2ChIP_binSize100kb.tsv"))
h2aw7_chip <- read_tsv(file.path(dir_chip, "H2AW7_WT_seedling_log2ChIP_binSize100kb.tsv"))

## dname landscape
mC_col <- read_tsv(file.path(dir_dname, "20230625_bsseq_migs_hta/bedg/tiled/window_100kb_step_100kb", "1-1_C_TAIR10_window100kb_step100kb.bedg"), col_names=c("chr", "window", "window_end", "mC.ratio"))
### add cumulative window
mC_col <- mC_col %>%
    mutate(cumwindow = case_when(chr == "Chr1" ~ window + tha.cum[1],
                                chr == "Chr2" ~ window + tha.cum[2],
                                chr == "Chr3" ~ window + tha.cum[2],
                                chr == "Chr4" ~ window + tha.cum[2],
                                chr == "Chr5" ~ window + tha.cum[2])) %>%
    relocate(cumwindow, .after=window_end) %>%
    filter(chr %in% paste0("Chr", 1:5))

## merge GBS results from the same genotype
wt_tsv <- dplyr::select(wt_i, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = wt_i$coInWindow + wt_r$coInWindow) %>%
    add_column(libSize = wt_i$libSize + wt_r$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

hta67_tsv <- dplyr::select(hta67_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = hta67_1$coInWindow + hta67_2$coInWindow) %>%
    add_column(libSize = hta67_1$libSize + hta67_2$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

hta6_tsv <- dplyr::select(hta6_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = hta6_1$coInWindow + hta6_2$coInWindow) %>%
    add_column(libSize = hta6_1$libSize + hta6_2$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

hta7_tsv <- dplyr::select(hta7_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = hta7_1$coInWindow + hta7_2$coInWindow) %>%
    add_column(libSize = hta7_1$libSize + hta7_2$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

hta12_tsv <- dplyr::select(hta12_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = hta12_1$coInWindow + hta12_2$coInWindow) %>%
    add_column(libSize = hta12_1$libSize + hta12_2$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

h2aw_tsv <- dplyr::select(h2aw_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = h2aw_1$coInWindow + h2aw_2$coInWindow + h2aw_3$coInWindow) %>%
    add_column(libSize = 287) %>%
    # add_column(libSize = h2aw_1$libSize + h2aw_2$libSize + h2aw_3$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

# co_tsv.list <- list(wt=col_tsv, h2aw=h2aw_tsv, suvh=suvh_tsv, dmigs=dmigs_tsv, cmt3=cmt3_tsv)
co_tsv.list <- list(wt=wt_tsv, hta6=hta6_tsv, hta7=hta7_tsv, hta12=hta12_tsv, hta67=hta67_tsv, h2aw=h2aw_tsv)


# cumulative coordinates of pericentromere
coord_pericen <- tibble(chr=1:5, north=north.start + tha.cum[1:5], south=south.end + tha.cum[1:5])

#===== MA smoothing ====
# k = 5 is recommended for drawing ChIP landscape

ma_width.co <- 7
ma_width.chip <- 5

mafilter <- function(dat, valueColumn, k){
        # k = arm width for MA smoothing
        filt <- 1/(2*k+1)
        filt <- rep(filt, 2*k+1)
        filt.collect <- NULL

        for(i in 1:5){
                chr <- dat %>%
                        filter(chr==paste0("Chr",i))
                filt.chr <- stats::filter(chr[[valueColumn]], filt)
                filt.chr[1:k] <- filt.chr[k+1]
                len <- length(filt.chr)
                filt.chr[(len-k+1):len] <- filt.chr[(len-k)]

                filt.collect <- c(filt.collect, filt.chr)
        }
        fin <- tibble(dat, smooth=filt.collect)
        return(fin)
}

# moving average genomeBin100kb
co_ma.list <- lapply(co_tsv.list, mafilter, "mean.coInWindow", ma_width.co)
co_ma <- bind_rows(co_ma.list, .id="genotype")
co_ma <- mutate(co_ma, genotype = factor(genotype, levels=c("wt", "hta6", "hta7", "hta12", "hta67", "h2aw")))  %>%
    mutate(cMMb = smooth *10 * 100/2)
    
mnase_ma <- mafilter(mnase_chip, "cov", ma_width.chip)
k9me2_ma <- mafilter(k9me2_chip, "cov", ma_width.chip)
h2aw6_ma <- mafilter(h2aw6_chip, "cov", ma_width.chip)
h2aw7_ma <- mafilter(h2aw7_chip, "cov", ma_width.chip)
h2aw_ma <- bind_rows(list(H2AW.6=h2aw6_ma, H2AW.7=h2aw7_ma), .id="feature")
h2aw_k9me2_ma <- bind_rows(list(H2AW.6=h2aw6_ma, H2AW.7=h2aw7_ma, H3K9me2=k9me2_ma), .id="feature")

# # CO PCA plot
# co_pca_input <- dplyr::select(co_ma, c("genotype", "cumwindow", "cMMb")) %>%
#     pivot_wider(names_from="cumwindow", values_from="cMMb") %>% as.data.frame

# pc <- prcomp(co_pca_input[,-1], center=T)
# attributes(pc)

# library(devtools)
# install_github("vqv/ggbiplot")
# library(ggbiplot)
# g <- ggbiplot(pc, obs.scale = 1, var.scale = 1, groups = co_pca_input$genotype, elipse=T, circle=F, var.axes=F)

# co_hclust_dist <- dist(co_pca_input[,-1])
# co_hclust <- hclust(co_hclust_dist, method = "average")
# plot(co_hclust)
# clusterCut <- cutree(co_hclust)

#===== DRAWING CO LANDSCAPE =====
## palette
pal_GBS <- c("black", pal_d3(palette = c("category10"))(5))
names(pal_GBS) <- c("wt","hta6", "hta67", "hta7", "h2aw", "hta12")
# names(pal_GBS) <- c("wt", "h2aw", "dmigs", "cmt3", "suvh")
pal_chip <- "blue"

# co landscape

## y-axis scale
# landscape_ylim <- c(0, max(co_ma$cMMb)*1.05) 
landscape_ylim <- c(0, 11.8)

drawCoLandscape <- function(dat, pal, prefix, plot_width, plot_height){

    ## convert CO per F2 per 100 kb into cM/Mb
    ## *10 : CO per F2 per 100kb --> CO per F2 per Mb
    ## *100: CO per F2 per Mb (M/Mb) --> cM/Mb

    # dat_cMMb <- mutate(dat, cMMb = smooth * 10 * 100 )
    # drawCoLandscape(co_ma, pal_GBS, "plots/mmH-hta", 3.3, 1.0)
    # dat <- filter(co_ma, genotype %in% c("wt", "hta6"))
    # pal <- pal_GBS


    p <- ggplot() +
        annotate("rect", xmin=coord_pericen$north[1], xmax=coord_pericen$south[1], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[2], xmax=coord_pericen$south[2], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[3], xmax=coord_pericen$south[3], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[4], xmax=coord_pericen$south[4], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[5], xmax=coord_pericen$south[5], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        # geom_segment(data=syri.cum, aes(x=start, xend=end, y=Inf, yend=Inf), size=5) +
        # geom_segment(data=nrz.cum, aes(x=start, xend=end, y=20, yend=20), size=3, colour="red") +
        geom_line(data=dat, aes(x=cumwindow, y=cMMb, colour=genotype), linewidth=0.2) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mnase_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
        scale_y_continuous(limits=landscape_ylim) +
        scale_colour_manual(values = pal) +
        geom_vline(xintercept=tha.cum[2:5], colour="Black", linewidth=0.2) +
        geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", linewidth=0.2) +
        labs(y="cM/Mb") +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
                legend.title=element_text(size=6),
                legend.text=element_text(size=6)) +
        theme(legend.title=element_blank()) +
        theme(legend.position = c(1, 1)) +
        theme(legend.justification = c(1,0.8)) +
        theme(legend.box.margin = margin(0.5,0.5,0.5,0.5)) +
        theme(text=element_text(size=6, family="Arial", colour="black"),
        axis.text=element_text(colour="black", size=6)) +
        theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))

        # pdf(file="hta6_co_landscape_with_SV_CEN_annotated.pdf", width=, 4, height=1.5)
        # print(p)
        # dev.off()

        pdf(file=paste0(prefix, "_co_landscape_ma-co_", ma_width.co,".pdf"), width=plot_width, height=plot_height)
        print(p)
        dev.off()
        png(file=paste0(prefix, "_co_landscape_ma-co_", ma_width.co,".png"), width=plot_width, height=plot_height, unit="in", res=300)
        print(p)
        dev.off()
}

## all merged
drawCoLandscape(co_ma, pal_GBS, "plots/hta", 3.3, 1.0)
## without cmt3
# drawCoLandscape(filter(co_ma, genotype != "cmt3"), pal_GBS, "plots/h2aw_dmigs", 5.7, 2.0)

## one mutant vs WT
# drawCoLandscape(filter(co_ma, genotype %in% c("wt", "h2aw")), pal_GBS, "plots/h2aw", 5.7, 2.0)
drawCoLandscape(filter(co_ma, genotype %in% c("wt", "hta6")), pal_GBS, "plots/hta6", 3.3, 1)
drawCoLandscape(filter(co_ma, genotype %in% c("wt", "hta7")), pal_GBS, "plots/hta7", 3.3, 1)
drawCoLandscape(filter(co_ma, genotype %in% c("wt", "hta12")), pal_GBS, "plots/hta12", 3.3, 1)
drawCoLandscape(filter(co_ma, genotype %in% c("wt", "hta67")), pal_GBS, "plots/hta67", 3.3, 1)
drawCoLandscape(filter(co_ma, genotype %in% c("wt", "h2aw")), pal_GBS, "plots/h2aw", 3.3, 1)

drawCoDiffLandscape <- function(dat, gt_test, gt_ctrl, linecolour, y_label, prefix, plot_width, plot_height, y_min, y_max){
    # dat <- co_ma
    # linecolour <- pal_GBS[c("cmt3", "h2aw")]
    # dat_cMMb <- mutate(dat, cMMb = smooth * 10 * 100 )

    # dat <- filter(co_ma, genotype %in% c("wt", "hta6"))
    # linecolour <- pal_GBS
    # gt_test <- "hta6"
    # gt_ctrl <- "wt"
    # y_label <- "-WT"
    # y_min <- -10.5
    # y_max <- 10

    diff_landscape_ylim <- c(y_min, y_max)

    dat_ctrl <- filter(dat, genotype == gt_ctrl)
    dat_diff <- list()

    for(i in 1:length(gt_test)){
        dat_test <- filter(dat, genotype == gt_test[i])
        dat_diff[[i]] <- dplyr::select(dat_ctrl, !c("genotype", "smooth", "cMMb"))
        dat_diff[[i]] <- add_column(dat_diff[[i]], difference =  dat_test$cMMb - dat_ctrl$cMMb)
    }
    names(dat_diff) <- gt_test
    dat_diff_bind <- bind_rows(dat_diff, .id="genotype")

    if(y_label == "wt"){
        y_label <- "WT"
    }
    y_label <- bquote(Delta ~ " cM/Mb (" ~ .(y_label) ~ ")")

    p <- ggplot() +
        annotate("rect", xmin=coord_pericen$north[1], xmax=coord_pericen$south[1], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[2], xmax=coord_pericen$south[2], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[3], xmax=coord_pericen$south[3], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[4], xmax=coord_pericen$south[4], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[5], xmax=coord_pericen$south[5], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        # geom_segment(data=syri.cum, aes(x=start, xend=end, y=Inf, yend=Inf), size=5) +
        # geom_segment(data=nrz.cum, aes(x=start, xend=end, y=8, yend=8), size=3, colour="red") +
        geom_line(data=dat_diff_bind, aes(x=cumwindow, y=difference, colour=genotype), linewidth=0.2) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mnase_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
        scale_y_continuous(limits=diff_landscape_ylim) +
        scale_colour_manual(values=linecolour) +
        geom_vline(xintercept=tha.cum[2:5], colour="Black", linewidth=0.2) +
        geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", linewidth=0.2) +
        geom_hline(yintercept=0, linewidth=0.2, colour="Black") +
        labs(y=y_label) +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6)) +
        theme(legend.title=element_blank()) +
        theme(legend.position = c(1, 1)) +
        theme(legend.justification = c(1,0.8)) +
        theme(legend.box.margin = margin(0.5,0.5,0.5,0.5)) +
        theme(text=element_text(size=6, family="Arial", colour="black"),
        axis.text=element_text(colour="black", size=6)) +
        theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))

        # pdf(file="hta6-wt_co_landscape_with_SV_CEN_annotated.pdf", width=, 4, height=1.5)
        # print(p)
        # dev.off()

        pdf(file=paste0(prefix, "_co-diff_landscape_ma-co_", ma_width.co,".pdf"), width=plot_width, height=plot_height)
        print(p)
        dev.off()
        png(file=paste0(prefix, "_co-diff_landscape_ma-co_", ma_width.co,".png"), width=plot_width, height=plot_height, unit="in", res=300)
        print(p)
        dev.off()
}

drawCoDiffLandscape(co_ma, "hta6", "wt", pal_GBS, "-WT", "plots/hta6-wt", 3.3, 1, -4, 4)
drawCoDiffLandscape(co_ma, "hta7", "wt", pal_GBS, "-WT", "plots/hta7-wt", 3.3, 1, -4, 4)
drawCoDiffLandscape(co_ma, "hta12", "wt", pal_GBS, "-WT", "plots/hta12-wt", 3.3, 1, -4, 4)
drawCoDiffLandscape(co_ma, "hta67", "wt", pal_GBS, "-WT", "plots/hta67-wt", 3.3, 1, -4, 4)
drawCoDiffLandscape(co_ma, "h2aw", "wt", pal_GBS, "-WT", "plots/h2aw-wt", 3.3, 1, -4, 4)
# drawCoDiffLandscape(co_ma, c("hta6", "hta7", "hta12", "hta67"), "wt", pal_GBS, "-WT", "plots/hta-wt", 5.7, 2, -4, 4)
# drawCoDiffLandscape(co_ma, "suvh", "wt", pal_GBS["suvh"], "-WT", "plots/suvh-wt", 5.7, 2, -10.5, 10)
# drawCoDiffLandscape(co_ma, "dmigs", "wt", pal_GBS["dmigs"], "-WT", "plots/dmigs-wt", 5.7, 2, -10.5, 10)
# drawCoDiffLandscape(co_ma, "dmigs", "h2aw", pal_GBS["dmigs"], "-h2aw", "plots/dmigs-h2aw", 5.7, 2, -10.5, 10)

#==== DRAW CHROMATIN LANDSCAPE =====
drawChipLandscape <- function(dat, prefix, y_label, plot_width, plot_height){
#    dat <- h2aw_ma
    pal <- colour("high contrast")(length(unique(dat$feature)))
    names(pal) <- NULL

    p <- ggplot() +
        annotate("rect", xmin=coord_pericen$north[1], xmax=coord_pericen$south[1], ymin=-Inf, ymax=Inf, fill="lightskyblue1") +
        annotate("rect", xmin=coord_pericen$north[2], xmax=coord_pericen$south[2], ymin=-Inf, ymax=Inf, fill="lightskyblue1") +
        annotate("rect", xmin=coord_pericen$north[3], xmax=coord_pericen$south[3], ymin=-Inf, ymax=Inf, fill="lightskyblue1") +
        annotate("rect", xmin=coord_pericen$north[4], xmax=coord_pericen$south[4], ymin=-Inf, ymax=Inf, fill="lightskyblue1") +
        annotate("rect", xmin=coord_pericen$north[5], xmax=coord_pericen$south[5], ymin=-Inf, ymax=Inf, fill="lightskyblue1") +
        annotate("rect", xmin=coord_subpericen$ColdPeri.start[1], xmax=coord_subpericen$ColdPeri.end[1], ymin=-Inf, ymax=Inf, fill="grey85") +
        annotate("rect", xmin=coord_subpericen$ColdPeri.start[2], xmax=coord_subpericen$ColdPeri.end[2], ymin=-Inf, ymax=Inf, fill="grey85") +
        annotate("rect", xmin=coord_subpericen$ColdPeri.start[3], xmax=coord_subpericen$ColdPeri.end[3], ymin=-Inf, ymax=Inf, fill="grey85") +
        annotate("rect", xmin=coord_subpericen$ColdPeri.start[4], xmax=coord_subpericen$ColdPeri.end[4], ymin=-Inf, ymax=Inf, fill="grey85") +
        annotate("rect", xmin=coord_subpericen$ColdPeri.start[5], xmax=coord_subpericen$ColdPeri.end[5], ymin=-Inf, ymax=Inf, fill="grey85") +
        geom_line(data=dat, aes(x=cumwindow, y=smooth, colour=feature), linewidth=0.3) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mnase_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
        scale_colour_manual(values = pal) +
        geom_vline(xintercept=tha.cum[2:5], colour="Black", linewidth=0.2) +
        geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", linewidth=0.2) +
        labs(y=y_label) +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
        legend.title=element_text(size=7),
        legend.text=element_text(size=7)) +
        theme(legend.title=element_blank()) +
        theme(legend.position = c(1, 1)) +
        theme(legend.justification = c(1,0.8)) +
        theme(legend.box.margin = margin(0.5,0.5,0.5,0.5)) +
        theme(text=element_text(size=9, family="Helvetica", colour="black"),
        axis.text=element_text(colour="black", size=7))

        pdf(file=paste0(prefix, "_chip_landscape.pdf"), width=plot_width, height=plot_height)
        print(p)
        dev.off()
        png(file=paste0(prefix, "_chip_landscape.png"), width=plot_width, height=plot_height, unit="in", res=300)
        print(p)
        dev.off()
}

drawChipLandscape(h2aw_ma, "H2AW_seedling", "log2[ChIP/input]", 3.3, 1.2)
drawChipLandscape(k9me2_ma, "H3K9me2_bud", "log2[ChIP/input]", 3.3, 1.2)
drawChipLandscape(h2aw_k9me2_ma, "H2AW-seedling_H3K9me2-bud", "log2[ChIP/input]", 3.3, 1.2)



#' 7 TEL-CEN plotting
# distFromTel() : calculate the distance from telomere in two directions (north or south).
# It takes *.cotable as input

#-----Note-----
# This code is kind of trikcy at a first glance.
# The first idea of TEL-CEN scaling that pop in my mind is below:
#   1. scale the CO coordinate into proportion
#   2. tile the proportion into bin
#   However, this does not work because when you scale the coordinate before tiling,
# you will eventually mixup the CO coordinates completely. It hides the pattern.
#   Therefore, before converting into proportion, just tile first so that the CO sites fromthe same chromosome can be tied together a little. Maybe it's the point where the method needs to be modified. 
#---------------


# distFromTel() : convert coordinate into proportion of distance from cent to tel
distFromTel <- function(dat, binsize, spline_df, target_column){
		# dat <- co_tsv.list$wt
        # spline_df <- 11
        # target_column <- "mean.coInWindow"
		# binsize <- 0.01
        # co_telcen.list <- lapply(co_tsv.list, distFromTel, 0.01, spline_df, "mean.coInWindow")

        dat <- arrange(dat, chr, window)
        dat$cent <- rep(centromeres, times=table(dat$chr))
        dat$end <- rep(chr.ends, times=table(dat$chr))
        
        left <- filter(dat, window < cent) %>%
                add_column(arm="north") %>%
                mutate(coord.prop=window/cent) %>%
				filter(chr %in% c("Chr1", "Chr3", "Chr5"))

		# Excluded north arm of Chr2 and Chr4 as these regions is too short for fair comparison with other chromosomes
        right <- filter(dat, window >= cent) %>%
                add_column(arm="south") %>%
                mutate(coord.prop=(end-window)/(end-cent))
        
        prop <- bind_rows(left, right) %>%
            arrange(coord.prop)

        # bin.spline <- smooth.spline(prop$coord.prop, prop$mean.coInWindow, nknots=10)

        # plot(prop$coord.prop, prop$mean.coInWindow)
        # lines(prop$coord.prop, c(bin.spline$y, 0), lty=2, col=2, lwd=2)
        
        collect_bin <- NULL
        wins <- seq(0, 1, by=binsize)

        for(j in 1:(length(wins)-1)){
                bin <- mean(prop[which(prop$coord.prop >= wins[j] & prop$coord.prop < wins[j+1]),][[target_column]])
                collect_bin <- c(collect_bin, bin)
        }

        # collect_bin.smooth <- smooth.spline(wins[2:length(wins)], collect_bin, nknots=n_knots)$y
        collect_bin.smooth <- smooth.spline(wins[2:length(wins)], collect_bin, df=spline_df)$y
        collect_bin.cmmb <- collect_bin.smooth/2*10*100

        # k=mafilt.size
        # filt <- 1/(2*k+1)
        # filt <- rep(filt, 2*k+1)
        # filt.dat <- stats::filter(collect_bin, filt)
        # filt.dat <- loess(formula = paste(target_column, "coord.prop", sep="~"), span=0.1, degree=1, data=prop)
        # filt.dat <- filt.dat$fitted

        res <- tibble(prop=wins[2: length(wins)],
                        bin=collect_bin,
                        bin.smooth=collect_bin.smooth,
                        bin.cmmb=collect_bin.cmmb
                    #     bin=collect_bin,
                    #   bin.ma=filt.dat,
                      )
        
        return(res)
}

spline_df <- 11
co_telcen.list <- lapply(co_tsv.list, distFromTel, 0.01, spline_df, "mean.coInWindow")
co_telcen.bind <- bind_rows(co_telcen.list, .id="sample") %>%
    mutate(sample = factor(sample, levels=c("wt", "hta6", "hta7", "hta12", "hta67", "h2aw")))
k9me2_telcen <- distFromTel(k9me2_chip, 0.01, spline_df, "cov")
h2aw6_telcen <- distFromTel(h2aw6_chip, 0.01, spline_df, "cov")
h2aw7_telcen <- distFromTel(h2aw7_chip, 0.01, spline_df, "cov")
mC_col_telcen <- distFromTel(mC_col, 0.01, spline_df, "mC.ratio")


## define centromeric-pericentromeric region
## region where DNA methylation level is higher than the average
define_telcen_pericen <- mC_col_telcen %>%
    mutate(distFromAvg = abs(bin.cmmb - mean(bin.cmmb))) %>%
    arrange(distFromAvg)
telcen_pericen_border <- define_telcen_pericen$prop[1]

## define tel-cen subpericentromeric region
define_subpericen_border <- co_telcen.list[["wt"]] %>%
    mutate(distFromAvg = abs(bin.cmmb - mean(bin.cmmb))) %>%
    filter(prop > telcen_pericen_border) %>%
    # mutate(meanco=mean(bin.cmmb))
    arrange(distFromAvg)

subpericen_border <- define_subpericen_border$prop[1]

drawTelCenProfile <- function(dat_co, dat_chromatin, chrom_label, prefix, plot_width, plot_height){
    # dat_co <- filter(co_telcen.bind, sample != "cmt3")
    # dat_chromatin <- mC_col_telcen
    # dat_chromatin <- k9me2_telcen
    # chrom_label <- "mC ratio"
    # chrom_label <- "k9me2"



    pal <- pal_GBS

    # scale two different y-axis
    range_co <- summary(dat_co$bin.cmmb)["Max."] - summary(dat_co$bin.cmmb)["Min."]
    range_chromatin <- summary(dat_chromatin$bin.smooth)["Max."] - summary(dat_chromatin$bin.smooth)["Min."]
    scale_y <- range_co/range_chromatin * 0.8
    # If smoothed mean CO is negative value, offset that negative value from y-scale of chromatin by adding the minimum smoothed mean CO
    offset_chromatin = 0
    if(min(dat_co$bin.cmmb) < 0){
        offset_chromatin <- min(dat_co$bin.cmmb)
    }
    # shift_y <- (min(dat_co$bin.smooth, na.rm=T) - min(dat_chromatin$bin.smooth, na.rm=T)*scale_y) - offset_chromatin
    # shift_y <- -offset_chromatin
    if(min(dat_chromatin$bin.smooth) < 0){
        shift_y <- -min(dat_chromatin$bin.smooth)
    } else {
        shift_y <- 0
    }

    # re-scaled value of chromatin = bin.smooth * scale_y + shift_y - offset_chromatin
    # min(dat_co$bin.smooth, na.rm=T) - min(dat_chromatin$bin.smooth, na.rm=T)*scale_y

    mean_cmmb <- dat_co %>%
            group_by(sample) %>%
            summarise(mean.cmmb=mean(bin, na.rm=TRUE)*500)
    
    yscale_lim <- c(min(dat_co$bin.cmmb, na.rm=T)*0.9, max(dat_co$bin.cmmb, na.rm=T)*1.1)
    # p <- ggplot() +
    #     # annotate("rect", xmin=telcen_pericen_border, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey90") +
    #     geom_vline(xintercept=telcen_pericen_border, size=0.2, linetype = "dashed", colour="grey60") +
    #     geom_vline(xintercept=subpericen_border, size=0.2, linetype = "dashed", colour="grey60") +
    #     geom_line(data=dat_co, aes(x=prop, y=bin.smooth, colour=sample), linewidth=0.2) +
    #     # geom_smooth(data=dat_co, aes(x=prop, y=bin, colour=sample), size=0.4, method="lm", formula = y ~ splines::bs(x, 3), se=FALSE) +
    #     geom_hline(data=mean_cmmb, aes(yintercept=mean.cmmb, colour=sample), linetype="dashed", size=0.2) +
    #     geom_area(data=dat_chromatin, aes(x=prop, y=(bin*scale_y) + shift_y), fill="grey80", alpha=0.5) +
    #     scale_y_continuous(name = "cM/Mb", sec.axis=sec_axis(~(. - shift_y)/scale_y, name=chrom_label)) +
    #     coord_cartesian(ylim=yscale_lim) +
    #     scale_x_continuous(name="Distance from the telomere (ratio)",
    #                        breaks=seq(0, 1, by=0.25),
    #                        label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
    #     scale_colour_manual(values=pal) +
    #     theme_classic() +
    #     theme(legend.key.size=unit(0.1, "inches"),
    #           legend.title=element_text(size=6),
    #           legend.text=element_text(size=6),
    #           legend.position=c(0,1),
    #           legend.justification=c(-0.2, 0.8)) +
    #     theme(legend.title = element_blank()) +
    #     theme(text=element_text(size=6, family="Arial", colour="black"),
    #     axis.text=element_text(size=6, colour="black"))
    p <- ggplot() +
        geom_vline(xintercept=telcen_pericen_border, size=0.2, linetype = "dashed", colour="black") +
        geom_line(data=dat_co, aes(x=prop, y=bin.cmmb, colour=sample), linewidth=0.2) +
        geom_hline(data=mean_cmmb, aes(yintercept=mean.cmmb, colour=sample), linetype="dashed", size=0.2) +
        geom_area(data=dat_chromatin, aes(x=prop, y=(bin.smooth + shift_y) * scale_y), fill="grey80", alpha=0.5) +
        scale_y_continuous(name = "cM/Mb", sec.axis=sec_axis(~(./scale_y) - shift_y, name=chrom_label)) +
        coord_cartesian(ylim=yscale_lim) +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        scale_colour_manual(values=pal) +
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=6),
              legend.text=element_text(size=6),
              legend.position=c(0,1),
              legend.justification=c(-0.2, 0.8)) +
        theme(legend.title = element_blank()) +
        theme(text=element_text(size=6, family="Arial", colour="black"),
        axis.text=element_text(size=6, colour="black")) +
        theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))
    
    pdf(file=paste0("plots/telcen_", prefix, "_spline-df_", spline_df, ".pdf"), plot_width, plot_height)
    print(p)
    dev.off()

    png(file=paste0("plots/telcen_", prefix, "_spline-df_", spline_df, ".png"), plot_width, plot_height, unit="in", res=300)
    print(p)
    dev.off()
}
drawTelCenProfile(filter(co_telcen.bind), h2aw6_telcen, "log2[H2AW.6/input]", "mmH-hta_with_H2AW6", 2.4, 1.4)
drawTelCenProfile(filter(co_telcen.bind), h2aw7_telcen, "log2[H2AW.7/input]", "mmH-hta_with_H2AW7", 2.4, 1.4)
drawTelCenProfile(filter(co_telcen.bind), k9me2_telcen, "log2[H3K9me2/input]", "mmH-hta_with_K9me2", 2.4, 1.4)
drawTelCenProfile(filter(co_telcen.bind), mC_col_telcen, "DNA methylation (ratio)", "mmH-hta_with_dname", 2.3, 2.08)

drawTelCenDiffProfile <- function(dat_co, ctrl_group, prefix, plot_width, plot_height){
    # dat_co <- filter(co_telcen.bind, sample != "cmt3")
    # ctrl_group <- "wt"
    # dat_chromatin <- mC_col_telcen
    # dat_chromatin <- k9me2_telcen
    # chrom_label <- "mC ratio"
    # chrom_label <- "k9me2"

    dat_co.ctrl <- filter(dat_co, sample == ctrl_group)
    dat_co.diff <- filter(dat_co, sample != ctrl_group) %>%
        add_column(bin.cmmb.wt = rep(dat_co.ctrl$bin.cmmb, times=length(unique(.$sample)))) %>%
        mutate(bin.cmmb.diff = bin.cmmb - bin.cmmb.wt)

    pal <- pal_GBS

    # yscale_lim <- c(min(dat_co.diff$bin.cmmb.diff, na.rm=T)*0.9, max(dat_co.diff$bin.cmmb.diff, na.rm=T)*1.1)
    yscale_lim <- c(-2, 1.5)

    if(ctrl_group == "wt"){
        ctrl_group = "WT"
    }
    y_label <- paste0("-", ctrl_group)
    y_label <- bquote(Delta ~ " cM/Mb (" ~ .(y_label) ~ ")")



    p <- ggplot() +
        # annotate("rect", xmin=telcen_pericen_border, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey90") +
        geom_vline(xintercept=telcen_pericen_border, size=0.2, linetype = "dashed", colour="black") +
        geom_line(data=dat_co.diff, aes(x=prop, y=bin.cmmb.diff, colour=sample), size=0.4, linewidth=0.2) +
        geom_hline(yintercept=0, linewidth=0.2, colour="Black") +
        # geom_smooth(data=dat_co, aes(x=prop, y=bin, colour=sample), size=0.4, method="lm", formula = y ~ splines::bs(x, 3), se=FALSE) +
        # geom_hline(data=mean_cos, aes(yintercept=meanCO, colour=sample), linetype="dashed", size=0.3) +
        # geom_area(data=dat_chromatin, aes(x=prop, y=(bin*scale_y) + shift_y), fill="grey90", alpha=0.5) +
        scale_y_continuous(name = y_label) +
        # coord_cartesian(ylim=yscale_lim) +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        scale_colour_manual(values=pal) +
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=6),
              legend.text=element_text(size=6),
              legend.position=c(0,1),
              legend.justification=c(-0.2, 0.8)) +
        theme(legend.title = element_blank()) +
        theme(text=element_text(size=6, family="Arial", colour="black"),
        axis.text=element_text(size=6, colour="black")) +
        theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))
    
    pdf(file=paste0("plots/telcen-diff_", prefix, "_spline-df_", spline_df, ".pdf"), plot_width, plot_height)
    print(p)
    dev.off()

    png(file=paste0("plots/telcen-diff_", prefix, "_spline-df_", spline_df, ".png"), plot_width, plot_height, unit="in", res=300)
    print(p)
    dev.off()
}
drawTelCenDiffProfile(co_telcen.bind, "wt", "ctrl=wt", 2.2, 2.08)
# drawTelCenDiffProfile(co_telcen.bind, "hta6", "ctrl=hta6", 1.8, 1.4)
# drawTelCenDiffProfile(co_telcen.bind, "hta7", "ctrl=hta7", 1.8, 1.4)
# drawTelCenDiffProfile(co_telcen.bind, "hta12", "ctrl=hta12", 1.8, 1.4)
# drawTelCenDiffProfile(co_telcen.bind, "hta67", "ctrl=hta67", 1.8, 1.4)


# ========================== ARCHIVE ========================
delH2AW <- left_join(cos.all.list.bin[[1]], cos.all.list.bin[[2]], by=c("chrs", "bin.start", "bin.end", "cum.start", "cum.end")) %>%
    mutate(delCO=bin.y - bin.x) %>%
    dplyr::select(-c("bin.x", "bin.y"))

delH2AW_mnase <- dplyr::select(diff_h2aw, c(chr, cumwindow)) %>%
    add_column(diff.h2aw=diff_h2aw$smooth.diff) %>%
    add_column(smooth.mnase=mnase_ma$smooth)

north.start.cum <- north.start + tha.cum[1:5]
south.end.cum <- south.end + tha.cum[1:5]

delH2AW_mnase_arm <- filter(delH2AW_mnase, !(cumwindow > north.start.cum[1] & cumwindow < south.end.cum[1])) %>%
                    filter(!(cumwindow > north.start.cum[2] & cumwindow < south.end.cum[2])) %>%
                    filter(!(cumwindow > north.start.cum[3] & cumwindow < south.end.cum[3])) %>%
                    filter(!(cumwindow > north.start.cum[4] & cumwindow < south.end.cum[4])) %>%
                    filter(!(cumwindow > north.start.cum[5] & cumwindow < south.end.cum[5]))
delH2AW_mnase_peri <- filter(delH2AW_mnase, (cumwindow > north.start.cum[1] & cumwindow < south.end.cum[1]) |
                    (cumwindow > north.start.cum[2] & cumwindow < south.end.cum[2]) |
                    (cumwindow > north.start.cum[3] & cumwindow < south.end.cum[3]) |
                    (cumwindow > north.start.cum[4] & cumwindow < south.end.cum[4]) |
                    (cumwindow > north.start.cum[5] & cumwindow < south.end.cum[5]))

                    filter(!(cum.start > north.start.cum[2] & cum.end < south.end.cum[2])) %>%
                    filter(!(cum.start > north.start.cum[3] & cum.end < south.end.cum[3])) %>%
                    filter(!(cum.start > north.start.cum[4] & cum.end < south.end.cum[4])) %>%
                    filter(!(cum.start > north.start.cum[5] & cum.end < south.end.cum[5]))

lm.delCO_by_mnase_arm <- lm(smooth.mnase ~ diff.h2aw, data=delH2AW_mnase_arm)

p_delCO_by_mnase <- ggplot(delH2AW_mnase_arm, aes(x=smooth.mnase, y=diff.h2aw)) +
    geom_point() +
    geom_smooth(method='lm') +
    scale_x_continuous(name="Nucleosome [MNase/gDNA]") +
    scale_y_continuous(name="Diff CO [migsH2AW - WT]")

pdf(file="h2aw-diffCO_by_mnase.pdf")
print(p_delCO_by_mnase)
dev.off()

# TEL-CEN proportional distance smoothed by MOVING AVERAGE.
distFromTel <- function(dat, binsize, mafilt.size, target_column){
#		dat <- cos.all.list.bin100$hcr2
#		binsize <- 0.01
#		mafilt.size <- 9
        dat <- filter(co_ma, genotype == "wt")
        # dat <- k9me2_col
        target_column <- "mean.coInWindow"
        # target_column <- "cov"
        binsize <- 0.01
        # mafilt.size <- 7
        n_knots <- 10

        dat <- arrange(dat, chr, window)
        dat$cent <- rep(centromeres, times=table(dat$chr))
        dat$end <- rep(chr.ends, times=table(dat$chr))
        
        left <- filter(dat, window < cent) %>%
                add_column(arm="north") %>%
                mutate(coord.prop=window/cent) %>%
				filter(chr %in% c("Chr1", "Chr3", "Chr5"))

		# Excluded north arm of Chr2 and Chr4 as these regions is too short for fair comparison with other chromosomes
        right <- filter(dat, window >= cent) %>%
                add_column(arm="south") %>%
                mutate(coord.prop=(end-window)/(end-cent))
        
        prop <- bind_rows(left, right) %>%
            arrange(coord.prop)

        # bin.spline <- smooth.spline(prop$coord.prop, prop$mean.coInWindow, nknots=10)

        # plot(prop$coord.prop, prop$mean.coInWindow)
        # lines(prop$coord.prop, c(bin.spline$y, 0), lty=2, col=2, lwd=2)
        
        collect_bin <- NULL
        wins <- seq(0, 1, by=binsize)

        for(j in 1:(length(wins)-1)){
                bin <- mean(prop[which(prop$coord.prop >= wins[j] & prop$coord.prop < wins[j+1]),][[target_column]])
                collect_bin <- c(collect_bin, bin)
        }

        collect_bin.smooth <- smooth.spline(wins[2:length(wins)], collect_bin, nknots=n_knots)$y

        # k=mafilt.size
        # filt <- 1/(2*k+1)
        # filt <- rep(filt, 2*k+1)
        # filt.dat <- stats::filter(collect_bin, filt)
        # filt.dat <- loess(formula = paste(target_column, "coord.prop", sep="~"), span=0.1, degree=1, data=prop)
        # filt.dat <- filt.dat$fitted

        res <- tibble(prop=wins,
                        bin=collect_bin,
                        bin.smooth <- collect_bin.smooth
                    #     bin=collect_bin,
                    #   bin.ma=filt.dat,
                      )
        
        return(res)
}

pTelCen.meth.ma9 <- ggplot() + 
        geom_line(data=filter(prop_co, !(group %in% c("snp", "mC"))), aes(x=prop, y=bin, colour=group), size=0.4) +
        geom_area(data=filter(prop_co, group=="mC"), aes(x=prop, y=bin*0.1), fill="yellowgreen", alpha=0.5) +
        geom_hline(data=mean_cos, aes(yintercept=meanco, colour=group), linetype="dashed", size=0.3) +
        scale_y_continuous(name="Crossovers per F2", sec.axis=sec_axis(~.*10, name="DNA methylation (mC/C)")) +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        scale_colour_manual(values=pal) +
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=7),
              legend.text=element_text(size=7),
              legend.position="top") +
        theme(text=element_text(size=9, family="helvetica", colour="black"),
        axis.text=element_text(size=7, colour="black"))

        
#svg(file=paste0(prefix, "_snp_telcen-ma9.svg"), width=3.25, height=2.25)
#pTelCen.ma9
#dev.off()
svg(file=paste0(prefix, "_meth_telcen-ma9.svg"), width=3.25, height=2.25)
pTelCen.meth.ma9
dev.off()