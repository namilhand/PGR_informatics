# Description 
library(tidyverse)
# library(khroma)
library(scales)
library(extrafont)
library(ggsci)
# library(RColorBrewer)
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

coord_subpericen <- read_csv("data/tair10_cumulative_coord.csv")

# INPUT DATA =============================================================

# wd <- /your/working/directory/here"
wd <- "/home/namilhand/01_Projects/PGR_informatics/GBS_analysis/genomewide_profile"
setwd(wd)
dir.create(file.path(wd, "results/plots"), recursive=T)

dir_gbs <- "data/tsv"

## GBS landscape
wt <- read_tsv(file.path(dir_gbs, "Selected_inhouse_WT_v2_genomeBin100kb.tsv"))
h2aw_1 <- read_tsv(file.path(dir_gbs, "h2aw-null_set1_genomeBin100kb.tsv"))
h2aw_2 <- read_tsv(file.path(dir_gbs, "h2aw-null_set2_genomeBin100kb.tsv"))
h2aw_3 <- read_tsv(file.path(dir_gbs, "h2aw-null_set3_genomeBin100kb.tsv"))

## dname landscape
mC_col <- read_tsv(file.path("data", "WT_mC_TAIR10_window100kb_step100kb.bedg"), col_names=c("chr", "window", "window_end", "mC.ratio"))
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
wt_tsv <- mutate(wt, mean.coInWindow = coInWindow/libSize)

h2aw_tsv <- dplyr::select(h2aw_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = h2aw_1$coInWindow + h2aw_2$coInWindow + h2aw_3$coInWindow) %>%
    add_column(libSize = 287) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

co_tsv.list <- list(wt=wt_tsv, h2aw=h2aw_tsv)


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
co_ma <- mutate(co_ma, genotype = factor(genotype, levels=c("wt", "h2aw")))  %>%
    mutate(cMMb = smooth *10 * 100/2)
    
#===== DRAWING CO LANDSCAPE =====
## palette
pal_GBS <- c("blue", "red")
names(pal_GBS) <- c("wt", "h2aw")

# co landscape

## y-axis scale
landscape_ylim <- c(0, max(co_ma$cMMb)*1.1) 

drawCoLandscape <- function(dat, pal, prefix, plot_width, plot_height){

    ## convert CO per F2 per 100 kb into cM/Mb
    ## *10 : CO per F2 per 100kb --> CO per F2 per Mb
    ## *100: CO per F2 per Mb (M/Mb) --> cM/Mb

    p <- ggplot() +
        annotate("rect", xmin=coord_pericen$north[1], xmax=coord_pericen$south[1], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[2], xmax=coord_pericen$south[2], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[3], xmax=coord_pericen$south[3], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[4], xmax=coord_pericen$south[4], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        annotate("rect", xmin=coord_pericen$north[5], xmax=coord_pericen$south[5], ymin=-Inf, ymax=Inf, fill="#D8EDFF") +
        geom_line(data=dat, aes(x=cumwindow, y=cMMb, colour=genotype), linewidth=0.2) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(co_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
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

        return(p)
}

gbs_landscape <- drawCoLandscape(co_ma, pal_GBS)

pdf(file=file.path("results", "co_landscape.pdf"), width=3.3, height=1.0)
print(gbs_landscape)
dev.off()

png(file=file.path("results", "co_landscape.png"), width=3.3, height=1.0, unit="in", res=300)
print(gbs_landscape)
dev.off()


drawCoDiffLandscape <- function(dat, gt_test, gt_ctrl, linecolour, y_label, y_min, y_max){

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
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(co_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
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

        return(p)

}

## adjust y_min and y_max accordingly
diff_landscape <- drawCoDiffLandscape(co_ma, "h2aw", "wt", pal_GBS, "-WT", -4, 4.5)

pdf(file="results/diff_landscape.pdf", width=3.3, height=1)
print(diff_landscape)
dev.off()

png(file="results/diff_landscape.png", width=3.3, height=1, unit="in", res=300)
print(diff_landscape)
dev.off()

#========= TEL-CEN plotting ==============
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

        collect_bin.smooth <- smooth.spline(wins[2:length(wins)], collect_bin, df=spline_df)$y
        collect_bin.cmmb <- collect_bin.smooth/2*10*100

        res <- tibble(prop=wins[2: length(wins)],
                        bin=collect_bin,
                        bin.smooth=collect_bin.smooth,
                        bin.cmmb=collect_bin.cmmb
                      )
        
        return(res)
}

spline_df <- 11
co_telcen.list <- lapply(co_tsv.list, distFromTel, 0.01, spline_df, "mean.coInWindow")
co_telcen.bind <- bind_rows(co_telcen.list, .id="sample") %>%
    mutate(sample = factor(sample, levels=c("wt", "h2aw")))
mC_col_telcen <- distFromTel(mC_col, 0.01, spline_df, "mC.ratio")


## define centromeric-pericentromeric region
## region where DNA methylation level is higher than the average
define_telcen_pericen <- mC_col_telcen %>%
    mutate(distFromAvg = abs(bin.cmmb - mean(bin.cmmb))) %>%
    arrange(distFromAvg)
telcen_pericen_border <- define_telcen_pericen$prop[1]
## 0.76

drawTelCenProfile <- function(dat_co, dat_chromatin, chrom_label, prefix){
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
    
}


telcen_profile <- drawTelCenProfile(co_telcen.bind, mC_col_telcen, "DNA methylation (ratio)")

pdf(file="results/telcen_profile.pdf", width=2.4, height=1.4)
print(telcen_profile)
dev.off()
png(file="results/telcen_profile.png", width=2.4, height=1.4, unit="in", res=300)
print(telcen_profile)
dev.off()

drawTelCenDiffProfile <- function(dat_co, ctrl_group, prefix){
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
    
}

diffTelCen <- drawTelCenDiffProfile(co_telcen.bind, "wt", "ctrl=wt")

pdf(file="results/telcen-diff_profile.pdf", width=2.2, height=2.08)
print(diffTelCen)
dev.off()

png(file="results/telcen-diff_profile.png", width=2.2, height=2.08, unit="in", res=300)
print(diffTelCen)
dev.off()