# Description 
library(tidyverse)
library(khroma)
library(scales)
library(extrafont)
library(ggsci)
library(ggbeeswarm)
font_import(pattern="Arial", prompt=FALSE)
loadfonts(device="pdf")

# INPUT DATA =============================================================
wd <- "/home/namilhand/01_Projects/PGR_informatics/Lesson05_various_CO_analysis/2_Col-CEN_landscape"
setwd(wd)
dir.create(file.path(wd, "plots"), recursive=T)

dir_gbs <- file.path(wd, "data/crossovers/genomeBin")

## GBS landscape
wt_i <- read_tsv(file.path(dir_gbs, "Selected_inhouse_WT_Col-CEN_genomeBin100kb.tsv"))
wt_r <- read_tsv(file.path(dir_gbs, "Selected_Rowan_WT_Col-CEN_genomeBin100kb.tsv"))

mj3_1 <- read_tsv(file.path(dir_gbs, "pJ3-mJ3_set1_Col-CEN_genomeBin100kb.tsv"))
mj3_2 <- read_tsv(file.path(dir_gbs, "pJ3-mJ3_set2_Col-CEN_genomeBin100kb.tsv"))
mj3_3 <- read_tsv(file.path(dir_gbs, "pJ3-mJ3_set3_Col-CEN_genomeBin100kb.tsv"))

h2aw_1 <- read_tsv(file.path(dir_gbs, "h2aw-null_set1_Col-CEN_genomeBin100kb.tsv"))
h2aw_2 <- read_tsv(file.path(dir_gbs, "h2aw-null_set2_Col-CEN_genomeBin100kb.tsv"))
h2aw_3 <- read_tsv(file.path(dir_gbs, "h2aw-null_set3_Col-CEN_genomeBin100kb.tsv"))

## merge GBS results from the same genotype
wt_tsv <- dplyr::select(wt_i, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = wt_i$coInWindow + wt_r$coInWindow) %>%
    add_column(libSize = wt_i$libSize + wt_r$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

mj3_tsv <- dplyr::select(mj3_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = mj3_1$coInWindow + mj3_2$coInWindow + mj3_3$coInWindow) %>%
    add_column(libSize = mj3_1$libSize + mj3_2$libSize + mj3_3$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

h2aw_tsv <- dplyr::select(h2aw_1, -c(coInWindow, libSize)) %>%
    add_column(coInWindow = h2aw_1$coInWindow + h2aw_2$coInWindow + h2aw_3$coInWindow) %>%
    add_column(libSize = h2aw_1$libSize + h2aw_2$libSize + h2aw_3$libSize) %>%
    mutate(mean.coInWindow = coInWindow/libSize)

co_tsv.list <- list(wt=wt_tsv, mj3=mj3_tsv, h2aw=h2aw_tsv)



#========= Setting for plots =========
# Figure theme setting {
my_theme <- theme_classic() +
    theme(text = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.text = element_text(size=6, family="Arial", colour="black")) +
    theme(plot.title = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.title = element_text(size=6, family="Arial", colour="black")) +
    theme(legend.margin = margin(0.5,0.5,0.5,0.5)) +
    theme(strip.text = element_text(size=6, family="Arial", colour="black")) +
    theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))
    # theme(legend.spacing.y = unit(1, "cm"))
theme_set(my_theme)

pal_GBS <- c("black", pal_d3(palette = c("category10"))(2))
names(pal_GBS) <- c("wt","mj3", "h2aw")
# show_col(pal_GBS)

#====== Col-CEN genome setting ======

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

pericen_left <- c(11442327, 1029534, 10397911, 1167024, 8901152)
pericen_right <- c(20412241, 9828983, 19004143, 9609618, 18023365)

nrz_mid <- (nrz_left + nrz_right)/2
cen_mid <- (cen_left + cen_right)/2

col_cen.coord <- col_cen %>%
    add_column(nrz_left, nrz_right) %>%
    add_column(nrz_left.cum = nrz_left + tha.cum[1:5]) %>%
    add_column(nrz_right.cum = nrz_right + tha.cum[1:5]) %>%
    add_column(lrz_left.cum = lrz_left + tha.cum[1:5]) %>%
    add_column(lrz_right.cum = lrz_right + tha.cum[1:5]) %>%
    add_column(cen_left.cum = cen_left + tha.cum[1:5]) %>%
    add_column(cen_right.cum = cen_right + tha.cum[1:5]) %>%
    add_column(pericen_left.cum = pericen_left + tha.cum[1:5]) %>%
    add_column(pericen_right.cum = pericen_right + tha.cum[1:5])

## cumulative coordinates of pericentromere
coord_pericen <- tibble(chr=1:5, north=pericen_left + tha.cum[1:5], south=pericen_right + tha.cum[1:5])

# SV
syri_out <- read_tsv(file.path(wd, "data/colcen/syri.out"), col_names=F)

collect_sv <- function(sv_type){
    # sv_type <- "NOTAL"
    result <- filter(syri_out, X11 == sv_type) %>%
        dplyr::select(c("X1", "X2", "X3", "X9", "X10", "X11"))

    colnames(result) <- c("chr", "start", "stop", "uniqID", "parentID", "type")
    result <- filter(result, start != "-") %>%
        mutate(start = as.numeric(start), stop = as.numeric(stop)) %>%
        mutate(width = stop - start + 1, .after="stop") %>%
        mutate(cum.start = start + tha.cum[as.numeric(gsub("Chr", "", chr))]) %>%
        mutate(cum.stop = stop + tha.cum[as.numeric(gsub("Chr", "", chr))])

    #print(sum(result$width)/1000)
    return(result)
}
notal <- collect_sv("NOTAL")
# 20.65 Mb
trans <- collect_sv("TRANS")
# 0.45 Mb
inv <- collect_sv("INV")
# 2.51 Mb
invtr <- collect_sv("INVTR")
# 0.02 Mb

#===== MA smoothing ====
ma_width.co <- 7

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
co_ma <- mutate(co_ma, genotype = factor(genotype, levels=c("wt", "mj3", "h2aw")))  %>%
    mutate(cMMb = smooth *10 * 100/2)

#===== DRAWING CO LANDSCAPE =====

ylim_max <- max(co_ma$cMMb) * 1.1
landscape_ylim <- c(0, ylim_max)

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
        # geom_rect(data=dat_pval , aes(xmin = cum.start, xmax = cum.end, ymin=landscape_ylim[2]*0, ymax=landscape_ylim[2]*0.95), fill="#ffff6b", alpha=0.7) +
        geom_rect(data=notal, aes(xmin=cum.start, xmax=cum.stop, ymin=landscape_ylim[2]* 0.95, ymax=landscape_ylim[2]), fill="cyan", colour=NA) +
        geom_rect(data=trans, aes(xmin=cum.start, xmax=cum.stop, ymin=landscape_ylim[2] * 0.95, ymax=landscape_ylim[2]), fill="blue", colour=NA) +
        geom_rect(data=inv, aes(xmin=cum.start, xmax=cum.stop, ymin=landscape_ylim[2] * 0.95, ymax=landscape_ylim[2]), fill="hotpink", colour=NA) +
        geom_rect(data=col_cen.coord, aes(xmin=cen_left.cum, xmax=cen_right.cum, ymin=landscape_ylim[2] * 0.95, ymax=landscape_ylim[2]), fill="grey30", colour=NA) +
        geom_line(data=dat, aes(x=cumwindow, y=cMMb, colour=genotype), linewidth=0.2) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(tha.cum), 20*10^6), 132000000), expand=c(0.01, 0.01)) +
        scale_y_continuous(limits=landscape_ylim) +
        scale_colour_manual(values = pal) +
        geom_vline(xintercept=tha.cum[2:5], colour="Black", linewidth=0.2) +
        # geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", linewidth=0.2) +
        labs(y="cM/Mb") +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
                legend.title=element_text(size=6),
                legend.text=element_text(size=6)) +
        theme(legend.title=element_blank()) +
        # theme(legend.position = c(1, 1)) +
        theme(legend.position = "none") +
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
drawCoLandscape(co_ma, pal_GBS, "plots/co_lancscape", 4.5, 1.3)

drawCoDiffLandscape <- function(dat, gt_test, gt_ctrl, linecolour, y_label, prefix, plot_width, plot_height, y_min, y_max){

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
        geom_rect(data=notal, aes(xmin=cum.start, xmax=cum.stop, ymin=diff_landscape_ylim[2]* 0.95, ymax=diff_landscape_ylim[2]), fill="cyan", colour=NA) +
        geom_rect(data=trans, aes(xmin=cum.start, xmax=cum.stop, ymin=diff_landscape_ylim[2] * 0.95, ymax=diff_landscape_ylim[2]), fill="blue", colour=NA) +
        geom_rect(data=inv, aes(xmin=cum.start, xmax=cum.stop, ymin=diff_landscape_ylim[2] * 0.95, ymax=diff_landscape_ylim[2]), fill="hotpink", colour=NA) +
        geom_rect(data=col_cen.coord, aes(xmin=cen_left.cum, xmax=cen_right.cum, ymin=diff_landscape_ylim[2] * 0.95, ymax=diff_landscape_ylim[2]), fill="grey30", colour=NA) +
        geom_line(data=dat_diff_bind, aes(x=cumwindow, y=difference, colour=genotype), linewidth=0.2) +
        scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(tha.cum), 20*10^6), 132000000), expand=c(0.01, 0.01)) +
        scale_y_continuous(limits=diff_landscape_ylim) +
        scale_colour_manual(values=linecolour) +
        geom_vline(xintercept=tha.cum[2:5], colour="Black", linewidth=0.2) +
        # geom_vline(xintercept=centromeres.cum, colour="Black", linetype="dashed", linewidth=0.2) +
        geom_hline(yintercept=0, linewidth=0.2, colour="Black") +
        labs(y=y_label) +
        theme_classic() +
        theme(legend.key.size=unit(0.2, "inches"),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6)) +
        theme(legend.title=element_blank()) +
        # theme(legend.position = c(1, 1)) +
        theme(legend.position = "none") +
        theme(legend.justification = c(1,0.8)) +
        theme(legend.box.margin = margin(0.5,0.5,0.5,0.5)) +
        theme(text=element_text(size=6, family="Arial", colour="black"),
        axis.text=element_text(colour="black", size=6)) +
        theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))

        pdf(file=paste0(prefix, "_co-diff_landscape_ma-co_", ma_width.co,".pdf"), width=plot_width, height=plot_height)
        print(p)
        dev.off()
        png(file=paste0(prefix, "_co-diff_landscape_ma-co_", ma_width.co,".png"), width=plot_width, height=plot_height, unit="in", res=300)
        print(p)
        dev.off()
}

drawCoDiffLandscape(co_ma, c("h2aw", "mj3"), "wt", pal_GBS, "-WT", "plots/all-wt", 4.5, 1.3, -3, 9)


#' 7 TEL-CEN plotting
# distFromTel() : convert coordinate into proportion of distance from cent to tel
distFromTel <- function(dat, binsize, spline_df, target_column){
		# dat <- tsv_male.list2$wt_gbs
        # spline_df <- 11
        # target_column <- "relco"
		# binsize <- 0.01
        # co_telcen.list <- lapply(co_tsv.list, distFromTel, 0.01, spline_df, "mean.coInWindow")

        dat <- arrange(dat, chr, window)
        dat$cen_mid <- rep(cen_mid, times=table(dat$chr))
        dat$end <- rep(chr.ends, times=table(dat$chr))
        
        left <- filter(dat, window < cen_mid) %>%
                add_column(arm="north") %>%
                mutate(coord.prop=window/cen_mid) %>%
				filter(chr %in% c("Chr1", "Chr3", "Chr5"))

		# Excluded north arm of Chr2 and Chr4 as these regions is too short for fair comparison with other chromosomes
        right <- filter(dat, window >= cen_mid) %>%
                add_column(arm="south") %>%
                mutate(coord.prop=(end-window)/(end-cen_mid))
        
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

        res <- tibble(prop=wins[2: length(wins)],
                        bin=collect_bin,
                        bin.smooth=collect_bin.smooth
                      )
        
        return(res)
}

telcen_pericen_boarder <- tibble(chr = 1:5, north = pericen_left, south = pericen_right) %>%
    add_column(end=chr.ends) %>%
    add_column(cen_mid=cen_mid) %>%
    pivot_longer(cols=c("north", "south"), names_to="region", values_to="pericen_boarder") %>%
    filter(!(chr == 2 & region == "north")) %>%
    filter(!(chr == 4 & region == "north")) %>%
    mutate(pericen_boarder_prop = case_when(
        region == "north" ~ pericen_boarder/cen_mid,
        region == "south" ~ (end - pericen_boarder)/(end - cen_mid),
        TRUE ~ NA))

# spline_df <- 11 
spline_df <- 20
co_telcen.list <- lapply(co_tsv.list, distFromTel, 0.01, spline_df, "mean.coInWindow")
co_telcen.bind <- bind_rows(co_telcen.list, .id="sample") %>%
    mutate(sample = factor(sample, levels=c("wt", "mj3", "h2aw")))


## define centromeric-pericentromeric region
mean_telcen_pericen_border <- mean(telcen_pericen_boarder$pericen_boarder_prop)

drawTelCenProfile <- function(dat_co, prefix, plot_width, plot_height){

    # dat_co <- co_telcen.bind
    # chrom_label <- "cM/Mb"

    p <- ggplot() +
        geom_line(data=dat_co, aes(x=prop, y=bin.smooth * 100 / 2, colour=sample), linewidth=0.2) +
        geom_vline(xintercept=mean_telcen_pericen_border, linewidth=0.2, linetype = "dashed") +
        scale_colour_manual(values=pal_GBS) +
        scale_y_continuous(name = "cM/Mb") +
        scale_x_continuous(name="Distance from the telomere (ratio)",
                           breaks=seq(0, 1, by=0.25),
                           label=c("TEL", 0.25, 0.5, 0.75, "CEN"))+
        theme_classic() +
        theme(legend.key.size=unit(0.1, "inches"),
              legend.title=element_text(size=6),
              legend.text=element_text(size=6),
              legend.position=c(0,1),
              legend.justification=c(-0.2, 0.8)) +
        theme(legend.title = element_blank()) +
        theme(legend.position = "none") +
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

drawTelCenProfile(co_telcen.bind, "WT_mJ3_h2aw", 2.1, 1.4)

drawTelCenDiffProfile <- function(dat_co, ctrl_group, prefix, plot_width, plot_height){

    dat_co <- mutate(dat_co, bin.cmmb= bin.smooth * 10 / 2 * 100)

    dat_co.ctrl <- filter(dat_co, sample == ctrl_group)
    dat_co.diff <- filter(dat_co, sample != ctrl_group) %>%
        add_column(bin.cmmb.wt = rep(dat_co.ctrl$bin.cmmb, times=length(unique(.$sample)))) %>%
        mutate(bin.cmmb.diff = bin.cmmb - bin.cmmb.wt)

    pal <- pal_GBS

    yscale_lim <- c(-2, 27.5)

    if(ctrl_group == "wt"){
        ctrl_group = "WT"
    }
    y_label <- paste0("-", ctrl_group)
    y_label <- bquote(Delta ~ " cM/Mb (" ~ .(y_label) ~ ")")



    p <- ggplot() +
        # annotate("rect", xmin=telcen_pericen_border, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey90") +
        geom_vline(xintercept=mean_telcen_pericen_border, linewidth=0.2, linetype = "dashed", colour="black") +
        geom_line(data=dat_co.diff, aes(x=prop, y=bin.cmmb.diff, colour=sample), size=0.4, linewidth=0.2) +
        geom_hline(yintercept=0, linewidth=0.2, colour="Black") +
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
drawTelCenDiffProfile(co_telcen.bind, "wt", "ctrl=wt", 1.8, 1.4)
