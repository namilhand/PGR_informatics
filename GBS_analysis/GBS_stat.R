# Description 
library(tidyverse)
# library(khroma)
library(scales)
library(extrafont)
# library(RColorBrewer)
library(ggsci)
library(multcompView)
library(rstatix)
font_import(pattern="Arial", prompt=F)
loadfonts(device="pdf")

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
theme_set(my_theme)

pal_GBS <- c("blue", "red")
names(pal_GBS) <- c("wt","h2aw")

# ===== INPUT DATA =====
wd <- "/home/namilhand/01_Projects/PGR_informatics/GBS_analysis"
setwd(wd)
dir.create(file.path(wd, "example/results"), recursive=T)

dir_cotable <- "example/data/cotable"

## cotable
wt <- read_csv(file.path(dir_cotable, "Selected_inhouse_WT_cotable.txt"))
h2aw_1 <- read_csv(file.path(dir_cotable, "h2aw-null_set1_cotable.txt"))
h2aw_2 <- read_csv(file.path(dir_cotable, "h2aw-null_set2_cotable.txt"))
h2aw_3 <- read_csv(file.path(dir_cotable, "h2aw-null_set3_cotable.txt"))


cotable_list <- list("wt" = wt, "h2aw_1"=h2aw_1, "h2aw_2"=h2aw_2, "h2aw_3"=h2aw_3)

# ===== GENOME INFO SETTING =====
#' 0. chromosome info 
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)

tha.cum <- c(0, cumsum(chr.ends))
tha.tot <- tha.cum[length(tha.cum)]
centromeres <- centromeres + tha.cum[1:5]
# pericentromeric + centromeric region
# (Ref.: Underwood, 2018, Genome Research)
north.start <- c(11420001, 910001, 10390001, 1070001, 8890001)
south.end <- c(18270000, 7320000, 16730000, 6630000, 15550000)

# cumulative coordinates of pericentromere
coord_pericen <- tibble(chr=1:5, north=north.start + tha.cum[1:5], south=south.end + tha.cum[1:5])

# Read arabidopsis cumulative genome coord
tha.cum_coord <- read_csv("example/data/tair10_cumulative_coord.csv")


# ===== INPUT PROCESSING =====
## trimming library name
trimLibName <- function(dat){
    output <- mutate(dat, lib = str_replace(lib, "results/06_tiger/lib", "")) %>%
            mutate(lib = str_replace(lib, "_MappedOn_tair10", "")) %>%
            mutate(lib = as.numeric(lib)) %>%
            arrange(lib, chrs, start)
    return(output)
}



cotable_list.libname_trimmed <- lapply(cotable_list, trimLibName)

## bind cotables from the same genotype
cotable_wt <- cotable_list.libname_trimmed[["wt"]]

cotable_h2aw <- cotable_list.libname_trimmed[c("h2aw_1", "h2aw_2", "h2aw_3")]
cotable_h2aw[["h2aw_2"]]$lib = cotable_h2aw[["h2aw_2"]]$lib + max(unique(cotable_h2aw[["h2aw_1"]]$lib))
cotable_h2aw[["h2aw_3"]]$lib = cotable_h2aw[["h2aw_3"]]$lib + max(unique(cotable_h2aw[["h2aw_2"]]$lib))
cotable_h2aw <- bind_rows(cotable_h2aw)

## merge the cotables from different genotypes
cotable_all.list <- list("wt"=cotable_wt, "h2aw"=cotable_h2aw)
cotable_all <- bind_rows(cotable_all.list, .id="genotype") %>%
    mutate(genotype = factor(genotype, levels=c("wt", "h2aw")))

## INFO: no.of libraries for each genotype
lib_genotypes <- c("wt", "h2aw")
lib_count <- c(240, 288)
lib_meta <- tibble(genotype = factor(lib_genotypes, levels=lib_genotypes), nlib = lib_count)

## count the number of CO for each individual
cos.count.backbone <- tibble(genotype=factor(levels=lib_genotypes), lib=numeric())
for(i in 1:nrow(lib_meta)){
    # i <- 1
    cos.count.backbone.gt <- rep(lib_meta$genotype[i], times=lib_meta$nlib[i])
    cos.count.backbone.lib <- 1:lib_meta$nlib[i]
    cos.count.backbone.temp <- tibble(genotype=cos.count.backbone.gt, lib=cos.count.backbone.lib)
    cos.count.backbone <- bind_rows(cos.count.backbone, cos.count.backbone.temp)
}
cos.count.pre <- cotable_all %>%
        group_by(genotype) %>%
        count(lib)
cos.count <- left_join(cos.count.backbone, cos.count.pre, by=c("genotype", "lib")) %>%
    mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n))

## filtering F2
### Note: h2aw includes 2 F2 library with no CO site: 3, 136
### Note: lib 106 of h2aw set (lib10 of Set2) is outlier as a result of low sequencing depth. Filter out.

cos.count <- filter(cos.count, n!=0) %>%
    bind_rows(tibble(genotype=factor("h2aw", "h2aw"), lib=c(3, 136), n=c(0,0))) %>%
    filter(!(genotype == "h2aw" & lib == 106)) %>%
    arrange(genotype, lib)

### update lib_meta
lib_meta$nlib <- table(cos.count$genotype)

# ===== CO METAPROFILE: CO MEAN, CO WIDTH, MWU, T-TEST =====
## wilcox test
two_group_wilcox <- function(test_group, ctrl_group){
    # test_group <- "dmigs"
    # ctrl_group <- "wt"
    wilcox_dat <- filter(cos.count, genotype %in% c(ctrl_group, test_group)) %>%
        dplyr::select(c(genotype, n))
    wilcox_pvalue <- wilcox.test(n ~ genotype, data=wilcox_dat)$p.value
    return(wilcox_pvalue)
}
wilcox_pvalue.h2aw <- two_group_wilcox("h2aw", "wt")

wilcox_pvalue <- c(1, wilcox_pvalue.h2aw)

## t test
t_test.pvalue <- c(1)
for(gt in c("h2aw")){ 
    pval <- t.test(filter(cos.count, genotype == "wt")$n, filter(cos.count, genotype == gt)$n)$p.value
    t_test.pvalue <- c(t_test.pvalue, pval)
}

## CO summary
cos.width <- cotable_all %>%
    group_by(genotype) %>%
    summarise(mean.width = mean(width))
cos.width <- cos.width$mean.width

cos.summary <- cos.count %>%
    group_by(genotype) %>%
    summarise(co_total = sum(n), nlib=length(unique(lib)), co_mean = mean(n), co_sd = sd(n)) %>%
    add_column(wilcox.vs_wt=wilcox_pvalue) %>%
    add_column(t_test.vs_wt=t_test.pvalue) %>%
    add_column(cos.width)

write_csv(cos.summary, file="CO_meta_summary.csv")


## CO frequency histogram
p.co_freq <- ggplot() +
    geom_histogram(data = cos.count, aes(x=n), binwidth=1, fill="white", colour="black", linewidth=0.2) +
    geom_vline(data=cos.summary, aes(xintercept=co_mean), linewidth=0.4, linetype="dashed", colour="red") +
    facet_grid(rows="genotype", scales="free_y") +
    labs(y="Frequency", x="Crossovers")

dir.create("plots", recursive=T)
pdf(file=file.path("plots", "CO_frequency_hist.pdf"), width=1.8, height=1.5)
print(p.co_freq)
dev.off()
png(file=file.path("plots", "CO_frequency_hist.png"), width=1.8, height=1.5, unit="in", res=300)
print(p.co_freq)
dev.off()

# == MEAN CO BY ARM AND PERICEN ==
cotable_armperi <- cotable_all %>%
    add_column(region = "dummy")

## annotate each CO site by the arm/pericen region
for(i in 1:5){
    # i <- 1
    chr_n <- paste0("Chr", i)
    # north <- coord_pericen$north[i]
    # south <- coord_pericen$south[i]
    # north <- north.start[i]
    # south <- south.end[i]

    north <- tha.cum_coord$pericen.start[i] - tha.cum_coord$chr.start[i]
    ColdPeri.start <- tha.cum_coord$ColdPeri.start[i] - tha.cum_coord$chr.start[i]
    south <- tha.cum_coord$pericen.end[i] - tha.cum_coord$chr.start[i]
    ColdPeri.end <- tha.cum_coord$ColdPeri.end[i] - tha.cum_coord$chr.start[i]

    # cotable_armperi <- cotable_armperi %>%
    #     mutate(region = case_when((chrs == chr_n) & (cos >= north) & (cos <= south) ~ "peri",
    #                                 chrs == chr_n ~ "arm",
    #                                 TRUE ~ region))
    cotable_armperi <- cotable_armperi %>%
        mutate(region = case_when((chrs == chr_n) & (cos >= north) & (cos <= south) ~ "peri",
                                    (chrs == chr_n) ~ "arm",
                                    TRUE ~ region))
    
}
cotable_armperi <- mutate(cotable_armperi, region = factor(region, levels=c("arm", "peri"))) %>%
    arrange(genotype, lib)

## summarise arm/peri CO count of each library

## initialize cotable_armperi_count
cotable_armperi_count.backbone <- tibble(genotype = factor(levels=lib_genotypes), lib=numeric(), region=factor(levels=c("arm", "peri")))
for(i in 1:nrow(lib_meta)){
    # i <- 1
    region_levels_length <- length(levels(cotable_armperi_count.backbone$region))
    selected_genotype <- lib_meta$genotype[i]

    cotable_armperi_count.backbone.gt <- rep(lib_meta$genotype[i], times=lib_meta$nlib[i]*region_levels_length)
    cotable_armperi_count.backbone.lib <- rep(filter(cos.count, genotype == selected_genotype)$lib, each=region_levels_length)
    cotable_armperi_count.backbone.region <- factor(rep(levels(cotable_armperi_count.backbone$region), times=lib_meta$nlib[i]), levels=levels(cotable_armperi_count.backbone$region))
    cotable_armperi_count.backbone.temp <- tibble(genotype=cotable_armperi_count.backbone.gt, lib=cotable_armperi_count.backbone.lib, region=cotable_armperi_count.backbone.region)

    cotable_armperi_count.backbone <- bind_rows(cotable_armperi_count.backbone, cotable_armperi_count.backbone.temp)
}

cotable_armperi_count1 <- cotable_armperi %>%
    group_by(genotype, lib, region) %>%
    count(genotype, lib, region, .drop=FALSE) %>%
    ungroup()

cotable_armperi_count <- left_join(cotable_armperi_count.backbone, bind_rows(cotable_armperi_count1), by=c("genotype", "lib", "region")) %>%
    mutate(n = case_when(is.na(n) ~ 0, TRUE ~ n))


## summarise arm/peri CO count of each genotype
cotable_armperi_summary <- cotable_armperi_count %>%
    group_by(genotype, region) %>%
    summarise(nlib = length(unique(lib)), co_mean=mean(n), co_sd=sd(n), co_se=sd(n)/sqrt(length(n))) %>%
    arrange(region, genotype)

## STAT test (t-test)
armperi_ttest <- function(ctrl_group, test_group, roi, althyp){
    # ctrl_group <- "wt"
    # test_group <- "mmH"
    # roi <- "ColdPeri"
    # althyp <- "greater"
    # roi = region of interest
    t.test(filter(cotable_armperi_count, genotype == test_group & region == roi)$n, filter(cotable_armperi_count, genotype == ctrl_group & region == roi)$n, alternative = althyp)$p.value
}
armperi_wilcox <- function(ctrl_group, test_group, roi, althyp){
    # ctrl_group <- "wt"
    # test_group <- "mmH"
    # roi <- "ColdPeri"
    # althyp <- "two.sided"
    # roi = region of interest
    wilcox.test(filter(cotable_armperi_count, genotype == test_group & region == roi)$n, filter(cotable_armperi_count, genotype == ctrl_group & region == roi)$n, alternative = althyp)$p.value
}

armperi_pval_twoside <- c()
armperi_wilcox_twoside <- c()
# armperi_pval_oneside <- c()
for(i in 1:nrow(cotable_armperi_summary)){
    ctrl_group <- "wt"
    test_group <- cotable_armperi_summary$genotype[i]
    roi <- cotable_armperi_summary$region[i]

    ttest_pval <- armperi_ttest(ctrl_group, test_group, roi, "two.sided")
    wilcox_pval <- armperi_wilcox(ctrl_group, test_group, roi, "two.sided")
    # oneside_pval <- armperi_ttest(ctrl_group, test_group, roi, "greater")

    armperi_pval_twoside <- c(armperi_pval_twoside, ttest_pval)
    armperi_wilcox_twoside <- c(armperi_wilcox_twoside, wilcox_pval)
    # armperi_pval_oneside <- c(armperi_pval_oneside, oneside_pval)
}
## update cotable_armperi_summary with pvalue

cotable_armperi_summary <- cotable_armperi_summary %>%
    add_column(pval.ttest = armperi_pval_twoside) %>%
    add_column(pval.wilcox = armperi_wilcox_twoside)

write_tsv(cotable_armperi_summary, file="ArmPeri_CO-rate_mean_sd_se_ttest.tsv")

## STAT test (chi-sqaure)
## chisq test (proportion of pericentromeric CO)
cotable_armperi_count_sum <- cotable_armperi_count %>%
    group_by(genotype, region) %>%
    summarise(n=sum(n), nlib=n()) %>%
    mutate(prop = n/sum(n))

cotable_armperi_count_sum_wide <- cotable_armperi_count_sum  %>% 
    pivot_wider(names_from=region, values_from=c(n, prop)) %>%
    mutate(mean_co_arm = n_arm/nlib, mean_co_peri = n_peri/nlib)

chisqInRegion <- function(region){
    chisq_data_all <- dplyr::select(cotable_armperi_count_sum_wide, c("genotype", "n_arm", "n_peri")) %>%
        mutate(total = sum(n_arm, n_peri))
    colnames(chisq_data_all) <- c("genotype", "arm", "peri", "total")
    # region = "arm"
    ctrl_gt = "wt"
    pvalues <- c()
    for(i in 1:nrow(cotable_armperi_count_sum_wide)){
        # i <- 2
        test_gt = as.character(cotable_armperi_count_sum_wide$genotype[i])
        chisq_table <- tibble(genotype = chisq_data_all$genotype, test_region = chisq_data_all[[region]], outside = chisq_data_all[["total"]] - chisq_data_all[[region]]) %>%
            filter(genotype %in% c(ctrl_gt, test_gt)) %>%
            column_to_rownames(var = "genotype")
        pvalues <- c(pvalues, chisq.test(chisq_table)$p.value)
    }
    return(pvalues)
}
cotable_armperi_count_sum_wide$chisq.arm = chisqInRegion("arm")
cotable_armperi_count_sum_wide$chisq.peri = chisqInRegion("peri")

write_tsv(cotable_armperi_count_sum_wide, file="hta_co_summary_with_chisq.tsv")

###################
# ANOVA: COs compared in pericentromere/arm
###################

cotable_armperi_count.peri <- filter(cotable_armperi_count, region == "peri") %>%
    dplyr::select(c(genotype, n))
cotable_armperi_count.arm <- filter(cotable_armperi_count, region == "arm") %>%
    dplyr::select(c(genotype, n))

# Test for homogeneity of variance
levene.peri <- levene_test(cotable_armperi_count.peri, n ~ genotype)
levene.arm <- levene_test(cotable_armperi_count.arm, n ~ genotype)
write.csv(file="cotable_armperi_levene_test.csv", bind_rows(list(peri=levene.peri, arm=levene.arm), .id="region"))

# pericentromere
## welch's one-way anova test
peri.aov <- welch_anova_test(n ~ genotype, data=cotable_armperi_count.peri)
## Pairwise comparisons (Games-Howell)
peri.gh <- games_howell_test(cotable_armperi_count.peri, n ~ genotype)
peri.result <- peri.gh$p.adj <= .05
names(peri.result) <- str_c(peri.gh$group1, peri.gh$group2, sep="-")
peri.letter <- multcompLetters(peri.result)

# arm
## welch's one-way anova test
arm.aov <- welch_anova_test(n ~ genotype, data=cotable_armperi_count.arm)
## Pairwise comparisons (Games-Howell)
arm.gh <- games_howell_test(cotable_armperi_count.arm, n ~ genotype)
arm.result <- arm.gh$p.adj <= .05
names(arm.result) <- str_c(arm.gh$group1, arm.gh$group2, sep="-")
arm.letter <- multcompLetters(arm.result)

# write the anova result to file
armperi.aov <- bind_rows(list(peri=peri.aov, arm=arm.aov), .id="region")
armperi.gh <- bind_rows(list(peri=peri.gh, arm=arm.gh), .id="region")
armperi.result <- list(peri=peri.letter, arm=arm.letter)

write_csv(armperi.aov, file="armperi_anova_result.csv")
write_csv(armperi.gh, file="armperi_games-howell_result.csv")
sink(file="armperi_games-howell_letter.txt")
print(armperi.result)
sink()

#===== Plot CO/F2 summary and distribution in different genomic regions =====
plotCoDistInRegion <- function(dat, regions){
    wt_mean <- group_by(filter(dat, genotype=="wt"), region) %>%
        summarise(meanco=mean(n)) %>%
        filter(region %in% c("arm", "peri"))
    
    p <- ggplot(data=filter(dat, region %in% regions), aes(x=genotype, y=n)) +
        facet_grid(cols=vars(region)) +
        geom_violin(aes(colour=genotype), linewidth=0.2, adjust=1.2) +
        geom_jitter(aes(colour=genotype), width=0.3, height=0.1, size=0.1, stroke=0.1) +
        stat_summary(fun = "mean", geom="crossbar", linewidth=0.1, width=0.6, colour="red") +
        scale_colour_manual(values = pal_GBS) +
        scale_x_discrete(labels=c("wt"="WT", "h2aw"="h2a.w")) +
        theme(axis.title.x = element_blank()) +
        theme(legend.position = "none") +
        theme(text=element_text(size=6, family="Arial", colour="black"),
        axis.text=element_text(size=6, family = "Arial", colour="black")) +
        # theme(axis.text.x = element_text(size=5, vjust=1, hjust=1, angle=20, face="italic")) +
        theme(axis.text.x = element_text(size=6,face="italic")) +
        labs(y="CO per F2") +
        theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))

    return(p)
}

p.codist.armperi <- plotCoDistInRegion(cotable_armperi_count, c("arm", "peri"))
p.codist.arm <- plotCoDistInRegion(cotable_armperi_count, c("arm"))
p.codist.peri <- plotCoDistInRegion(cotable_armperi_count, c("peri"))

cairo_pdf(file="plots/CO_dist_in_armperi.pdf", width=2.5, height=1.5)
print(p.codist.armperi)
dev.off()
png(file="plots/CO_dist_in_armperi.png", width=2.5, height=1.5, unit="in", res=300)
print(p.codist.armperi)
dev.off()

cairo_pdf(file="plots/CO_dist_in_arm.pdf", width=2.5, height=1.5)
print(p.codist.arm)
dev.off()
png(file="plots/CO_dist_in_arm.png", width=2.5, height=1.5, unit="in", res=300)
print(p.codist.arm)
dev.off()

cairo_pdf(file="plots/CO_dist_in_peri.pdf", width=2.5, height=1.5)
print(p.codist.peri)
dev.off()
png(file="plots/CO_dist_in_peri.png", width=2.5, height=1.5, unit="in", res=300)
print(p.codist.peri)
dev.off()
