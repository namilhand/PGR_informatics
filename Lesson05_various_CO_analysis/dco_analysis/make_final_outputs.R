library(tidyverse)
library(extrafont)
font_import(pattern="Arial", prompt=F)
loadfonts(device="pdf")

dir_dco <- "/home/namilhand/01_Projects/H2AW/h2aw_server/dco_analysis/v1/02_DCO_analysis"
dir_dco_v2 <- "/home/namilhand/01_Projects/H2AW/h2aw_server/dco_analysis/v2/02_DCO_analysis"
dir_results <- file.path(dir_dco_v2, "results")
dir.create(dir_results, recursive=T)

col_dco <- read_csv(file.path(dir_dco, "col_permutation_test.csv"))
hta6_dco <- read_csv(file.path(dir_dco, "hta6_permutation_test.csv"))
hta7_dco <- read_csv(file.path(dir_dco, "hta7_permutation_test.csv"))
hta12_dco <- read_csv(file.path(dir_dco, "hta12_permutation_test.csv"))
hta67_dco <- read_csv(file.path(dir_dco, "hta67_permutation_test.csv"))
h2aw_dco <- read_csv(file.path(dir_dco, "h2aw_permutation_test.csv"))
mmH_dco <- read_csv(file.path(dir_dco, "mmH_permutation_test.csv"))
mmS_dco <- read_csv(file.path(dir_dco, "mmS_permutation_test.csv"))
mmHS_dco <- read_csv(file.path(dir_dco, "mmHS_permutation_test.csv"))
recq4_dco <- read_csv(file.path(dir_dco_v2, "recq4_permutation_test.csv"))

all_dco_distance <- c(col_dco$distance, hta6_dco$distance, hta7_dco$distance, hta12_dco$distance, hta67_dco$distance, h2aw_dco$distance, mmH_dco$distance, mmS_dco$distance, mmHS_dco$distance)

x_limits <- c(min(all_dco_distance) * 0.9, max(all_dco_distance) * 1.1)
x_limits.recq4 <- c(min(recq4_dco$distance) * 0.9, max(recq4_dco$distance) * 1.1)


calcPval <- function(dat){
    # dat <- hta67_dco
    pval <- sum(filter(dat, Group == "sim")$distance > filter(dat, Group == "obs")$distance)/2000
    return(pval)

    # obs.distance <- mean(filter(dat, Group == "obs")$distance)
    # mu.sim.distance <- mean(filter(dat, Group == "sim")$distance)
    # sd.sim.distance <- sd(filter(dat, Group == "sim")$distance)

    # #n_iter <- length(sim.distance.repeat)
    # se <- sd.sim.distance/sqrt(192)
    # z <- (obs.distance - mu.sim.distance)/sd.sim.distance
    # p <- pnorm(z, lower.tail=FALSE)
    # return(p)
}

pval.col <- calcPval(col_dco)
# 0
pval.hta6 <- calcPval(hta6_dco)
# 0
pval.hta7 <- calcPval(hta7_dco)
# 0.0055
pval.hta12 <- calcPval(hta12_dco)
# 0
pval.hta67 <- calcPval(hta67_dco)
# 0.002
pval.h2aw <- calcPval(h2aw_dco)
# 0
pval.mmH <- calcPval(mmH_dco)
# 0
pval.mmS <- calcPval(mmS_dco)
# 0
pval.mmHS <- calcPval(mmHS_dco)
# 0
pval.recq4 <- calcPval(recq4_dco)
# 0.779

drawPlot <- function(dat){
#    dat <- col_dco

    sim.tb <- filter(dat, Group == "sim")
    obs.tb <- filter(dat, Group == "obs")
    mu.sim.distance <- mean(sim.tb$distance)
    obs.distance <- mean(obs.tb$distance)

    p <- ggplot() +
         geom_histogram(data=sim.tb, aes(distance), bins=60, colour="grey60", fill="grey60", linewidth=0.1) +
         geom_vline(xintercept=mu.sim.distance, colour="blue", linewidth=0.2) +
         geom_vline(xintercept=obs.distance, colour="red", linewidth=0.2) +
         scale_x_continuous(limits=x_limits, labels=scales::label_number(scale = 1/1000000), breaks=c(seq(0, 10000000, 2*10^6), 12*10^6)) +
        # scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mnase_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
         labs(y="Frequency", x="Double crossover distance (Mb)") +
         theme_classic() +
         theme(text=element_text(size=8, family="Arial", colour="black")) +
         theme(axis.text = element_text(size=7, family="Arial", colour="black")) +
         theme(axis.title = element_text(size=7, family="Arial", colour="black")) +
         theme(axis.line = element_line(linewidth=0.2), axis.ticks = element_line(linewidth=0.2))
     return(p)
}

plot.col <- drawPlot(col_dco)
plot.hta6 <- drawPlot(hta6_dco)
plot.hta7 <- drawPlot(hta7_dco)
plot.hta12 <- drawPlot(hta12_dco)
plot.hta67 <- drawPlot(hta67_dco)
plot.h2aw <- drawPlot(h2aw_dco)
plot.mmH <- drawPlot(mmH_dco)
plot.mmS <- drawPlot(mmS_dco)
plot.mmHS <- drawPlot(mmHS_dco)
plot.recq4 <- drawPlot(recq4_dco) + 
         scale_x_continuous(limits=x_limits.recq4, labels=scales::label_number(scale = 1/1000000), breaks=c(seq(0, 10000000, 0.5*10^6), 12*10^6))

pdf(file=file.path(dir_results, "col_permutation_test.pdf"), width=2, height=1.5)
print(plot.col)
dev.off()
pdf(file=file.path(dir_results, "hta6_permutation_test.pdf"), width=2, height=1.5)
print(plot.hta6)
dev.off()
pdf(file=file.path(dir_results, "hta7_permutation_test.pdf"), width=2, height=1.5)
print(plot.hta7)
dev.off()
pdf(file=file.path(dir_results, "hta12_permutation_test.pdf"), width=2, height=1.5)
print(plot.hta12)
dev.off()
pdf(file=file.path(dir_results, "hta67_permutation_test.pdf"), width=2, height=1.5)
print(plot.hta67)
dev.off()
pdf(file=file.path(dir_results, "h2aw_permutation_test.pdf"), width=2, height=1.5)
print(plot.h2aw)
dev.off()
pdf(file=file.path(dir_results, "mmH_permutation_test.pdf"), width=2, height=1.5)
print(plot.mmH)
dev.off()
pdf(file=file.path(dir_results, "mmS_permutation_test.pdf"), width=2, height=1.5)
print(plot.mmS)
dev.off()
pdf(file=file.path(dir_results, "mmHS_permutation_test.pdf"), width=2, height=1.5)
print(plot.mmHS)
dev.off()
pdf(file=file.path(dir_results, "recq4_permutation_test.pdf"), width=2, height=1.5)
print(plot.recq4)
dev.off()

#=======================
# Fig. 3F
#=======================

fig3f <- function(dat){
#    dat <- col_dco

    sim.tb <- filter(dat, Group == "sim")
    obs.tb <- filter(dat, Group == "obs")
    mu.sim.distance <- mean(sim.tb$distance)
    obs.distance <- mean(obs.tb$distance)

    p <- ggplot() +
         geom_histogram(data=sim.tb, aes(distance), bins=60, colour="grey60", fill="grey60", linewidth=0.1) +
         geom_vline(xintercept=mu.sim.distance, colour="blue", linewidth=0.2) +
         geom_vline(xintercept=obs.distance, colour="red", linewidth=0.2) +
         scale_x_continuous(limits=x_limits, labels=scales::label_number(scale = 1/1000000), breaks=c(seq(0, 10000000, 2*10^6), 12*10^6)) +
        # scale_x_continuous(name = "Coordinates (Mb)", labels=scales::label_number(scale = 1/1000000), breaks=c(seq(1, max(mnase_ma$cumwindow), 20*10^6), 120000000), expand=c(0.01, 0.01)) +
         labs(y="Frequency", x="Distance (Mb)") +
         theme_classic() +
         theme(text=element_text(size=6, family="Arial", colour="black")) +
         theme(axis.text = element_text(size=6, family="Arial", colour="black")) +
         theme(axis.title = element_text(size=6, family="Arial", colour="black")) +
         theme(axis.line = element_line(linewidth=0.2), axis.ticks = element_line(linewidth=0.2))
     return(p)
}

fig3f.col <- fig3f(col_dco)
fig3f.h2aw <- fig3f(h2aw_dco)

pdf(file=file.path(dir_results, "Fig3F_col.pdf"), width=1.8, height=1.2)
print(fig3f.col)
dev.off()
pdf(file=file.path(dir_results, "Fig3F_h2aw.pdf"), width=1.8, height=1.2)
print(fig3f.h2aw)
dev.off()