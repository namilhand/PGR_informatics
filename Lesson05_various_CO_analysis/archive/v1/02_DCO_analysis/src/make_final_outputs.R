library(tidyverse)

dir_dco <- "/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/02_DCO_analysis"

col_dco <- read_csv(file.path(dir_dco, "col_permutation_test.csv"))
hta6_dco <- read_csv(file.path(dir_dco, "hta6_permutation_test.csv"))
hta7_dco <- read_csv(file.path(dir_dco, "hta7_permutation_test.csv"))
hta12_dco <- read_csv(file.path(dir_dco, "hta12_permutation_test.csv"))
hta67_dco <- read_csv(file.path(dir_dco, "hta67_permutation_test.csv"))
h2aw_dco <- read_csv(file.path(dir_dco, "h2aw_permutation_test.csv"))
mmH_dco <- read_csv(file.path(dir_dco, "mmH_permutation_test.csv"))
mmS_dco <- read_csv(file.path(dir_dco, "mmS_permutation_test.csv"))
mmHS_dco <- read_csv(file.path(dir_dco, "mmHS_permutation_test.csv"))


calcPval <- function(dat){
    #dat <- col_dco

    obs.distance <- mean(filter(dat, Group == "obs")$distance)
    mu.sim.distance <- mean(filter(dat, Group == "sim")$distance)
    sd.sim.distance <- sd(filter(dat, Group == "sim")$distance)

    #n_iter <- length(sim.distance.repeat)
    se <- sd.sim.distance/sqrt(nrow(filter(dat, Group == "sim")))
    z <- (obs.distance - mu.sim.distance)/se
    p <- pnorm(z, lower.tail=FALSE)
    return(p)
}

pval.col <- calcPval(col_dco)
pval.hta6 <- calcPval(hta6_dco)
pval.hta7 <- calcPval(hta7_dco)
pval.hta12 <- calcPval(hta12_dco)
pval.hta67 <- calcPval(hta67_dco)
pval.h2aw <- calcPval(h2aw_dco)
pval.mmH <- calcPval(mmH_dco)
pval.mmS <- calcPval(mmS_dco)
pval.mmHS <- calcPval(mmHS_dco)

drawPlot <- function(dat){
#    dat <- col_dco

    sim.tb <- filter(dat, Group == "sim")
    obs.tb <- filter(dat, Group == "obs")
    mu.sim.distance <- mean(sim.tb$distance)
    obs.distance <- mean(obs.tb$distance)

    p <- ggplot() +
         geom_histogram(data=sim.tb, aes(distance), colour="black") +
         geom_vline(xintercept=mu.sim.distance, colour="blue") +
         geom_vline(xintercept=obs.distance, colour="red") +
         theme_classic() +
         theme(text=element_text(size=8, family="Arial", colour="black")) +
         theme(axis.text = element_text(size=6, family="Arial", colour="black")) +
         theme(axis.title = element_text(size=6, family="Arial", colour="black")) +
         theme(axis.line = element_line(linewidth=0.2), axis.ticks = element_line(linewidth=0.2))
     return(p)
}

pdf(file=file.path(dirout, paste0(prefix, "_permutation_test.pdf")), width=2.75, height=2.75)
print(p.perm_test)
dev.off()

