library(tidyverse)
library(extrafont)
#font_import(pattern="Arial", prompt=F)
#loadfonts(device="pdf")

# Set dir_dco accordingly
dir_dco <- "/datasets/data_4/nison/PGR_informatics/Lesson05_various_CO_analysis/1_DCO_analysis/results/3_dco_analysis"
test_dco <- read_csv(file.path(dir_dco, "test_permutation_test.csv"))

dir_result <- file.path(dir_dco, "final_plot")
dir.create(dir_result, recursive=T)


all_dco_distance <- c(test_dco$distance)

x_limits <- c(min(all_dco_distance) * 0.9, max(all_dco_distance) * 1.1)


calcPval <- function(dat){
    pval <- sum(filter(dat, Group == "sim")$distance > filter(dat, Group == "obs")$distance)/2000
    return(pval)
}

pval <- calcPval(test_dco)
print(paste0("pvalue = ",pval))


drawPlot <- function(dat){
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

plot.dco <- drawPlot(test_dco)

pdf(file=file.path(dir_result, "test_permutation_test.pdf"), width=2, height=1.5)
print(plot.dco)
dev.off()
