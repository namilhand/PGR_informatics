# This script is for detecting recessive causal allele
library(tidyverse)

args <- commandArgs(trailingOnly=T)
inputFile <- args[1]
dirout <- args[2]
#inputFile <- "/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/04_result/hcr7_allvar.bulkedSeg.allvar.tsv"

prefix <- basename(inputFile) %>% str_replace(".se.tsv", "")
input <- read_delim(inputFile, col_names=T, delim="\t")

#tableColName <- c("chr", "pos", "ref", "alt", "type", "wt", "mut", "dp.wt", "dp.mut")
#tableColName <- c("chr", "pos", "ref", "alt", "wt.ad", "mut.ad", "wt.dp", "mut.dp", "ratio")
#colnames(input) <- tableColName

input.ctga <- input %>%
    filter(wt.alt != 0, mut.alt !=0) %>%
    filter((ref == "C" & alt == "T") | (ref == "G" & alt == "A")) %>%
    filter(chr != "ChrM", chr != "ChrC")

write_tsv(input.ctga, file.path(dirout, paste0(prefix, ".refined_markers.tsv")))


pEMS <- ggplot(input.ctga, aes(x=pos, y=ratio)) +
		geom_point(size=0.6) +
		facet_grid(rows=vars(chr), scales = "free_x") +
		geom_smooth(colour="blue", method="loess", span = 0.1) +
		scale_y_continuous(limits=c(-1, 1)) +
		scale_x_continuous(labels= ~ round(.x/1000000, digits=0)) +
		theme_bw() +
		labs(x="Position (Mb)", y="ratio")


pdf(file=file.path(dirout, paste0(prefix, ".AF_map.pdf")))
print(pEMS)
dev.off()
