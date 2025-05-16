library(qtl2)
library(tidyverse)
library(extrafont)
font_import(pattern="arial", prompt=F)
loadfonts(device="pdf")


#===== 0. read input files =====
physical_map <- "input/dig-col-f2_physical_map.csv"
input <- "dig-col-f2.yaml"
dirout <- "result"
dir.create(dirout, recursive=T)

yaml <- read_cross2(input)
# read input files configured in yaml

#===== 1. calculating genotype probabilites =====
# The first task is to calculate conditional genotype probabilities, given the observed marker data, at each putative QTL position. This is accomplished with the calc_genoprob() function.

# If we wish to perform QTL calculatings at positions between markers (so called "pseudomarkers"), we first need to insert such positions into the genetic map with the function insert_pseudomarkers().

# Here, we don't need pseudomarker as enough marker set is already prepared.

pr <- calc_genoprob(cross=yaml, map=yaml$gmap, cores=32)


#===== 2. performing a genome scan =====
# Calculate LOD score at each marker position using scan1 function
# If no kinship data -> single-QTL model = Haley-Knott regression
# With kinship data, single-QTL model =  ?
out <- scan1(genoprobs = pr, pheno = yaml$pheno, cores=32)

# performing a permutation test
## To perform permutation test to establish the statistical significance of the results of a genome scan, use the function scan1perm().
## In R/qtl, a single function, scanone(), was used for both performing a genome scan and for getting permutation-based significance thresholds, but in R/qtl2, we've decided to make two separate functions.

# permutate what?: permutate genotype so that none of the position become QTL
perm <- scan1perm(genoprobs = pr, pheno = yaml$pheno, cores=32, n_perm=1000)
thr <- summary(perm)
# 5% significance threashold from 1000 permutation = 3.77
# the default behaviour of summary() is to return the 5% significance threshold. You can adjust significance threshold via the "alpha" argument


#===== 3. finding lod peaks =====
# find_peaks() can be used to identify a set of LOD peaks that exceed some threshold.
lod_peak <- find_peaks(out, yaml$gmap, threshold=thr, drop=1.8)
# drop: Defines boundaries as regions where the LOD score drop by a specified value
# ci_lo: lower bound of interval
# ci_hi: upper bound of interval
## chr4: 11,652,333 --- 12,344,849 --- 14,151,468

#===== 4. draw plots =====
pmap.tb <- read_csv(physical_map, col_names=T)
scan.tb <- add_column(pmap.tb, lod=out[,1]) %>%
		mutate(chr=paste0("Chr",chr))
# genetic markers associated with LOD score

## function smoothLod: smooth lod score with loess {{{
smoothLod <- function(x){
		s <- split(x, x$chr)
		lll <- parallel::mclapply(s, function(x) {loess(x$lod ~ x$pos, degree=2, span=0.01)}, mc.cores=5)
		mmm <- lapply(lll, '[[', 'fitted')
		fit <- Reduce(c, mmm)
		fit.tbl <- tibble(x, fit)
		return(fit.tbl)
}
## }}}

scan.fit.tb <- smoothLod(scan.tb)
write_tsv(scan.fit.tb, file="result/dig-col-f2_rQTL_marker-lod.txt")

# A. Genome-wide LOD score
pLodDist <- ggplot(scan.fit.tb, aes(x=pos, y=fit)) +
		geom_line(colour="blue") +
		geom_hline(yintercept=thr[1], colour="red") +
		facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		#ylim(-0.5, 8) +
		scale_x_continuous(labels=~ ./1000000) +
		ylab("LOD score") +
		xlab("Position (Mb)")
pdf(file=file.path(dirout, "dig-col_lod-distribution.pdf"), width=7, height=2.5)
print(pLodDist)
dev.off()

png(file=file.path(dirout, "dig-col_lod-distribution.png"), width=7, height=2.5, unit="in", res=300)
print(pLodDist)
dev.off()


# B. Zoom in by chromosome
pLodDist4 <- ggplot(filter(scan.fit.tb, chr=="Chr4"), aes(x=pos, y=fit)) +
		geom_line(colour="blue") +
		geom_hline(yintercept=thr[1], colour="red") +
		geom_segment(aes(x=pos, xend=pos, y=0, yend=0.15), colour="black", size=0.05, alpha=0.05) +
		facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		scale_x_continuous(labels=~ ./1000000) +
		xlab("Position (Mb)") +
		ylab("LOD score")
pdf(file=file.path(dirout, "dig-col_lod-distribution_chr4.pdf"), width=3, height=2.5)
print(pLodDist4)
dev.off()

png(file=file.path(dirout, "dig-col_lod-distribution_chr4.png"), width=3, height=2.5, unit="in", res=300)
print(pLodDist4)
dev.off()

# C. Zoom in by peaks
## peaks: LOD peaks of rQTL result.
## chr:ci_lo---pos---ci_hi
## chr4: 11,652,333 --- 12,344,849 --- 14,151,468
peak1 <- filter(scan.fit.tb, chr == "Chr4" & pos > 11652333 & pos < 14151468)

pDist.peak1 <- ggplot(peak1) +
		#geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=taf4b[1]/1000000, xmax=taf4b[2]/1000000), colour="yellow", alpha=0.4) +
		geom_point(colour="black", alpha=0.5, aes(x=pos, y=lod), size=0.3) +
		geom_line(colour="blue", aes(x=pos, y=fit)) +
		geom_hline(yintercept=thr[1], colour="red") +
		geom_segment(aes(x=pos, xend=pos, y=0, yend=0.15), colour="black", size=0.05, alpha=0.05) +
		scale_x_continuous(labels=~ ./1000000) +
		#facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		#ylim(-0.5, 8) +
		#xlim(7.5, 15) +
		xlab("Position (Mb)") +
		ylab("LOD score")
pdf(file="result/dig-col_lod-distribution_peak1.pdf", width=2, height=2)
print(pDist.peak1)
dev.off()

png(file="result/dig-col_lod-distribution_peak1.png", width=2, height=2, unit="in", res=300)
print(pDist.peak1)
dev.off()
