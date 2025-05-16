# R/qtl2

- R/qtl2 package is used for this analysis. It can be installed from CRAN:
    - `install.packages("qtl2")`
- R/qtl2 (aka qtl2) is an R package for QTL analysis.
    - https://kbroman.org/qtl2/
    - https://link.springer.com/book/10.1007/978-0-387-92125-9 (We have this book in the bookshelf in room #305!)

- Here, I used phenotype and genotype data of Col x Di-G F2 (Provided by Seula Lee). The mean number of seed in a selique was measured for each F2 plant.

- For QTL analysis we need:
    1. genotype data of each F2 plant at each genetic markers
    2. Phenotype of the F2 plant (here, mean seed number)
    3. Genetic map of markers (here, SNPs were used as marker)

Below is examples of input files prepared for qtl2 package.

- Use conda environment “qtl”. The recipe yaml file is included in the [**Lesson02**](https://github.com/namilhand/PGR_informatics/tree/main/Lesson02_FilteringSNPs)

# Inputs required for qtl2

1. Genotype data

| F2_id | 1_1 | 1_2 | … | 5_1530 | … |
| --- | --- | --- | --- | --- | --- |
| f2_1 | LL | LL | LL | CC |  |
| f2_2 | LL | LL | CC | CC |  |
| f2_3 | CC | CC | CC | CC |  |
| … |  |  |  |  |  |
- From column2 to the last are SNP markers
    - 1_1: Chr1, SNP1
    - 5_1530: Chr5, SNP1530
    - …

 

1. Phenotype data

| F2_id | seed |
| --- | --- |
| f2_1 | 57.83 |
| f2_2 | 66.33 |
| f2_3 | 57.4 |
| … | … |

1. Physical map of markers (SNPs)

| marker | chr | pos |
| --- | --- | --- |
| 1_1 | 1 | 711 |
| 1_2 | 1 | 892 |
| 1_3 | 1 | 10904 |
| … |  |  |

Then, how to prepare the inputs? Let’s move on to the next section.

# Preparing inputs

```r
library(qtl2)
library(tidyverse)

################
# For running qtl2, we need:
# 1. genotype data (from F2 sequencing)
# 2. genetic map (converting physical map to Mbp)
# - Since the exact physical position of each SNP marker is known,
#   we can skip using a geneticm map. Simply reformat the position from bp to Mbp
# 3. phenotypes
################

#===== 1. Genotype data =====

# Use genotype data of Di-G x Col F2 (239 F2 plants)
# The genotype data can be found in the server 141.223.132.213 at the location below
set123 <- "/datasets/data_4/nison/GBS/GBS_marker_v2/20240329_Di-G_set123/results/06_tiger"

cohort.column_order <- paste0("f2_", 1:239)

# readHMMbulk: Read HMM-output in bulk (i.e. *.rough.co from TIGER output).
## HMM inferred genotype from each file is saved as an element in a list
## The list is then combined into a single dataframe.

## Why HMM-inferred genotype? Genotype fidelity is low in heterozygous regions,
## and HMM significantly improves accuracy.

## Be aware that TIGER has an error that HMM-corrected genotype is switched in that "LL" to "CL" and "CL" to "LL", except for the last chromosome
## readHMMbulk function correct the mislabeled genotype in Chr1-4
readHMMbulk <- function(dirin){
		files <- list.files(dirin)
		files.hmm <- files[grep("rough.co$", files)]
		lib.nums <- length(files.hmm)

		cohort.gt <- NULL 

		for(k in 1:lib.nums){
				hmm <- file.path(dirin, files.hmm[k])
				dat <- read_tsv(file=hmm, col_names=F)
				colnames(dat) <- c("lib", "chr", "pos", "basecall", "hmm_1", "hmm_2", "ref", "ref.count", "alt", "alt.count")

				# Be aware that TIGER has an error that HMM-corrected genotype is switched in that "LL" to "CL" and "CL" to "LL", except for the last chromosome
				dat_1 <- filter(dat, chr %in% 1:4) %>%
					mutate(hmm_1 = case_when(hmm_1 == "LL" ~ "CL",
											 hmm_1 == "CL" ~ "LL",
											 hmm_1 == "CC" ~ "CC",
											 TRUE ~ hmm_1)) %>%
					mutate(hmm_2 = case_when(hmm_2 == "LL" ~ "CL",
											 hmm_2 == "CL" ~ "LL",
											 hmm_2 == "CC" ~ "CC",
											 TRUE ~ hmm_2))
				dat_2 <- filter(dat, chr == 5)

				dat_corrected <- bind_rows(dat_1, dat_2)

				gt <- dat_corrected$hmm_1
				lib_name <- str_replace(dat_corrected$lib[1], "results/06_tiger/lib", "")
                lib_name <- str_replace(lib_name, "_MappedOn_tair10", "")
				lib_name <- paste0("f2_", lib_name)

				cohort.gt[[k]] <- gt
				names(cohort.gt)[k] <- lib_name
		}
		cohort.gt.df <- as.data.frame(cohort.gt)
		# arrrange by library number
		cohort.gt.df <- cohort.gt.df[,cohort.column_order]
		cohort.gt.df <- bind_cols(dat[,2:3], cohort.gt.df)
		return(cohort.gt.df)
		}
gt.set123 <- readHMMbulk(set123)
cohort.gt <- gt.set123

gt <- cohort.gt[,3:241]
gt <- rbind(colnames(gt), gt)
gt_col1 <- c("id", phys_map_id)
gt_id <- cbind(gt_col1, gt)
gt_t <- t(gt_id)
gt_df <- as_tibble(data.frame(gt_t))

write_csv(gt_df, file="input/dig-col-f2_genotypes.csv", col_names=F)

#===== 2. physical map of markers =====
phys_map <- cohort.gt[,1:2]

phys_map_chr_table <- table(phys_map$chr)
phys_map_id <- c(paste0("1_", 1:phys_map_chr_table[1]),
				 paste0("2_", 1:phys_map_chr_table[2]),
				 paste0("3_", 1:phys_map_chr_table[3]),
				 paste0("4_", 1:phys_map_chr_table[4]),
				 paste0("5_", 1:phys_map_chr_table[5]))
phys_map <- add_column(phys_map, marker=phys_map_id, .before="chr")

write_csv(phys_map, file="input/dig-col-f2_physical_map.csv", col_names=T)

#===== 3. phenotype =====
phenotype <- read_csv("data/dig-col_f2_ovule_phenotype.CSV", col_names=T) %>%
    dplyr::select(c(1, 3)) %>%
	mutate(sample_no = paste0("f2_", sample_no))
colnames(phenotype) <- c("id", "seed")

write_csv(phenotype, file="input/dig-col-f2_phenotype.csv", col_names=T)
```

Now that all inputs are prepared, we need to write an configuration file (.yaml) for qtl2 package as below.

```yaml
crosstype: f2
geno: dig-col-f2_genotypes.csv
pheno: dig-col-f2_phenotype.csv
gmap: dig-col-f2_physical_map.csv
alleles:
- C
- L
genotypes:
  CC: 1
  CL: 2
  LL: 3
```

### Note

- Always be aware that TIGER pipeline contains an error where HMM-corrected genotype are switched: "LL" is mislabeled as "CL" and "CL" as "LL", except for the last chromosome.
The mislabeled genotyep of Chr1 to Chr4 are corrected in the above script.
- **How to get reliable SNP list? → See [Lesson02](https://github.com/namilhand/PGR_informatics/tree/main/Lesson02_FilteringSNPs)**

Now, all is set. Let’s run qtl2 package!

## Running qtl2

### 1 . Setting packages and reading input files

```r
library(qtl2)
library(tidyverse)
library(org.At.tair.db)
library(TxDb.Athaliana.BioMart.plantsmart22)
library(extrafont)
font_import(pattern="arial", prompt=F)
loadfonts(device="pdf")

# For annotation
txdb <- TxDb.Athaliana.BioMart.plantsmart22
seqlevels(txdb) <- "5"
genes5 <- as_tibble(genes(txdb)) %>%
		mutate(seqnames=paste0("Chr", seqnames))
colnames(genes5)[1] <- "chr"

#===== 0. read input files =====
physical_map <- "input/dig-col-f2_physical_map.csv"
input <- "dig-col-f2.yaml"
dirout <- "result"
gtfpath <- "data/Araport11_GTF_genes_transposons.current.gtf"
dir.create(dirout, recursive=T)

yaml <- read_cross2(input)
# read input files configured in yaml

gtf <- read_tsv(gtfpath, col_names=F) %>%
		dplyr::select(c(X1, X3, X4, X5, X7, X9))
colnames(gtf) <- c("chr", "type", "start", "end", "strand", "id")
```

- `read_cross2` reads input files configured in `dig-col-f2.yaml`
- `gtf` is for the annotating genes after QTL mapping

- To see the input summary, type `yaml` and you will see:

| Category | Count |
| --- | --- |
| Total individuals | 239 |
| Genotyped individuals | 239 |
| Phenotyped individuals | 239 |
| With both geno & pheno | 239 |
| Number of phenotypes | 1 |
| Number of covariates | 0 |
| Phenotype covariates | 0 |

| Chromosome | Number of Markers |
| --- | --- |
| 1 | 92,789 |
| 2 | 63,941 |
| 3 | 69,695 |
| 4 | 59,445 |
| 5 | 87,378 |
| Total | 373,248 |

### Calculating genotype probabilities

The first task is to calculate conditional genotype probabilities, given the observed marker data, at each putative QTL position. This is accomplished with the `calc_genoprob` function.

If we wish to perform QTL calculatings at positions between markers (so called "pseudomarkers"), we first need to insert such positions into the genetic map with the function `insert_pseudomarkers`.

Here, we don't need pseudomarker as enough marker set is already prepared.

```r
#===== 1. calculating genotype probabilites =====
# The first task is to calculate conditional genotype probabilities, given the observed marker data, at each putative QTL position. This is accomplished with the calc_genoprob() function.

# If we wish to perform QTL calculatings at positions between markers (so called "pseudomarkers"), we first need to insert such positions into the genetic map with the function insert_pseudomarkers().

# Here, we don't need pseudomarker as enough marker set is already prepared.

pr <- calc_genoprob(cross=yaml, map=yaml$gmap, cores=32)
```

Genotype probability of f2 plants at each marker is calculated and stored in `pr`. It doesn’t have any practical meaning in our case actually, since we’ve already corrected and confirmed genotypes using Hidden-Markov model (`readHMMbulk`). Here marker 1 of Chr 1 is shown at below.

|  | CC | CL | LL |
| --- | --- | --- | --- |
| f2_1 | 2.496541e-05 | 9.999501e-01 | 2.496541e-05 |
| f2_2 | 2.496541e-05 | 9.999501e-01 | 2.496541e-05 |
| f2_3 | 2.496541e-05 | 9.999501e-01 | 2.496541e-05 |
| f2_4 | 9.998603e-01 | 9.478041e-05 | 4.491844e-05 |
| f2_5 | 4.491844e-05 | 9.478041e-05 | 9.998603e-01 |
| f2_6 | 9.998603e-01 | 9.478041e-05 | 4.491844e-05 |
| f2_7 | 4.491844e-05 | 9.478041e-05 | 9.998603e-01 |
| f2_8 | 9.998603e-01 | 9.478041e-05 | 4.491844e-05 |
| f2_9 | 2.496541e-05 | 9.999501e-01 | 2.496541e-05 |
| f2_10 | 9.998603e-01 | 9.478041e-05 | 4.491844e-05 |
| f2_11 | 2.496541e-05 | 9.999501e-01 | 2.496541e-05 |

### Calculaing LOD score at each marker position (genome scan)

```r
#===== 2. performing a genome scan =====
# Calculate LOD score at each marker position using scan1 function
# If no kinship data -> single-QTL model = Haley-Knott regression
# With kinship data, single-QTL model =  ?
out <- scan1(genoprobs = pr, pheno = yaml$pheno, cores=32)

```

output: LOD score at each marker position

| marker | LOD score |
| --- | --- |
| 1_1 | 1.680702 |
| 1_2 | 1.680702 |
| 1_3 | 1.680702 |
| 1_4 | 1.680702 |
| 1_5 | 1.680702 |
| 1_6 | 1.680702 |

What is LOD score?

- Log10(likelihood that a locus is linked with phenotype/likelihood of a locus is not lnked with phenotype)
- If LOD=3 → likelihood that a locus is linked with phenotype is 1000-fold higher than the opposite.
- For detail explanation? See “A guide to QTL mapping with R/qtl” by Karl Broman, pg.76-77

### Finding LOD score threshold

Since LOD scores are calculated at millions of markers, some loci may may have high LOD score by chance.

We perform a permutation test (`scan1perm`) to identify the LOD score threshold exceeded by only 5% of loci due to random chance. You can change the threshold by adjusting `alpha` argument in `scan1perm` function.

```r
# performing a permutation test
## To perform permutation test to establish the statistical significance of the results of a genome scan, use the function scan1perm().
## In R/qtl, a single function, scanone(), was used for both performing a genome scan and for getting permutation-based significance thresholds, but in R/qtl2, we've decided to make two separate functions.

# permutate what?: permutate genotype so that none of the position become QTL
perm <- scan1perm(genoprobs = pr, pheno = yaml$pheno, cores=32, n_perm=1000)
thr <- summary(perm)
# 5% significance threashold from 1000 permutation = 3.77
# the default behaviour of summary() is to return the 5% significance threshold. You can adjust significance threshold via the "alpha" argument
```

LOD threshold = 3.77

### Identifying LOD peaks

`find_peaks` can be used to identify a set of LOD peaks that exceed threshold.

`drop`  defines boundaries as regions where the LOD score drop by a specified value

```r
#===== 3. finding lod peaks =====
# find_peaks() can be used to identify a set of LOD peaks that exceed some threshold.
lod_peak <- find_peaks(out, yaml$gmap, threshold=thr, drop=1.8)
# drop: Defines boundaries as regions where the LOD score drop by a specified value
# ci_lo: lower bound of interval
# ci_hi: upper bound of interval
## chr4: 11,652,333 --- 12,344,849 --- 14,151,468

```

`find_peaks` found 4 LOD peaks as below. Among the peaks, peak at Chr4:12344849 is the most promising one.

output:

| lodindex | lodcolumn | chr | pos | lod | ci_lo | ci_hi |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | seed | 1 | 21897271 | 5.112728 | 3774050 | 22892932 |
| 1 | seed | 3 | 9435545 | 6.894646 | 9434566 | 9452964 |
| 1 | seed | 4 | 12344849 | 27.393547 | 11652333 | 14151468 |
| 1 | seed | 5 | 17603686 | 7.213653 | 17596399 | 17606520 |

### Drawing plots

Now let’s draw the QTL mapping result.

*Smoothing LOD score*

```r
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
```

example of `scan.fit.tb`

| marker | chr | pos | lod | fit |
| --- | --- | --- | --- | --- |
| 1_1 | Chr1 | 711 | 1.68 | 1.89 |
| 1_2 | Chr1 | 892 | 1.68 | 1.89 |
| 1_3 | Chr1 | 10904 | 1.68 | 1.89 |
| 1_4 | Chr1 | 32210 | 1.68 | 1.90 |
| 1_5 | Chr1 | 37388 | 1.68 | 1.90 |
| 1_6 | Chr1 | 38554 | 1.68 | 1.90 |
| 1_7 | Chr1 | 71326 | 1.68 | 1.90 |
| 1_8 | Chr1 | 71348 | 1.68 | 1.90 |
| 1_9 | Chr1 | 88300 | 1.68 | 1.91 |
| 1_10 | Chr1 | 90571 | 1.68 | 1.91 |

*codes for drawing LOD landscape*

```r
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

# B. Zoom in by chromosome
pLodDist1 <- ggplot(filter(scan.fit.tb, chr=="Chr1"), aes(x=pos, y=fit)) +
		geom_line(colour="blue") +
		geom_hline(yintercept=thr[1], colour="red") +
		geom_segment(aes(x=pos, xend=pos, y=0, yend=0.15), colour="black", size=0.05, alpha=0.05) +
		facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		ylim(-0.5, 8) +
		xlab("Position (Mb)") +
		ylab("LOD score")
pdf(file=file.path(dirout, "dig-col_lod-distribution_chr1.pdf"), width=3, height=2.5)
print(pLodDist1)
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
		#facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		#ylim(-0.5, 8) +
		#xlim(7.5, 15) +
		xlab("Position (Mb)") +
		ylab("LOD score")
pdf(file="result/dig-col_lod-distribution_peak1.pdf", width=2, height=2)
print(pDist.peak1)
dev.off()
```

1. Genome-wide LOD landscape
![Fig.1](https://github.com/namilhand/PGR_informatics/blob/main/Lesson03_QTLmapping/img/dig-col_lod-distribution.png)

2. LOD landscape at Chr. 4
![Fig.2](https://github.com/namilhand/PGR_informatics/blob/main/Lesson03_QTLmapping/img/dig-col_lod-distribution_chr4.png)

3. LOD landscape at the peak
![Fig.3](https://github.com/namilhand/PGR_informatics/blob/main/Lesson03_QTLmapping/img/dig-col_lod-distribution_peak1.png)



# Full code

```r
library(qtl2)
#library(qtl2ggplot)
library(tidyverse)
library(org.At.tair.db)
library(TxDb.Athaliana.BioMart.plantsmart22)
library(extrafont)
font_import(pattern="Arial", prompt=F)
loadfonts(device="pdf")

# For annotation
txdb <- TxDb.Athaliana.BioMart.plantsmart22
seqlevels(txdb) <- "5"
genes5 <- as_tibble(genes(txdb)) %>%
		mutate(seqnames=paste0("Chr", seqnames))
colnames(genes5)[1] <- "chr"

#===== 0. read input files =====
physical_map <- "input/dig-col-f2_physical_map.csv"
input <- "dig-col-f2.yaml"
dirout <- "result"
gtfpath <- "/home/nison/work/refgenome/araport11/Araport11_GTF_genes_transposons.current.gtf"
dir.create(dirout, recursive=T)

yaml <- read_cross2(input)
# read input files configured in yaml

gtf <- read_tsv(gtfpath, col_names=F) %>%
		dplyr::select(c(X1, X3, X4, X5, X7, X9))
colnames(gtf) <- c("chr", "type", "start", "end", "strand", "id")

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

# B. Zoom in by chromosome
pLodDist1 <- ggplot(filter(scan.fit.tb, chr=="Chr1"), aes(x=pos, y=fit)) +
		geom_line(colour="blue") +
		geom_hline(yintercept=thr[1], colour="red") +
		geom_segment(aes(x=pos, xend=pos, y=0, yend=0.15), colour="black", size=0.05, alpha=0.05) +
		facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		ylim(-0.5, 8) +
		xlab("Position (Mb)") +
		ylab("LOD score")
pdf(file=file.path(dirout, "dig-col_lod-distribution_chr1.pdf"), width=3, height=2.5)
print(pLodDist1)
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
		#facet_wrap(facet="chr", nrow=1, ncol=5, scales="free_x") +
		#ylim(-0.5, 8) +
		#xlim(7.5, 15) +
		xlab("Position (Mb)") +
		ylab("LOD score")
pdf(file="result/dig-col_lod-distribution_peak1.pdf", width=2, height=2)
print(pDist.peak1)
dev.off()
```

