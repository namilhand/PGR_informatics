library(qtl2)
library(tidyverse)

################
# For running qtl2, we need:
# 1. genotypes (from F2 sequencing)
# 2. genetic map (converting physical map to Mbp)
# - Since the exact physical position of each SNP marker is known,
#   we can skip using a geneticm map. Simply reformat the position from bp to Mbp
# 3. phenotypes
# 4. control file
################


#===== 1. Genotype data =====

# Use genotype data of Di-G x Col F2 (239 F2 plants)
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

