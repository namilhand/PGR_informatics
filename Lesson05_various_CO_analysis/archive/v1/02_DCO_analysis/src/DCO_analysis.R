library(tidyverse)

#dirin <- "/datasets/data_4/nison/GBS/0_analysis/hcr3-male-female/1_cotables"
#wt_male <- read_csv(file.path(dirin, "wt_male_cotable.csv"), col_names=T)
#input <- wt_male
#wt_female <- read_csv(file.path(dirin, "wt_female_cotable.csv"), col_names=T)
#
#pJ3_mJ3_male <- read_csv(file.path(dirin, "pJ3-mJ3-male_cotable.csv"), col_names=T)
#pJ3_mJ3_female <- read_csv(file.path(dirin, "pJ3-mJ3-female_cotable.csv"), col_names=T)
#setwd("/datasets/data_4/nison/GBS/0_analysis/hcr3-male-female/3_DCO_analysis")


args <- commandArgs(trailingOnly=T)
input <- args[1]
prefix <- args[2]
dirout <- args[3]
n_iter <- as.numeric(args[4])
chr_size=c("Chr1"=30427671, "Chr2"=19698289, "Chr3"=23459830, "Chr4"=18585056, "Chr5"=26975502)

######
#dirin <- "/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/01_dco_cotable/dco_cotable/combine_by_genotype"
#input <- file.path(dirin, "mmH_dco.cotable.txt")
#prefix <- "mmH"
#dirout <- "/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/02_DCO_analysis"
######

input <- read_tsv(input, col_names=T)

####################
# Define functions
####################

# simulateCOS: simulate CO site corresponding to each chromosome of a library
simulateCOS <- function(dat){
		lib <- dat[1] %>%
				as.numeric()
		chr <- dat[2]
		nco <- dat[3]
		cos <- round(runif(nco, min=1, max=chr_size[chr]))
		sim.cotable <- tibble(lib=lib, chr=chr, cos=cos)
		return(sim.cotable)
}

# simulateF2: simulate cotable with same CO distribution
simulateF2 <- function(dat){
		dat <- dat %>%
				arrange(lib, chr)

		cotable.summary <- dat %>%
				group_by(lib, chr) %>%
				summarise(nco=n())
		f2 <- apply(cotable.summary, 1, simulateCOS) %>%
				bind_rows() %>%
				arrange(lib, chr, cos)
		return(f2)
}


# calcDCOdistance: calculate mean DCO distance of cotable
calcDCOdistance <- function(dat){
        dat <- arrange(dat, lib, chr, cos)
        dco_end <- filter(dat, row_number()%%2 == 0)
        dco_start <- filter(dat, row_number()%%2 == 1)
        dco_distance <- dco_end$cos - dco_start$cos
        dco_distance.mean <- mean(dco_distance)
		return(dco_distance.mean)
}

#========================
# calculate DCO distance of observed data and permuted data
#========================
obs <- input
obs.distance <- calcDCOdistance(obs)

sim.distance.repeat <- c()
for(i in 1:n_iter){
    sim <- simulateF2(obs)
    sim.distance <- calcDCOdistance(sim)
    sim.distance.repeat <- c(sim.distance.repeat, sim.distance)
    message(paste0("iteration ", i))
}

mu.sim.distance <- mean(sim.distance.repeat)
sd.sim.distance <- sd(sim.distance.repeat)
#n_iter <- length(sim.distance.repeat)
se <- sd.sim.distance/sqrt(n_iter)
z <- (obs.distance - mu.sim.distance)/se
p <- pnorm(z, lower.tail=FALSE)

obs.tb <- tibble(Group="obs", distance=obs.distance)
sim.tb <- tibble(Group="sim", distance=sim.distance.repeat)
perm_test.tb <- bind_rows(obs.tb, sim.tb)


write_csv(perm_test.tb, file=file.path(dirout, paste0(prefix, "_permutation_test.csv")), col_names=T)


p.perm_test <- ggplot() +
    geom_histogram(data=sim.tb, aes(distance), colour="black") +
    geom_vline(xintercept=mu.sim.distance, colour="blue") +
    geom_vline(xintercept=obs.distance, colour="red") +
    theme(text=element_text(size=10)) +
    theme_classic()

pdf(file=file.path(dirout, paste0(prefix, "_permutation_test.pdf")), width=2.75, height=2.75)
print(p.perm_test)
dev.off()

