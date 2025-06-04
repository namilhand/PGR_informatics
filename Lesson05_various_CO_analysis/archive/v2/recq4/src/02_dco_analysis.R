library(tidyverse)

# 1. calculate dco distance of provided dco_table
# 2. Make permuted data of the cotable which corresponds to the dco_table
# 3. Calculate mean co-co distance of the permuted cotable
# 4. Iterate step 2-3 for permutation test
dirout <- "/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/v2/02_dco_analysis"
prefix <- "recq4"

# cotable
recq4.cotable.file <- "/datasets/data_4/nison/GBS/GBS_marker_v2/20210322_recq4/results/recq4_2021_cotable.txt"

recq4.cotable <- read_csv(recq4.cotable.file) %>%
        mutate(lib = str_replace(lib, "results/06_tiger/lib", "")) %>%
        mutate(lib = str_replace(lib, "_MappedOn_tair10", ""))

colnames(recq4.cotable)[2] <- "chr"

recq4.cotable.summary <- group_by(recq4.cotable, lib, chr) %>%
    summarize(ncos=n())

# dco_table
recq4.dcotable <- read_tsv("/datasets/data_4/nison/h2aw/git_h2aw_analysis/dco_analysis/v2/01_dco_table/recq4_dco.cotable.txt")

chr_size=c("Chr1"=30427671, "Chr2"=19698289, "Chr3"=23459830, "Chr4"=18585056, "Chr5"=26975502)


####################
# Define functions
####################

# simulateCOCOdist: simulate mean CO-CO distance of each chromosome of a library
simulateCOCOdist <- function(dat){
		lib <- dat[1] %>%
				as.numeric()
		chr <- dat[2]
		nco <- dat[3]
		cos <- sort(round(runif(nco, min=1, max=chr_size[chr])))
        coco_distance <- cos[2:length(cos)] - cos[1:(length(cos)-1)]
        mean_coco_distance <- mean(coco_distance)
		#sim.cotable <- tibble(lib=lib, chr=chr, cos=cos)
		return(mean_coco_distance)
}

# simulateCOCOdistPop: simulate CO-CO distance of the permuted cotable
# input = a cotable to make corresponding permuted cotable
simulateCOCOdistPop <- function(dat){
		dat <- dat %>%
				arrange(lib, chr)

		cotable.summary <- dat %>%
				group_by(lib, chr) %>%
				summarise(nco=n())
		f2 <- mean(apply(cotable.summary, 1, simulateCOCOdist), na.rm=T)
				#bind_rows() %>%
				#arrange(lib, chr, cos)
		return(f2)
}

# calcDCOdistance: calculate mean DCO distance of dco.table
calcDCOdistance <- function(dat){
        dat <- arrange(dat, lib, chr, cos)
#        dco_end
        dco_end <- filter(dat, row_number()%%2 == 0)
        dco_start <- filter(dat, row_number()%%2 == 1)
        dco_distance <- dco_end$cos - dco_start$cos
        dco_distance.mean <- mean(dco_distance)
		return(dco_distance.mean)
}

#========================
# calculate DCO distance of observed data and permuted data
#========================
obs <- recq4.dcotable
obs.distance <- calcDCOdistance(obs)


n_iter <- 2000
sim.distance.repeat <- c()
for(i in 1:n_iter){
    sim <- simulateCOCOdistPop(recq4.cotable)
    #sim <- simulateF2(obs)
    #sim.distance <- calcDCOdistance(sim)
    sim.distance.repeat <- c(sim.distance.repeat, sim)
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



#for(i in 1:nrow(recq4.cotable.summary)){
#    i <- 1
#    lib <- recq4.cotable.summary$lib[i]
#    chrs <- recq4.cotable.summary$chrs[i]
#    nchr <- as.numeric(str_replace(chrs, "Chr", ""))
#
#    random_cos <- runif(recq4.cotable.summary$ncos[i], min=1, max=chr.ends[nchr])
#
#    temp.cotable <- 
#    ctrl.cotable <- bind_rows(ctrl.cotable, )
#
#
#for(nlib in 1:length(recq4.cotable$lib)){
#    ctr.cotable.lib <- tibble()
#
#    for(nchr in 1:5){
#
#
#
#write_tsv(all.dco_table, file=file.path(dirout, paste0(libname, "_dco.smooth.txt")))
