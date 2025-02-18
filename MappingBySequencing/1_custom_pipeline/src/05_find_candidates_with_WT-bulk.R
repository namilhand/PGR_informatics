library(tidyverse)

# Set the arguments below

allvar=read_tsv("../results/04_annotate_vcf/hcr7_allvar.se.tsv", col_names=T)
dirout="../results/05_find_candidates"
dir.create(dirout, recursive=T)

#=============================
# 1. Draw landscape
#=============================
markers <- filter(allvar, (ref == "C" & alt == "T") | (ref == "G" & alt == "A")) %>%
	filter(chr != "ChrM", chr != "ChrC") %>%
	filter(wt.dp > 8 & mut.dp > 8) %>%
	filter(wt.alt != 0, mut.alt !=0)

pMarker <- ggplot(markers, aes(x=pos, y=ratio)) +
	geom_point(size=0.6) +
	facet_grid(rows=vars(chr), scales="free_x") +
	geom_smooth(colour="blue", method="loess", span=0.1, linewidth=0.5) +
	scale_y_continuous(limits=c(-1, 1)) +
	scale_x_continuous(labels= ~ round(.x/1000000, digits=0)) +
	theme_bw() +
	labs(x="Position (Mb)", y="ratio")

pdf(file=file.path(dirout, paste0("marker_ratio_landscape.pdf")), width=4, height=4)
print(pMarker)
dev.off()

png(file=file.path(dirout, paste0("marker_ratio_landscape.png")), width=4, height=4, res=300, unit="in")
print(pMarker)
dev.off()


# Causal locus
# Chr3: 18,000,000 - End


#=============================
# 2. Find candidate
#=============================

effect <- c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "stop_lost", "start_lost", "stop_gained", "missense_variant", "coding_sequence_variant", "inframe_insertion", "disruptive_inframe_insertion", "inframe_deletion", "disruptive_inframe_deletion", "exon_variant", "exon_loss_variant", "duplication", "inversion", "frameshift_variant", "frameshift_ablation", "gene_fusion", "bidirectional_gene_fusion", "rearranged_at_DNA_level", "miRNA", "initiator_codon_variant", "start_retained")

cand <- separate(allvar, col=mutation_effect, into=c("effect1", "effect2"), sep="&", fill="right") %>%
	filter(effect1 %in% effect | effect2 %in% effect) %>%
	filter(type == "INDEL" | (ref == "C" & alt == "T") | (ref == "G" & alt == "A")) %>%
	filter(chr == "Chr3" & pos > 18000000) %>%
	filter(ratio > 0.3)

write_tsv(cand, file.path(dirout, paste0("hcr7_candidate.tsv")))
