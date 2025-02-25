#!/bin/bash

export srcPlot="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/src/draw_plot.R"
export dirSnpEff="/home/nison/src/snpEff"
export dirVcf="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/03_vcf"
export dirout="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/04_result"

function returnCand {
	input=$1
	dirout=$2
	base=$(basename $input)
	prefix=${base%.vcf}

	#java -Xmx4g -jar $dirSnpEff/snpEff.jar Arabidopsis_thaliana -c $dirSnpEff/snpEff.config -s $dirout/${prefix}.summary.html $input > $dirout/${prefix}.se.vcf
	# #-Xmx: maximum heap size

	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "wt.ref" "wt.alt" "mut.ref" "mut.alt" "wt.dp" "mut.dp" "ratio" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" > $dirout/${prefix}.se.tsv

	bcftools query $dirout/${prefix}.se.vcf \
			-f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/ANN[\t%AD][\t%DP]\n' |
	awk 'BEGIN{IFS="\t" ; OFS="\t";}{print}' | tr -s "\t" " " | \
	awk 'BEGIN{FS=" "; OFS=" "}{split($7, ad_wt, ","); split($8, ad_mut, ","); ratio=((ad_wt[1] + 0.1)/(ad_wt[1] + ad_wt[2] + 0.1) - (ad_mut[1] + 0.1)/(ad_mut[1] + ad_mut[2] + 0.1)); print $0, ratio}' | \
	awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,wt_ad,","); split($8,mut_ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, wt_ad[1], wt_ad[2], mut_ad[1], mut_ad[2], $9, $10, $11, a[2], a[4], a[5], a[10], a[11]}' >> $dirout/${prefix}.se.tsv


#	bcftools query $dirout/${prefix}.se.vcf \
#			-f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/ANN[\t%AD][\t%DP]\n' |
#			-o $dirout/${prefix}.se.table
#
#
#	awk 'BEGIN{IFS="\t" ; OFS="\t";}{print}' $dirout/${prefix}.se.table | tr -s "\t" " " | \
#	awk 'BEGIN{FS=" "; OFS=" "}{split($7, ad_wt, ","); split($8, ad_mut, ","); ratio=((ad_wt[1] + 0.1)/(ad_wt[1] + ad_wt[2] + 0.1) - (ad_mut[1] + 0.1)/(ad_mut[1] + ad_mut[2] + 0.1)); print $0, ratio}' | \
#	awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,wt_ad,","); split($8,mut_ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, wt_ad[1], wt_ad[2], mut_ad[1], mut_ad[2], $9, $10, $11, a[2], a[4], a[5], a[10], a[11]}' >> $dirout/${prefix}.bulkedSeg.allvar.tsv
#
    Rscript $srcPlot $dirout/${prefix}.se.tsv $dirout

	#awk 'BEGIN{IFS=" "; OFS=" ";}/splice_acceptor_variant|splice_donor_variant|splice_region_variant|stop_lost|start_lost|stop_gained|missense_variant|coding_sequence_variant|inframe_insertion|disruptive_inframe_insertion|inframe_deletion|disruptive_inframe_deletion|exon_variant|exon_loss_variant|duplication|inversion|frame_shift_variant|frameshift_ablation|gene_fusion|bidirectional_gene_fusion|rearranged_at_DNA_level|miRNA|initiator_codon_variant|start_retained/' $dirout/${prefix}.se.ctga.table |\
	#awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,wt.ad,","); split($8,mut.ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, $9, $10, $11, a[2], a[4], a[5], a[10], a[11]}' >> ${prefix}.bulkedSeg.result
	#Rscript $srcPlot $dirout/${prefix}.se.ctga.table $dirout/${prefix}.svg
	}

export -f returnCand
mkdir -p $dirout
returnCand $dirVcf/hcr7_allvar.vcf $dirout

