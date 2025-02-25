#!/bin/bash

#export srcPlot="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/src/draw_plot.R"
#export dirSnpEff="/home/nison/src/snpEff"
#export dirVcf="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/03_vcf"
#export dirout="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/04_result"

vcf="results/03_vcf/hcr7_allvar.vcf"
dirout="results/04_annotate_vcf"

# Does vcf file contain the WT_bulk information?
# yes or no
wt_bulk="yes"

mkdir -p $dirout


snpEff Arabidopsis_thaliana -s $dirout/hcr7_allvar.summary.html $vcf > $dirout/hcr7_allvar.se.vcf

if [[ $wt_bulk == "yes" ]]; then
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "type" "wt.ref" "wt.alt" "mut.ref" "mut.alt" "wt.dp" "mut.dp" "ratio" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" > $dirout/hcr7_allvar.se.tsv

	bcftools query $dirout/hcr7_allvar.se.vcf \
			-f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/ANN[\t%AD][\t%DP]\n' |
	awk 'BEGIN{IFS="\t" ; OFS="\t";}{print}' | tr -s "\t" " " | \
	awk 'BEGIN{FS=" "; OFS=" "}{split($7, ad_wt, ","); split($8, ad_mut, ","); ratio=(ad_wt[1]/(ad_wt[1] + ad_wt[2] + 0.001) - (ad_mut[1])/(ad_mut[1] + ad_mut[2] + 0.001)); print $0, ratio}' | \
	awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,wt_ad,","); split($8,mut_ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, $5, wt_ad[1], wt_ad[2], mut_ad[1], mut_ad[2], $9, $10, $11, a[2], a[4], a[5], a[10], a[11]}' >> $dirout/hcr7_allvar.se.tsv
elif [[ $wt_bulk == "no" ]]; then
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "type" "mut.ref" "mut.alt" "mut.dp" "ratio" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" > $dirout/hcr7_allvar.se.tsv

	bcftools query $dirout/hcr7_allvar.se.vcf \
			-f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/ANN[\t%AD][\t%DP]\n' |
	awk 'BEGIN{IFS="\t" ; OFS="\t";}{print}' | tr -s "\t" " " | \
	awk 'BEGIN{FS=" "; OFS=" "}{split($7, ad_mut, ","); ratio=ad_mut[1]/(ad_mut[1] + ad_mut[2] + 0.001); print $0, ratio}' | \
	awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,mut_ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, $5, mut_ad[1], mut_ad[2], $8, $9, a[2], a[4], a[5], a[10], a[11]}' >> $dirout/hcr7_allvar.se.tsv
else
	print "wt_bulk argument should be either yer or no"
fi

