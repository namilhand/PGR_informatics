#!/bin/bash

wd="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/result"
dirout="$wd/01_syri"
mkdir -p $dirout

tair10="/datasets/data_3/genome_assemblies/TAIR10_Chr15.fasta"
lerhifi="/datasets/data_3/genome_assemblies/Ler-HiFi_Chr15.fasta"

thread=20



# 1. align Ler-HiFi to Col-Hifi
minimap2 -ax asm5 -g 100 --eqx -t $thread $tair10 $lerhifi | samtools sort -O BAM - > $dirout/lerhifi_to_tair10.bam \
	&& samtools index $dirout/lerhifi_to_tair10.bam

## minimap2 -a -a [-x preset] target.fa query.fa > output

# 2. Run SyRi
syri -c $dirout/lerhifi_to_tair10.bam -r $tair10 -q $lerhifi -k -F B --dir "$dirout"
## -k : keep intermeidate output files
## -F : input file type. S: SAM

# 3. Extracting syntenic snps from syri.out
awk 'BEGIN{FS="\t"; OFS="\t"}substr($10, 1, 3)=="SYN" && $11=="SNP"{printf "%s\t%s\t%s\n", $1, $2-1, $2}' $dirout/syri.out > $dirout/synsnps.bed

