#!/bin/bash

dirbam="results/02_bowtie"
ref="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
dirvcf="results/03_vcf"
dirlog="$dirvcf/log"

mkdir -p $dirlog

bcftools mpileup -f $ref -Ou -a FORMAT/AD,FORMAT/DP \
    $dirbam/hcr7_wt.sort.md.bam $dirbam/hcr7_mut.sort.md.bam \
    | bcftools call -mv -Ov -o $dirvcf/hcr7_allvar.vcf 2> $dirlog/hcr7_allvar_mpileup.log

# bcftools mpileup: Generate VCF/BCF containing genotype likelihoods for one or multiple alignment files.
## -Ou: output in uncompressed format --> Stream the output to the next pipe
## -a: annotate
## FORMAT/AD: Allelic depth
## FORMAT/DP: Number of high-quality bases

# bcftools call: Calling SNP/indel
# -m: default variant calling model.
# -v: variants-only
# -Ov: output in vcf format
