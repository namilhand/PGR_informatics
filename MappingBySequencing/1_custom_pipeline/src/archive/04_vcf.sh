#!/bin/bash

dirbam="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/02_bowtie"
ref="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
dirvcf="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/03_vcf"
dirlog=$dirvcf/log

mkdir -p $dirlog

#bcftools mpileup -f $ref -Ou -a FORMAT/AD,FORMAT/DP \
#    $dirbam/hcr7_WT_p_re.sort.md.bam $dirbam/hcr7_mut_p_re.sort.md.bam \
#    | bcftools call -mv -Ov -V indels -o $dirvcf/hcr7.vcf 2> $dirlog/hcr7_mpileup.log
bcftools mpileup -f $ref -Ou -a FORMAT/AD,FORMAT/DP \
    $dirbam/hcr7_WT_p_re.sort.md.bam $dirbam/hcr7_mut_p_re.sort.md.bam \
    | bcftools call -mv -Ov -o $dirvcf/hcr7_allvar.vcf 2> $dirlog/hcr7_allvar_mpileup.log
