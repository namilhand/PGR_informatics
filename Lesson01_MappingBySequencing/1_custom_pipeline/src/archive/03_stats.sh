#!/bin/bash

dirbam="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/02_bowtie"
dirstat=$dirbam/stats

mkdir -p $dirstat

samtools idxstats $dirbam/hcr7_WT_p_re.sort.bam > $dirstat/hcr7_WT_p_re.sort.idxstat &
samtools idxstats $dirbam/hcr7_mut_p_re.sort.bam > $dirstat/hcr7_mut_p_re.sort.idxstat
