#!/bin/bash

dirhome="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs"

# 1. Make sample map for GenomicsDBImport
dirvcf="$dirhome/result/02_haplotypecaller"
dirout="$dirhome/result/03_joint_genotyping"

mkdir -p $dirout/log

#for f in $dirvcf/*.md.rg.g.vcf.gz; do
#    index=$(basename $f)
#    index=${index%_MappedOn_tair10_sort.md.rg.g.vcf.gz}
#    index=${index#lib}
#    sample_name="wt20_${index}"
#    echo -e "$sample_name\t$f" >> $dirout/col_ler_f2.sample_map
#done

# 2. GenomicsDBImport

dir_db="$dirhome/result/03_joint_genotyping/genomics_db"
samplemap="$dirout/col_ler_f2.sample_map"
dir_interval="$dirhome/data"
interval_list_dbimport="$dir_interval/genomics_db.intervals"
interval_list_gvcfs="$dir_interval/intervals_for_genotypegvcfs.txt"


tair10="/datasets/data_3/genome_assemblies/TAIR10_Chr15.fasta"
tmpdir="$dirhome/tmp"
mkdir -p $tmpdir

#gatk --java-options "-Xmx20g -Xms20g" \
#		GenomicsDBImport \
#		--genomicsdb-workspace-path $dir_db \
#		--batch-size 50 \
#        -L $interval_list_dbimport \
#		--sample-name-map $samplemap \
#        --tmp-dir $tmpdir \
#        --max-num-intervals-to-import-in-parallel 30 \
#		--reader-threads 5 &> genomicsdbimport.log


# 3. Run GenotypeGVCFs in parallel
while read -r chr start end; do
    echo "start running interval ${chr}:${start}-${end}"
    gatk --java-options "-Xmx20g -Xms20g" GenotypeGVCFs \
        -R $tair10 \
        -V gendb://$dir_db \
        -L ${chr}:${start}-${end} \
        --only-output-calls-starting-in-intervals \
        -O $dirout/col_ler_f2_genotype_${chr}_${start}-${end}.vcf.gz &> $dirout/log/genotypegvcfs_${chr}_${start}-${end}.log &
done < $interval_list_gvcfs
wait;

echo finished cohort genotyping

