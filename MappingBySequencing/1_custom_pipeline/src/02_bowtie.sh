#!/bin/bash

dirin="results/01_cutadapt"
dirbam="results/02_bowtie"
dirlog="log/bowtie2"

mkdir -p $dirin $dirbam $dirlog

# input
mut_1="results/01_cutadapt/hcr7_mut_p_re_1.tr.fastq.gz"
mut_2="results/01_cutadapt/hcr7_mut_p_re_2.tr.fastq.gz"
wt_1="results/01_cutadapt/hcr7_WT_p_re_1.tr.fastq.gz"
wt_2="results/01_cutadapt/hcr7_WT_p_re_2.tr.fastq.gz"

# output
wt_bam="results/02_bowtie/hcr7_wt.bam"
mut_bam="results/02_bowtie/hcr7_mut.bam"

# bowtie2 parameters
n_thread=30
tair10_fasta="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
tair10_index="/home/nison/work/refgenome/TAIR10/bowtie2_index_test"

#===========================
# 1. bowtie2
#===========================

# Build Bowtie2 index
if [[ ! -d $tair10_index ]]; then
	echo "building bowtie2 index"
	mkdir -p $tair10_index
	cd $tair10_index
	bowtie2-build --thread 5 $tair10_fasta tair10
fi

# Run bowtie2
(bowtie2 --very-sensitive --threads $n_thread -x $tair10_index/tair10 \
	-1 $wt_1 -2 $wt_2 \
	| samtools view -bh -@ $n_thread -F 2308 -o $wt_bam - ) 2> $dirlog/hcr7_wt_bowtie2.log

(bowtie2 --very-sensitive --threads $n_thread -x $tair10_index/tair10 \
	-1 $mut_1 -2 $mut_2 \
	| samtools view -bh -@ $n_thread -F 2308 -o $mut_bam - ) 2> $dirlog/hcr7_mut_bowtie2.log

# [samtools paramters explained]
# -F 2308 = filter out below
## read unmapped
## not primary alignment
## supplementary alignment

# -bh = output in bam format with header included

echo finished bowtie2

#===========================
# 2. Sorting and Filtering
#===========================

samtools sort -@ $n_thread -m 5G -o ${wt_bam%.bam}.sort.bam ${wt_bam} 2> $dirlog/wt_sort.log
samtools sort -@ $n_thread -m 5G -o ${mut_bam%.bam}.sort.bam ${mut_bam} 2> $dirlog/mut_sort.log

echo finished sorting
##============================
## 3. markdup
##============================

picard MarkDuplicates -I ${wt_bam%.bam}.sort.bam \
	-O ${wt_bam%.bam}.sort.md.bam \
	-M $dirlog/wt_md.log \
	--REMOVE_DUPLICATES true

picard MarkDuplicates -I ${mut_bam%.bam}.sort.bam \
	-O ${mut_bam%.bam}.sort.md.bam \
	-M $dirlog/mut_md.log \
	--REMOVE_DUPLICATES true

samtools index ${wt_bam%.bam}.sort.md.bam
samtools index ${mut_bam%.bam}.sort.md.bam
