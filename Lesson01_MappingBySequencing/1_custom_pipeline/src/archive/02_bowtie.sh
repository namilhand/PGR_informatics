#!/bin/bash

export dirin="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/01_cutadapt"
export dirbam="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/02_bowtie"

export mut_1="hcr7_mut_p_re_1.tr.fastq.gz"
export mut_2="hcr7_mut_p_re_2.tr.fastq.gz"
export wt_1="hcr7_WT_p_re_1.tr.fastq.gz"
export wt_2="hcr7_WT_p_re_2.tr.fastq.gz"

export n_thread=30

export ref="/home/nison/work/refgenome/TAIR10/bowtie2_index/tair10"

#===========================
# 1. bowtie2
#===========================

function bt {
    input_1=$1
    input_2=$2
    dirout=$3

    output=$(basename $input_1)
    output=${output%_1.tr.fastq.gz}.bam

    log=${output%.bam}.log
    dirlog=$dirout/log

    mkdir -p $dirout $dirlog

    (bowtie2 --very-sensitive \
        --threads $n_thread \
        -x $ref -1 $input_1 -2 $input_2 \
        | samtools view -bh -@ $n_thread -F 2308 -o $dirout/$output - ) 2> $dirlog/$log
    }

export -f bt

#bt $dirin/$wt_1 $dirin/$wt_2 $dirbam
#bt $dirin/$mut_1 $dirin/$mut_2 $dirbam
#wait;
echo finished bowtie2

#===========================
# 2. filtering
#===========================

function sortbam {
    input=$1
    dirout="$dirbam"
    output=$(basename $input)
    output=${output%.bam}.sort.bam
    
    dirlog=$dirout/log
    log=${output%.bam}.log

    mkdir -p $dirlog

    #(samtools view -h $input -q 2 -u \
    #    | samtools sort -@ $n_thread -m 5G -o $dirout/$output - ) 2> $dirlog/$log

    samtools sort -@ $n_thread -m 5G -o $dirout/$output $input 2> $dirlog/log

    }

export -f sortbam

#sortbam $dirbam/hcr7_WT_p_re.bam
#sortbam $dirbam/hcr7_mut_p_re.bam
echo finished sorting

#============================
# 3. markdup
#============================

function markdup {
    input=$1
    dirout="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/results/02_bowtie"
    output=$(basename $input)
    output=${output%.bam}.md.bam
    index=${output}.bai

    dirlog=$dirout/log
    output_metric=${output%.bam}.metric.txt

    mkdir -p $dirlog

    picard MarkDuplicates -I $input \
        -O $dirout/$output \
        -M $dirlog/$output_metric \
        --REMOVE_DUPLICATES true;
    
    samtools index $dirout/$output -o $dirout/$index
    }

export -f markdup

markdup $dirbam/hcr7_WT_p_re.sort.bam
markdup $dirbam/hcr7_mut_p_re.sort.bam
