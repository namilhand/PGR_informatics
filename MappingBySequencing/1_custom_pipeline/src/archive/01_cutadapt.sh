#!/bin/bash

mkdir -p "/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results"
export dirraw="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/raw"
export dirout="/datasets/data_4/nison/NGS_libraries/hcr7_mapping/202408/results/01_cutadapt"

export mut_1="hcr7_mut_p_re_1.fastq.gz"
export mut_2="hcr7_mut_p_re_2.fastq.gz"
export wt_1="hcr7_WT_p_re_1.fastq.gz"
export wt_2="hcr7_WT_p_re_2.fastq.gz"

export n_thread=30

#===========================
# 1. trimming
#===========================

function trimming {
    input_1=$1
    input_2=$2
    wtmut=$3
    dir_out=$4
    input_name_1=$(basename $input_1)
    input_name_2=$(basename $input_2)
    output_name_1=${input_name_1%.fastq.gz}.tr.fastq.gz
    output_name_2=${input_name_2%.fastq.gz}.tr.fastq.gz

    dirqc=$dir_out/qc
    qcfile=${input_name_1%.fastq.gz}.cutadapt.qc.txt
    logfile=${input_name_1%.fastq.gz}.cutadapt.log

    mkdir -p $dir_out $dirqc


    # wt: P1=88
    # mut: P1=84
    if [[ $wtmut == "wt" ]]; then
        read2_tr=AGATAGAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG
    elif [[ $wtmut == "mut" ]]; then
        read2_tr=AGTAAGGTGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG
    else
        echo wtmut should be either wt or mut
    fi

    

    cutadapt -u 9 \
        -U 9 \
        -a AGAGGATAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA \
        -A $read2_tr \
        -O 3 \
        -q 19 \
        $input_1 $input_2 \
        -o $dir_out/$output_name_1 \
        -p $dir_out/$output_name_2 > $dirqc/$qcfile 2> $dirqc/$logfile
    }

export -f trimming

trimming $dirraw/$wt_1 $dirraw/$wt_2 wt $dirout &
trimming $dirraw/$mut_1 $dirraw/$mut_2 mut $dirout
wait;
echo finished trimming

