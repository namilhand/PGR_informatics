#!/bin/bash

export mut_1="raw/hcr7_mut_p_re_1.fastq.gz"
export mut_2="raw/hcr7_mut_p_re_2.fastq.gz"
export wt_1="raw/hcr7_WT_p_re_1.fastq.gz"
export wt_2="raw/hcr7_WT_p_re_2.fastq.gz"

# hcr7 mapping libraries were constructed using GBS protocol.
# The adapters used for libraries are:
# wt: P1=88
# mut: P1=84

export read1_tr="AGAGGATAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGA"
# reverse complement of P2 adapter sequence 

export wt_read2_tr="AGATAGAGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG"
export mut_read2_tr="AGTAAGGTGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG"
# reverse complement of P1 adapter sequence
# WT = P1-88
# mut = P1-84

export dirout="results/01_cutadapt"
export dirqc="qc"
export dirlog="log"

mkdir -p $dirout $dirqc $dirlog


function trimming {
    input_1=$1
    input_2=$2
	dirout=$3
	read2_tr=$4

    input_name_1=$(basename $input_1)
    input_name_2=$(basename $input_2)

    output_name_1=${input_name_1%.fastq.gz}.tr.fastq.gz
    output_name_2=${input_name_2%.fastq.gz}.tr.fastq.gz

	qcfile=${intput_name_1%.fastq.gz}.cutadapt.qc.txt
	logfile=${intput_name_1%.fastq.gz}.cutadapt.log

    cutadapt -u 9 \
        -U 9 \
        -a $read1_tr \
        -A $read2_tr \
        -O 3 \
        -q 19 \
        $input_1 $input_2 \
        -o $dirout/$output_name_1 \
        -p $dirout/$output_name_2 > $dirqc/$qcfile 2> $dirlog/$logfile

	# -u: Remove the first N bases from Read1
	# -U: Remove the first N bases from Read2
	# -a: Regular 3' adapter of Read1
	# -A: Regular 3' adapter of Read2
	# -O: minimum overlap
	# -q: base quality filter
}

export -f trimming

trimming $wt_1 $wt_2 $dirout $wt_read2_tr &
trimming $mut_1 $mut_2 $dirout $mut_read2_tr
wait;
echo finished trimming
