#!/bin/bash

export synsnps="result/01_syri/synsnps.bed"
export ref="/datasets/data_3/genome_assemblies/TAIR10_Chr15.fasta"
export rgpu="HHKJ7CCX2"
export dirbam="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/input"
export dirvcf="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/result/02_haplotypecaller"

mkdir -p $dirvcf

function addrg {
    input=$1
    file=$(basename $input)
    sample=${file%_MappedOn_tair10_sort.md.bam}
    gatk AddOrReplaceReadGroups \
        I=$input \
        O=${input%.bam}.rg.bam \
        RGID=$sample \
        RGLB=$sample \
        RGPL=ILLUMINA \
        RGPU=$rgpu \
        RGSM=$sample
    }
export -f addrg


function hc {
		input=$1
		filename=$(basename $input)
		index=${filename%_MappedOn_tair10_sort.md.rg.bam}
        dirout=$2

        dirlog=$dirout/log
        dirtmp=$dirout/tmp

        mkdir -p $dirout
        mkdir -p $dirlog
        mkdir -p $dirtmp

		if [[ ! -f ${ref%.fasta}.dict ]]; then
			echo ">>>>>>>>> CREATING FASTA DICTIONARY FILE <<<<<<<"
			gatk CreateSequenceDictionary -R $ref
		fi

		if [[ -f $dirout/${filename%.bam}.g.vcf.gz.tbi ]]; then
				echo ">>>>>>>> skip $index <<<<<<<<"
				exit;
		else
				echo "======== run $index ========"
		fi
		gatk --java-options "-Xmx1g -Djava.io.tmpdir=$dirtmp"\
				HaplotypeCaller \
				-R $ref \
				--intervals $synsnps \
				-ERC GVCF \
				-I $input \
				-O $dirout/${filename%.bam}.g.vcf.gz &> $dirlog/${index}.g.vcf.log
		}
export -f hc

#parallel --jobs 32 addrg {} ::: $dirbam/*.md.bam
#parallel --jobs 35 samtools index {} ::: $dirbam/*.md.rg.bam
parallel --jobs 25 hc {} $dirvcf ::: $dirbam/*.md.rg.bam



# -Xmx4g option in java: set maximum Java heap size as 4g
# --interval-padding100: added 100 bp padding to interval so that the caller sees enough context to reevaluate the call appropriately
# --spark-master local[32]: run on the local machine using 32 cores
