# pipeline

Objective: To retain highly reliable SNPs for GBS

Input data: bam files of WT GBS (2020)

# Process

## 1. Calling Syntenic SNPs using SyRi

script: 01_syri.sh

SyRi는 서로 다른 genome assembly를 비교해서 variant들을 찾아주는 프로그램이다. Col과 Ler assembly를 비교해 Syntenic region의 SNP들을 찾을 수 있다.

```bash
#!/bin/bash

# Replace 'wd' and 'thread' variable accordingly
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

```

`syri.out` 파일에는 두 assembly 사이의 variant들이 기록돼 있다. 이 중 Syntenic region의 SNP들만을 골라서 bed format으로 변환한다.

```bash

# 3. Extracting syntenic snps from syri.out
awk 'BEGIN{FS="\t"; OFS="\t"}substr($10, 1, 3)=="SYN" && $11=="SNP"{printf "%s\t%s\t%s\n", $1, $2-1, $2}' $dirout/syri.out > $dirout/synsnps.bed
```

## 2. HaplotypeCaller

script: 02_haplotypecaller.sh

SyRi를 이용해 찾은 syntenic region의 snp들 중에 genotyping marker로 신뢰할만한 SNP들을 추려내는 과정이 필요하다. Col x Ler F2 population에서 각 SNP의 segregation pattern등을 활용할 수 있다.

`bcftools mpileup` 과 마찬가지로 variant가 발견되는 position에서 sample의 genotype likelihood를 계산한다. bcftools 대신에 GATK HaplotypeCaller를 사용하는 이유는 GATK가 여러 sample의 genotype data들을 한번에 처리하기에 더 용이하기 때문이다. SNP filtering의 경우 수십-수백개 sample의 genotype 정보를 동시에 처리해야 하기 때문에 BCFtools보다 GATK를 사용하는 것이 더 수월하다.

### 1) AddOrReplaceReadGroups

GATK는 bam file에 RG (ReadGroup) tag가 달려있는 것을 가정한다. BWA를 사용할 경우 RG tag이 기본으로 달리지만 bowtie2는 그렇지 않다. 따라서 AddOrReplaceGroup을 이용해 RG tag을 따라 달아줘야한다.

RGID부터 RGSM까지의 tag들은 모두 하나의 read가 어느 sample, 어느 library, 어느 sequencing platform 에서 왔는지 등을 표시한다. 아무 의미 없는 문자열을 달아주면 된다.

```bash
#!/bin/bash

export rgpu="HHKJ7CCX2"
export dirbam="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/input"

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

parallel --jobs 32 addrg {} ::: $dirbam/*.md.bam

```

(addrg, index) haplotype caller

### 2) HaplotypeCaller

HaplotypeCaller를 돌리기 위해서는 bam파일들이 indexing되어 있어야 한다. `samtools index`를 이용할 수 있다.

GATK tool들은 reference genome의 fasta파일과 `.dict` 로 끝나는 dictionary file을 요구한다. `.dict`가 없다면 `gatk CreateSequenceDictionary` 명령어로 생성한다.

```bash
export ref="/datasets/data_3/genome_assemblies/TAIR10_Chr15.fasta"
export synsnps="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/result/01_syri/synsnps.bed"
export dirvcf="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/result/02_haplotypecaller"
mkdir -p $dirvcf

parallel --jobs 35 samtools index {} ::: $dirbam/*.md.rg.bam
# HaplotypeCaller requires bam files to be indexed

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
        
    # Create dictionary for reference sequence
    if [[ ! -f ${ref%.fasta}.dict ]]; then
		    echo ">>>>>>>>> CREATING FASTA DICTIONARY FILE <<<<<<<"
		    gatk CreateSequenceDictionary -R $ref
		    
		# Check if a bam file is already genotyped. If so, skip that file.
		if [[ -f $dirout/${filename%.bam}.g.vcf.gz.tbi ]]; then
				echo ">>>>>>>> skip $index <<<<<<<<"
				exit;
		else
				echo "======== run $index ========"
		fi
		
		# HaplotypeCaller
		gatk --java-options "-Xmx1g -Djava.io.tmpdir=$dirtmp"\
				HaplotypeCaller \
				-R $ref \
				--intervals $synsnps \
				-ERC GVCF \
				-I $input \
				-O $dirout/${filename%.bam}.g.vcf.gz &> $dirlog/${index}.g.vcf.log
		
		# [HaplotypeCaller options]
		# -ERC GVCF: output in genome-wide VCF format which contains information about
		# not only variant position but also non-variant position
		
		# --intervals: Genotype only selected intervals. Here, we only genotype the Col/Ler Syntenic SNP positions
		
		# -Xmx1g: (java option) set maximum Java heap size as 4G
		
		# -Djava.io.tmpdir=$tmpdir: (java option; very important!) default temporary directory is /tmp which has very limited space in our server.
		# It should be set in somewhere having enough space.
		}
export -f hc

parallel --jobs 25 hc {} $dirvcf ::: $dirbam/*.md.rg.bam
```

## joint genotyping

script: 03_joint_genotyping.sh

**GATK tools used in this script**

1. `gatk GenomicsDBImport`: 여러 sample들의 VCF 파일들에 대한 DB를 만든다. 추후에 sample이 더 늘어났을 때 DB에 추가할 수 있다. 이 때 sample map이 필요하다.
2. `gatk GenotypeGVCFs`: DB에 등록된 VCF 파일들에 대해 genotyping을 수행한다.

### 1) Preparing sample map

```bash
#!/bin/bash

dirhome="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs"

# 1. Make sample map for GenomicsDBImport
dirvcf="$dirhome/result/02_haplotypecaller"
dirout="$dirhome/result/03_joint_genotyping"

mkdir -p $dirout/log

for f in $dirvcf/*.md.rg.g.vcf.gz; do
    index=$(basename $f)
    index=${index%_MappedOn_tair10_sort.md.rg.g.vcf.gz}
    index=${index#lib}
    sample_name="wt20_${index}"
    echo -e "$sample_name\t$f" >> $dirout/col_ler_f2.sample_map
done

```

**example output**

### 2) GenomicsDBImport

```bash
# 2. GenomicsDBImport

dir_db="$dirhome/result/03_joint_genotyping/genomics_db"
samplemap="$dirout/col_ler_f2.sample_map"
dir_interval="$dirhome/data"
interval_list_dbimport="$dir_interval/genomics_db.intervals"
interval_list_gvcfs="$dir_interval/intervals_for_genotypegvcfs.txt"

tair10="/datasets/data_3/genome_assemblies/TAIR10_Chr15.fasta"
tmpdir="$dirhome/tmp"
mkdir -p $tmpdir

gatk --java-options "-Xmx20g -Xms20g" \
		GenomicsDBImport \
		--genomicsdb-workspace-path $dir_db \
		--batch-size 50 \
        -L $interval_list_dbimport \
		--sample-name-map $samplemap \
        --tmp-dir $tmpdir \
        --max-num-intervals-to-import-in-parallel 30 \
		--reader-threads 5 &> genomicsdbimport.log

```

### 3) GenotypeGVCFs

```bash
# 3. Run GenotypeGVCFs in parallel
while read -r chr start end; do
    echo "start running interval ${chr}:${start}-${end}"
    gatk --java-options "-Xmx20g -Xms20g" GenotypeGVCFs \
        -R $tair10 \
        -V gendb://$ws \
        -L ${chr}:${start}-${end} \
        --only-output-calls-starting-in-intervals \
        -O $dirout/col_ler_f2_genotype_${chr}_${start}-${end}.vcf.gz &> $dirout/log/genotypegvcfs_${chr}_${start}-${end}.log &
done < $interval_list_gvcfs
wait;

echo finished cohort genotyping
```

### Filtering variants

### Select variants

### result: Filtered syntenic snps

## References

[1] minimap2: https://lh3.github.io/minimap2/minimap2.html

[2] syri: [https://schneebergerlab.github.io/syri/](https://schneebergerlab.github.io/syri/)

[3] plotsr: https://schneebergerlab.github.io/syri/plotsr.html