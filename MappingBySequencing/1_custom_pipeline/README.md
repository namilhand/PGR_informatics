# custom pipeline

## Setting environment

```bash
mamba env create --name MapBySeq --file=environment.yaml

conda activate MapBySeq
```

## Input file

/datasets/data_3/PGR_NGS_backup/2024/20240815_GBS_VIP3_HSC6_hcr7_HCR1_126_s_CO_m/RawData/hcr7_{mut/WT}_p_re_{1/2}.fastq.gz

hcr7_mut_p_re_1/2.fastq.gz: bulk sequencing of pooled F2 plants with mutant phenotype.

hcr7_WT_p_re_1/2.fastq.gz: bulk sequencing of pooled F2 plants with WT phenotype.

WT bulk sequencing 결과는 있으면 좋고 없어도 괜찮다.

## Process

### 1. Trimming reads

Script: src/01_cutadapt.sh

```bash
**#!/bin/bash

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
# reverse complement of P1 adapter sequence (See Figure 1)
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
echo finished trimming**
```

![**Figure 1. 현재 사용중인 dual-indexed GBS library.** P1/P2로 PCR후 Read1/2 sequencing primer로 sequencing한다. sequencing 결과 첫 9-bp는 index sequence이다.](https://github.com/namilhand/PGR_informatics/blob/main/MappingBySequencing/1_custom_pipeline/images/image_1.png)

**Figure 1. 현재 사용중인 dual-indexed GBS library.** P1/P2로 PCR후 Read1/2 sequencing primer로 sequencing한다. sequencing 결과 첫 9-bp는 index sequence이다.

Library제작에 사용된 Kit와 adapter 정보를 확인해서 trimming에 사용되는 sequence를 지정한다. Illumina library 의 경우 아래 페이지에서 trimming할 adapter sequence를 알 수 있다.

https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html

### 2. Align to TAIR10 and filtering bam file

Script: src/02_bowtie.sh

Note:

- You need to prepare TAIR10.fasta file.

```bash
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

# ==samtools paramters explained - See ref [4]==
# -F: filter by FLAG. What is FLAG? See ref [3]
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

# -m 5G: allocate 5G of memory to this task

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
	
# Without "--REMOVE_DUPLICATES true", duplicated alignment will not be discarded but marked with MD tag

samtools index ${wt_bam%.bam}.sort.md.bam
samtools index ${mut_bam%.bam}.sort.md.bam
```

### 3. Variant calling

- Bam file로부터 variant들을 불러내고 variant position에서 발견되는 allele들의 depth를 측정한다.
- Script: src/03_vcf.sh

- WT bulk에 대한 sequencing결과가 없다면 `mpileup`에서 hcr7_wt.sort.md.bam을 빼고 실행한다.
- `bcftools mpileup` → `bcftools call` 에 여러 BAM file을 한번에 넣어줄 수 있다. 각 BAM file에 대한  genotype 관련 정보 중 어떤 것들을 출력할지는 `-a` 를 이용해서 지정할 수 있다.

```bash
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
```

**Output example**

![**Figure 2. hcr7_allvar.vcf**](https://github.com/namilhand/PGR_informatics/blob/main/MappingBySequencing/1_custom_pipeline/images/image_2.png)

- FORMAT은 각 bamfile에 대한 annotation field (마지막 두 field)가 어떤 형식으로 기록되어 있는지를 보여준다.
- 넣어준 input bam file들의 genotyping 정보를 각 genomic position마다 FORMAT의 형식으로 보여준다.

### 4. Annotating variants

- src/4_annotate_variants.sh

snpEff는 variant마다 gene expression이나 protein sequence에 어떤 영향을 미칠 수 있는지를 예측해준다.
https://pcingola.github.io/SnpEff/


```bash
#!/bin/bash

vcf="results/03_vcf/hcr7_allvar.vcf"
dirout="results/04_annotate_vcf"

# Does vcf file contain the WT_bulk information?
# yes or no
wt_bulk="yes"

mkdir -p $dirout

snpEff Arabidopsis_thaliana -s $dirout/hcr7_allvar.summary.html $vcf > $dirout/hcr7_allvar.se.vcf
# Annotate each variants with snpEff.
# snpEff predicts the effect of variants on genes

if [[ $wt_bulk == "yes" ]]; then
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "wt.ref" "wt.alt" "mut.ref" "mut.alt" "wt.dp" "mut.dp" "ratio" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" > $dirout/${prefix}.se.tsv

	bcftools query $dirout/hcr7_allvar.se.vcf \
			-f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/ANN[\t%AD][\t%DP]\n' |
	awk 'BEGIN{IFS="\t" ; OFS="\t";}{print}' | tr -s "\t" " " | \
	awk 'BEGIN{FS=" "; OFS=" "}{split($7, ad_wt, ","); split($8, ad_mut, ","); ratio=(ad_wt[1]/(ad_wt[1] + ad_wt[2] + 0.001) - (ad_mut[1])/(ad_mut[1] + ad_mut[2] + 0.001)); print $0, ratio}' | \
	awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,wt_ad,","); split($8,mut_ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, wt_ad[1], wt_ad[2], mut_ad[1], mut_ad[2], $9, $10, $11, a[2], a[4], a[5], a[10], a[11]}' >> $dirout/hcr7_allvar.se.tsv
elif [[ $wt_bulk == "no" ]]; then
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "mut.ref" "mut.alt" "mut.dp" "ratio" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" > $dirout/${prefix}.se.tsv

	bcftools query $dirout/hcr7_allvar.se.vcf \
			-f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/ANN[\t%AD][\t%DP]\n' |
	awk 'BEGIN{IFS="\t" ; OFS="\t";}{print}' | tr -s "\t" " " | \
	awk 'BEGIN{FS=" "; OFS=" "}{split($7, ad_mut, ","); ratio=ad_mut[1]/(ad_mut[1] + ad_mut[2] + 0.001); print $0, ratio}' | \
	awk 'BEGIN{IFS=" "; OFS="\t"} {split($6,a,"|"); split($7,mut_ad,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, mut_ad[1], mut_ad[2], $8, $9, a[2], a[4], a[5], a[10], a[11]}' >> $dirout/hcr7_allvar.se.tsv
else
	print "wt_bulk argument should be either yer or no"
fi
```
**Output example**

![**Figure 3. hcr7_allvar.se.vcf**](https://github.com/namilhand/PGR_informatics/blob/main/MappingBySequencing/1_custom_pipeline/images/image%203.png)

**Figure 3. hcr7_allvar.se.vcf**

wt.ref: reference count of wt bulk

wt.alt: alternative count of wt bulk

mut.ref: reference count of mut bulk

mut.alt: alternative count of mut bulk

wt.dp: reliable depth of the position

mut.dp: reliable depth of the position

**Calculation of ratio**

ratio: {mut.alt/(mut.alt + mut.ref)} - {wt.alt/(wt.alt + wt.ref)} 

- If a variant is the recessive causal mutation, expected value is …
    - mut.alt/(mut.alt + mut.ref) → 1
    - wt.alt/(wt.alt + wt.ref) → 1/3
    - **Therefore, the ideal ratio is 2/3 (~0.66)**
    - **In this way you can exclude background mutations.**

Otherwise, you can just use allele frequency of mutant bulk. In that case,

- ratio = mut.alt/(mut.ref + mut.alt) ~ 1

### 5. Find candidates

- Script: src/05_find_candidates_[with/without]_WT-bulk.R
    - WT-bulk data가 있는 경우 → 05_find_candidates_with_WT-bulk.R
- 아래에는 WT-bulk data가 있는 경우를 기준으로 설명한다.

1. Setting input and output directory

```r
library(tidyverse)

# Set the arguments below

allvar=read_tsv("../results/04_annotate_vcf/hcr7_allvar.se.tsv", col_names=T)
dirout="../results/05_find_candidates"
dir.create(dirout, recursive=T)
```

1. Find locus

Reliable한 Marker를 찾고 marker들의 alternative allele frequency가 높은 지역을 찾는다.  

Reliable marker:

1. C>T or G>A mutation
2. Not in organell genome
3. Higher depth (wt.dp and mut.dp) ensure less variation
4. wt.alt(WT bulk의 alternative allele)이나 mut.alt(mutant bulk의 alternative allele)의 depth가 0일 수 없다.

sequencing depth에 따라 wt.dp, mut.dp 값은 조정이 필요하다.

```r
#=============================
# 1. Draw landscape
#=============================
markers <- filter(allvar, (ref == "C" & alt == "T") | (ref == "G" & alt == "A")) %>%
	filter(chr != "ChrM", chr != "ChrC") %>%
	filter(wt.dp > 8 & mut.dp > 8) %>%
	filter(wt.alt != 0, mut.alt !=0)

pMarker <- ggplot(markers, aes(x=pos, y=ratio)) +
	geom_point(size=0.6) +
	facet_grid(rows=vars(chr), scales="free_x") +
	geom_smooth(colour="blue", method="loess", span=0.1, linewidth=0.5) +
	scale_y_continuous(limits=c(-1, 1)) +
	scale_x_continuous(labels= ~ round(.x/1000000, digits=0)) +
	theme_bw() +
	labs(x="Position (Mb)", y="ratio")

pdf(file=file.path(dirout, paste0("marker_ratio_landscape.pdf")), width=4, height=4)
print(pMarker)
dev.off()

png(file=file.path(dirout, paste0("marker_ratio_landscape.png")), width=4, height=4, res=300, unit="in")
print(pMarker)
dev.off()

```

**Output example**

![**Figure 4. Alternative allele frequencies of EMS marker positions. Chr3의 south arm 끝 부분이 causal locus이다.** Chr3: 18 Mb-End 내에서 candidate mutation을 찾아보자.]([image%203.png](https://github.com/namilhand/PGR_informatics/blob/main/MappingBySequencing/1_custom_pipeline/images/image%204.png))

**Figure 4. Alternative allele frequencies of EMS marker positions. Chr3의 south arm 끝 부분이 causal locus이다.** Chr3: 18 Mb-End 내에서 candidate mutation을 찾아보자.

1. Find candidates

```r

# Causal locus
# Chr3: 18,000,000 - End

#=============================
# 2. Find candidate
#=============================

effect <- c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "stop_lost", "start_lost", "stop_gained", "missense_variant", "coding_sequence_variant", "inframe_insertion", "disruptive_inframe_insertion", "inframe_deletion", "disruptive_inframe_deletion", "exon_variant", "exon_loss_variant", "duplication", "inversion", "frameshift_variant", "frameshift_ablation", "gene_fusion", "bidirectional_gene_fusion", "rearranged_at_DNA_level", "miRNA", "initiator_codon_variant", "start_retained")

cand <- separate(allvar, col=mutation_effect, into=c("effect1", "effect2"), sep="&", fill="right") %>%
	filter(effect1 %in% effect | effect2 %in% effect) %>%
	filter(type == "INDEL" | (ref == "C" & alt == "T") | (ref == "G" & alt == "A")) %>%
	filter(chr == "Chr3" & pos > 18000000) %>%
	filter(ratio > 0.3)

write_tsv(cand, file.path(dirout, paste0("hcr7_candidate.tsv")))
```

- effect: snpEff annotation 중에 protein sequence나 gene expresion에 영향을 줄 수 있는 effect들.
- allvar 중에서 INDEL이거나 C>T or G>A mutation만을 고려한다.
- Causal locus안에서 ratio가 0.3 이상인 variant들만을 고려한다.


# References

[1] https://cutadapt.readthedocs.io/en/stable/index.html

[2]https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

[3] SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf

[4] https://broadinstitute.github.io/picard/explain-flags.html

[5] https://www.htslib.org/

[6] https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard
