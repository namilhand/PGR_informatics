# Filtering SNPs for GBS

**Objective: To retain highly reliable SNPs for GBS**

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
		fi
		    
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

## 3. joint genotyping

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

## 4. Filtering variants

src:

Reliable한 variant를 골라내는 과정. 분석에 사용되는 F2 개체수와 sequencing depth등에 따라 cutoff를 조정해야 한다.

### 1) Filter genotyping result by GQ

Genotype 결과 중에 GQ score가 충분히 높지 않은 것들을 제거한다.

```r
#!/bin/bash

export wd="/datasets/data_4/nison/PGR_informatics/Lesson02_FilteringSNPs/result"
export dirvcf="$wd/03_joint_genotyping"
export dir_result="$wd/04_filtering_variants"
mkdir -p $dir_result

cd $wd

# 1. Filter genotyping result by GQ
# Passing variants are annotated as PASS
# Failing variants are annotated with the name of the filter they failed
for f in $dirvcf/*.vcf.gz; do
    vcf=$f
	input_name=$(basename $vcf)
    filtered_vcf=${input_name%.vcf.gz}.VF.vcf
    gatk VariantFiltration \
        -V $vcf \
        -O $dir_result/$filtered_vcf \
        --genotype-filter-expression "GQ < 10" \
        --genotype-filter-name "GQ10" &
done

wait;
echo "finished filtering vcf"
```

### 02) Filter by allele frequency

allele frequency가 0.1보다 낮은 variant들을 제거한다.

```r
# 2. Filter by allele frequency
# Filter SNP positions by allele frequency
# If allele frequency is lower than 0.1, mark those positions

now=$(date +"%T")
echo "Current time: $now"
for f in $dir_result/*.VF.vcf; do
    input=$f
    output=$(basename $input)
    output=${output%.vcf}.AFfilter.vcf.gz

    bcftools filter $input --exclude 'INFO/AF[0] < 0.1' -O z -o $dir_result/$output &
    echo filtering $input by AF
done
wait;
echo finished AF filtering
```

### 03) SelectVariants

Variant중 SNP만 골라내고, GQ에 의해 filtering된 sample들의 genotype을 nocall (./.)로 바꾼다. 이 과정을 거치고 나면 GQ filter를 통과한 sample들만 mendelian segregation test에 사용된다.

```r
# 3. Select variants
# Select a subset of variants from a vcf file

ref="/datasets/data_3/genome_assemblies/TAIR10_Chr15.fasta"

echo "Current time: $now"
cd $dir_result
for f in *.VF.AFfilter.vcf.gz; do
    gatk IndexFeatureFile \
        -I $f &
done
wait;
echo "Current time: $now"
echo finished indexing
# gatk SelectVariants require vcf file to be indexed.

echo "Current time: $now"
for f in *.VF.AFfilter.vcf.gz; do
    gatk SelectVariants \
        -R $ref \
        -V $f \
        --select-type-to-include SNP \
        --set-filtered-gt-to-nocall TRUE \
        -O ${f%.vcf.gz}.SelVar.vcf &
    done
wait;
echo "Current time: $now"

```

### 04) variants to table

vcf 형식의 파일을 원하는 format으로 reformatting한다. 이후에 데이터 처리를 용이하게 하기 위해서 진행한다.

```r
mkdir -p $dir_result/variant_tables

function vartotab {
    input=$1
    output=$(basename $input)
    output_prefix=${output%.vcf}
    

    gatk VariantsToTable \
        -V $input \
        -F CHROM -F POS -F REF -F ALT -F TYPE -F DP -F AF -F NSAMPLES -F HOM-REF -F HET -F HOM-VAR -F NCALLED -F NO-CALL -F VAR \
        --split-multi-allelic \
        -O $dir_result/variant_tables/${output_prefix}.temp

	# -F NSAMPLES: number of samples
	# -F HOM-REF: count of homozygous reference genotypes
	# -F HET: count of heterozygous genotypes
	# -F HOM-VAR: count of homozygous variant genotypes
	# -F NCALLED: number of called samples
	# -F NO-CALL: count of no-call genotypes
	# -F VAR: count of non-reference genotypes

    awk 'BEGIN{FS="\t"; OFS="\t"}NR==1{print; next} $7 > 0.1 {print}' $dir_result/variant_tables/${output_prefix}.temp > $dir_result/variant_tables/${output_prefix}.tsv
	# discard SNPs with allele frequency lower than or equal to 0.1

    rm $dir_result/variant_tables/${output_prefix}.temp
    }

export -f vartotab

parallel --jobs 32 vartotab {} ::: $dir_result/*.SelVar.vcf
```

## 5. Mendelian segregation test

각 variant마다 fisher’s exact test를 이용해 hom-ref/het/hom-alt 비율이 1/2/1과 같은지 test한다.

아래와 같은 기준으로 최종 SNP list를 결정한다. 기준은 상황에 따라 수정할 수 있다.

1. The number of samples with genotype called (NCALLED) > 20
2. allele frequency 0.4-0.6
3. exact test p-value > 0.1

exact test를 모든 variant에 적용하는데 시간이 걸리기 때문에 `foreach`와  `doParallel`을 사용해서 여러 CPU를 동시에 사용한다.

(https://blog.naver.com/wjddudwo209/220616435793)

```r
library(tidyverse)
library(foreach)
library(doParallel)

#====== for parallel processing ======
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

#====== input ======
dirtsv <- "PGR_informatics/Lesson02_FilteringSNPs/result/04_filtering_variants/variant_tables"

tsv_files <- list.files(dirtsv, pattern=".tsv")
tsv_list <- lapply(file.path(dirtsv, tsv_files), read_tsv)
## total 334795 variants

markers <- bind_rows(tsv_list) %>% arrange(CHROM, POS) %>%
    filter(NCALLED != 0)
## variants with NCALLED = 0 result in errors during fisher's exact test

colnames(markers) <- c("chr", "pos", "ref", "alt", "type", "dp", "af", "nsamples", "homref", "het", "homvar", "ncalled", "nocall", "var")

# NOTES: how are these markers obtained
# 1. allele frequency higher than 0.1
# 2. only snp

#====== output ======
dirout="result/05_final_snp"
dir.create(dirout, recursive=T)

#====== mendelian segregation test =======

mendelTest <- function(dat){
    # dat <- markers[1,]
    observed <- c(RR=dat$homref, RA=dat$het, AA=dat$homvar)
    expected <- c(RR=1/4, RA=1/2, AA=1/4) * sum(observed)
    chisq.table <- rbind(observed, expected)
    p <- fisher.test(chisq.table)$p.value

    # why fisher's exact test not chi-square test?
    return(p)
}

mendel_test_result <- foreach(i=1:nrow(markers), .combine=c) %dopar%{
    mendelTest(markers[i,])
}
markers$pval <- mendel_test_result

#======= exploratory analysis on SNPs ======
# 334,795 variants

my_theme <- theme_classic() +
    theme(text = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.text = element_text(size=6, family="Arial", colour="black")) +
    theme(plot.title = element_text(size=6, family="Arial", colour="black")) +
    theme(axis.title = element_text(size=6, family="Arial", colour="black")) +
    theme(legend.margin = margin(0.5,0.5,0.5,0.5)) +
    theme(strip.text = element_text(size=6, family="Arial", colour="black")) +
    theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme(axis.line = element_line(linewidth=0.2), axis.ticks=element_line(linewidth=0.2))
theme_set(my_theme)

# 1. NCALLED

mean.ncalled <- mean(markers$ncalled) # 39.05

p.hist_ncalled <- ggplot() +
    geom_histogram(data=markers, aes(x=ncalled), binwidth=1) +
    geom_vline(xintercept=20, colour="red", linewidth=0.4)

png(file=file.path(dirout, "hist_ncalled.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.hist_ncalled)
dev.off()

## ncalled > 10 : 332601
## ncalled > 15 : 326175
## ncalled > 20 : 311,608

# 2. allele frequency

mean.af <- mean(markers$af) # 0.47

p.hist_af <- ggplot() +
    geom_histogram(data=markers, aes(x=af), binwidth=0.01) +
    geom_vline(xintercept=0.4, colour="red", linewidth=0.4) +
    geom_vline(xintercept=0.6, colour="red", linewidth=0.4)

png(file=file.path(dirout, "hist_allele_frequency.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.hist_af)
dev.off()

## AF 0.4 - 0.6 : 235,004 variants

# 3. fisher's exact test p-value

p.pval <- ggplot() +
    geom_histogram(data=markers, aes(x=pval), binwidth=0.01) +
    geom_vline(xintercept=0.1, colour="red", linewidth=0.4)

png(file=file.path(dirout, "hist_exact-test.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.pval)
dev.off()

## variants with pvalue < 0.05 : 34,036
## variants with pvalue < 0.10 : 52,228

# 4. pvalue by allele frequency

p.pval_af <- ggplot() +
    geom_point(data=markers, aes(x=af, y=pval), shape=21, colour="NA", fill="black", alpha=0.3, size=0.1) +
    geom_vline(xintercept=c(0.4, 0.6), colour="red", linewidth=0.4) +
    geom_hline(yintercept=0.1, colour="red", linewidth=0.4) +
    scale_y_continuous(limits=c(0, 1)) +
    scale_x_continuous(limits=c(0, 1))

png(file=file.path(dirout, "plot_pvalue-by-af.png"), width=2.5, height=2.5, unit="in", res=300)
print(p.pval_af)
dev.off()

#====== Final SNP sets ======
# 1. NCALLED > 20
# 2. AF 0.4 - 0.6
# 3. chi-square p-value > 0.1

markers.fin <- filter(markers, ncalled > 20, af > 0.4 & af < 0.6, pval > 0.1)
# 223,870 snps

markers.fin.tsv <- dplyr::select(markers.fin, c(chr, pos, ref, alt))

write_tsv(markers.fin.tsv, file="result/05_final_snp/coller_gbs_markers.tsv", col_names=F)
```

**Result**

- Total 223,870 Col/Ler Syntenic SNPs were retained

## References

[1] minimap2: https://lh3.github.io/minimap2/minimap2.html

[2] syri: [https://schneebergerlab.github.io/syri/](https://schneebergerlab.github.io/syri/)

[3] plotsr: https://schneebergerlab.github.io/syri/plotsr.html