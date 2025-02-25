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
