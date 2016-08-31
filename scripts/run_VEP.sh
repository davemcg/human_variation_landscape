#!/bin/bash

module load VEP/84


input_vcf=$1
genome=$2
cores=$4
file_type=$3

if [ "$genome" == "GRCh38" ] || [ "$genome" == "GRCh37" ] && [ "$file_type" == "tab" ]; then
	variant_effect_predictor.pl -i $input_vcf --offline \
	--cache --dir_cache $VEPCACHEDIR \
	--fasta $VEPCACHEDIR/$genome.fa --species human --assembly $genome  \
	--output ${input_vcf%.vcf}.VEPnoPick.$genome.tab \
    --canonical \
	--coding_only \
		--sift s \
    --polyphen s \
    --symbol \
	--gmaf \
	--maf_1kg \
	--maf_esp \
	--maf_exac \
	--plugin CADD \
	--domains \
	--total_length \
	--tab --force_overwrite --fork $cores

elif [ "$genome" == "GRCh38" ] || [ "$genome" == "GRCh37" ] && [ "$file_type" == "vcf" ]; then
	variant_effect_predictor.pl -i $input_vcf --offline \
	--cache --dir_cache $VEPCACHEDIR \
	--fasta $VEPCACHEDIR/$genome.fa --species human --assembly $genome  \
	--output ${input_vcf%.vcf}.VEPnoPick.$genome.vcf \
    --canonical \
	--coding_only \
		--sift s \
    --polyphen s \
    --symbol \
	--gmaf \
	--maf_1kg \
	--maf_esp \
	--maf_exac \
	--plugin CADD \
	--domains \
	--total_length \
	--vcf --force_overwrite --fork $cores

else
    echo "Pick either GRCh38 or GRCh37 genomes and vcf or tab file output types"
fi
