#!/bin/sh
module load bedtools
module load bedops

# get Ensembl's GRCh38 v85 (latest as of 2016-08-17)  human variation file
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/GRCh38
wget http://ftp.ensembl.org/pub/release-85/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz
wget http://ftp.ensembl.org/pub/release-85/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz.tbi


# keep only coding variation
# http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
cat <(zcat /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz | head -n 1000 | grep ^#) <(zcat /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz | grep 'frameshift_variant\|inframe\|missense_variant\|start_lost\|stop_gained\|stop_lost\|synonymous_variant') | bgzip > /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences__codingOnly.vcf
tabix -p vcf /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences__codingOnly.vcf
# annotate with VEP/84, adding MAF from 1000G, ESP, ExAC as well as gene coding positions
sbatch --cpus-per-task 16 ~/git/human_variation_landscape/scripts/run_VEP.sh /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences__codingOnly.tab GRCh38 16
# reorder by transcript, then position
cat <(cat Homo_sapiens_incl_consequences__codingOnly.VEPnoPick.GRCh38.tab | head -n 1000 | grep ^#) <(grep -v ^# Homo_sapiens_incl_consequences__codingOnly.VEPnoPick.GRCh38.tab | sort -k5,5 -k2,2n) > Homo_sapiens_incl_consequences__codingOnly.VEPnoPick.GRCh38.k55.k22n.tab

# get gencode's v25 (latest as of 2016-08-17) gene annotation info
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
# convert to bed
gunzip gencode.v25.annotation.gtf.gz
gtf2bed --do-not-sort < gencode.v25.annotation.gtf > gencode.v25.annotation.bed
gzip gencode.v25.annotation.gtf
bgzip gencode.v25.annotation.bed

# keep only coding exons (CDS) from canonical genes
zcat gencode.v25.annotation.bed.gz | awk '$8=="CDS" {print $0}' | grep -i 'tag \"basic\"' | grep -i 'gene_type \"protein_coding\"' | grep -i 'appris_principal_1' | grep -i 'tag \"CCDS\"' > gencode.v25.annotation.CDS.protein-coding.principal.CCDS.bed.gz



#######################
# Computation of window counts with my own script.
# as of 2016-08-30, the script outputs allele counts in 100bp windows for all variants,
# and variants in VEP/snpeff HIGH, MODERATE, LOW categories
# b8655ec (git commit) 
sbatch --time=3:00:00 ~/git/human_variation_landscape/scripts/calculate_VL_stats.py -g gencode.v25.annotation.CDS.protein-coding.principal.CCDS.bed -v Homo_sapiens_incl_consequences__codingOnly.VEPnoPick.GRCh38.k55.k22n.h.tab -o Gencode_v25_Ensembl_v85_VEP_v84
########################

# sort, bgzip, tabix
for i in Gen*VEP*bed; do echo 'sort -k1,1 -k2,2n ' $i '> '${i%.bed}.s.bed'; bgzip ' ${i%.bed}.s.bed'; tabix -p bed ' ${i%.bed}.s.bed.gz ; done > sort_bgzip_tabix.swarm
swarm -f sort_bgzip_tabix.swarm --partition=quick












# the run_VEP script
sed -e 's/^/# /g' ~/git/human_variation_landscape/scripts/run_VEP.sh 
# #!/bin/bash
# 
# module load VEP/84
# 
# 
# input_vcf=$1
# genome=$2
# cores=$3
# 
# if [ "$genome" == "GRCh38" ] || [ "$genome" == "GRCh37" ]; then
# 	variant_effect_predictor.pl -i $input_vcf --offline \
# 	--cache --dir_cache $VEPCACHEDIR \
# 	--fasta $VEPCACHEDIR/$genome.fa --species human --assembly $genome  \
# 	--output ${input_vcf%.vcf}.VEPnoPick.$genome.tab \
#     --canonical \
# 	--coding_only \
# 		--sift s \
#     --polyphen s \
#     --symbol \
# 	--gmaf \
# 	--maf_1kg \
# 	--maf_esp \
# 	--maf_exac \
# 	--domains \
# 	--total_length \
# 	--tab --force_overwrite --fork $cores
# else
#     echo "Pick either GRCh38 or GRCh37 genomes"
# fi







