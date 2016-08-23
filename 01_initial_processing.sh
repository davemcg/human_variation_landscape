#!/bin/sh
module load bedtools
module load bedops

#####
# two (three) paths: 1 (1_sub) and 2
# 1 is to calculate amount of human variation in all transcripts/genes. 100bp windows. Sliding by 1bp.
#	- output are four bed files with different stats
#	- GRCh38 Ensembl v85 and gencode v25 (latest versions as of 2016-08-17)
# 	- 1_sub is the same, but done with GRCh37 Ensembl v83
#		and gencode v25 lifted onto GRCh37 (latest versions as of 2016-08-17)
# 2 is to get the human variation in 100bp windows around each known variant 
#	- output is a new VCF file with the four stats in the VCF

# data storage on biowulf2
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/

#####
#1 alter gencode gtf to give bed style coordinates for each transcript
# mostly uses bedops with a awk step to only keep transcripts
# download it. v25 latest as of 2016-08-17
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
gunzip gencode.v25.annotation.gtf.gz
gtf2bed --do-not-sort < gencode.v25.annotation.gtf > gencode.v25.annotation.bed
gzip gencode.v25.annotation.gtf
gzip gencode.v25.annotation.bed
zcat gencode.v25.annotation.bed.gz| awk '$8 == "transcript" {print $0}' | gzip -f > gencode.v25.annotation.transcriptsOnly.bed.gz #only keep transcripts

# get Ensembl's GRCh38 v85 (latest as of 2016-08-17)  human variation file
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/GRCh38
wget http://ftp.ensembl.org/pub/release-85/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz
wget http://ftp.ensembl.org/pub/release-85/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz.tbi

# move back to project location on biowulf2
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape

# set of piped commands to process the annotation info and create four bed files with 100bp windows (1bp steps) for stats
# on the number of variants in the 100bp vicinity
    # increase each transcript size by 1000 in each direction then merge overlapping transcripts
    bedtools slop -g /data/mcgaugheyd/genomes/GRCh38/hg38.chrom.sizes -b 1000 -i gencode.v25.annotation.transcriptsOnly.bed.gz | \
	# get chr pos sort 
	sort -k1,1 -k2,2n | \
    # merge overlaps
    bedtools merge -i - | \
    # convert gencode to ensembl chr notation (chr1 to 1)
    ~/git/ChromosomeMappings/./convert_notation.py -c ~/git/ChromosomeMappings/GRCh38_gencode2ensembl.txt -f - | \
	# get chr sort correct for new naming
	sort -k1,1 -k2,2n | \
	 # create 100bp windows, sliding by one
    bedtools makewindows -b - -w 100 -s 1  | \
    # create ID
    awk -v OFS='\t' '{key=$1"_"$2"_"$3; print $1, $2, $3, key}' | \
    # loj back onto the ensembl variation file
    bedtools intersect -b /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz -a - -loj -sorted -g /data/mcgaugheyd/genomes/GRCh38/GRCh38.ensembl.k11sort.chrom.sizes | \
    ~/git/human_variation_landscape/scripts/loj_groupby_gw.py -n Gencode_v25_Ensembl_v85 -l -


# can compress bed files down to bedGraph with scripts/bed_to_bedGraph.py
swarm -f ~/git/human_variation_landscape/scripts/run_bed_to_bedGraph.swarm

#########
#1_sub do the same as above, but with GRCh37 annotations/build
# get Ensembl's GRCh37 human variation file
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/GRCh37
wget http://ftp.ensembl.org/pub/grch37/release-83/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz
wget http://ftp.ensembl.org/pub/grch37/release-83/variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz.tbi

# move to data location on biowulf2
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/GRCh37

# grab GRCh37 gencode genes
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz
gunzip gencode.v25lift37.annotation.gtf.gz
gtf2bed --do-not-sort < gencode.v25lift37.annotation.gtf > gencode.v25lift37.annotation.bed #convert gtf to bed
# compress again
gzip gencode.v25lift37.annotation.gtf
gzip gencode.v25lift37.annotation.bed
# only keep transcripts
zcat gencode.v25lift37.annotation.bed.gz| awk '$8 == "transcript" {print $0}' | gzip -f > gencode.v25lift37.annotation.transcriptsOnly.bed.gz


# set of piped commands to process the annotation info and create four bed files with 100bp windows (1bp steps) for stats
# on the number of variants in the 100bp vicinity
    # increase each transcript size by 1000 in each direction then merge overlapping transcripts
    bedtools slop -g /data/mcgaugheyd/genomes/GRCh37/GRCh37.gencode.chrom.sizes -b 1000 -i gencode.v25lift37.annotation.transcriptsOnly.bed.gz | \
	# get chr pos sorted
	sort -k1,1 -k2,2n | \
    # merge overlaps
    bedtools merge -i - | \
    # convert gencode to ensembl chr notation (chr1 to 1)
    ~/git/ChromosomeMappings/./convert_notation.py -c ~/git/ChromosomeMappings/GRCh37_gencode2ensembl.txt -f - | \
 	# ensure sort is correct for new chr naming
	sort -k1,1 -k2,2n | \    
	# create 100bp windows, sliding by one
    bedtools makewindows -b - -w 100 -s 1  | \
    # create ID
    awk -v OFS='\t' '{key=$1"_"$2"_"$3; print $1, $2, $3, key}' | \
    # loj back onto the ensembl variation file
    bedtools intersect -b /data/mcgaugheyd/genomes/GRCh37/Homo_sapiens_incl_consequences.vcf.gz -a - -loj -sorted -g /data/mcgaugheyd/genomes/GRCh37/GRCh37.ensembl.k11_sort.chrom.sizes |
    ~/git/human_variation_landscape/scripts/loj_groupby_gw.py -n Gencode_v25lift37_Ensembl_v83 -l -





####
# 2. Annotate existing variation for surrouding variation
# skip header lines
zcat /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz | grep -v ^# | \
# create a 100bp window around each variant and make a CHR_POS_REF_ALT key
awk -v OFS='\t' '{key=$1"_"$2"_"$4"_"$5; print $1, $2-50, $2+50, key}' | \
# intersect against the original variation file and do a loj
bedtools intersect -b /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz -a - -loj -sorted | \
# send to my script which collapses the overlapping variants and calculates the four stats and writes
# a new VCF with the stats
~/git/human_variation_landscape/scripts/loj_groupby_vcf.py > 20160701_ensembl_homo_sapiens_variation.100bp_window.vcf

# compress
bgzip 20160701_ensembl_homo_sapiens_variation.100bp_window.vcf
tabix -p vcf 20160701_ensembl_homo_sapiens_variation.100bp_window.vcf.gz




