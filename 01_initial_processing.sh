#!/bin/sh
module load bedtools
module load bedops

#####
# two paths: 1 and 2
# 1 is to calculate amount of human variation in all transcripts/genes. 100bp windows. Sliding by 1bp.
# 2 is to get the human variation in 100bp windows around each known variant 

# data storage on biowulf2
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/

#####
#1 create sliding 100bp windows for all genes
#cd /home/mcgaugheyd/git/human_variation_landscape
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
gunzip gencode.v25.annotation.gtf.gz
gtf2bed < gencode.v25.annotation.gtf > gencode.v25.annotation.bed
gzip gencode.v25.annotation.gtf
gzip gencode.v25.annotation.bed
zcat gencode.v25.annotation.bed.gz| awk '$8 == "transcript" {print $0}' | gzip -f > gencode.v25.annotation.transcriptsOnly.bed.gz #only keep transcripts
gzcat 100bp_gene_windows_1bp_slide.20160701_ensembl_homo_sapiens_variation.loj.dat.gz | awk -v OFS='\t' '{key=$1"_"$2"_"$3; print $1, $2, $3, key, $4, $5, $6, $7, $8, $9, $10, $11}' | gzip -f > temp.gzi #should have added a key (bedtools has that option). Adding my own.
mv temp.gz 100bp_gene_windows_1bp_slide.20160701_ensembl_homo_sapiens_variation.loj.dat.gz
# increase each transcript size by 1000 in each direction then merge overlapping transcripts
bedtools slop -g /data/mcgaugheyd/genomes/GRCh38/hg38.chrom.sizes -b 1000 -i gencode.v25.annotation.transcriptsOnly.bed.gz | bedtools  merge -i - | gzip -f > gencode.v25.annotation.transcriptsOnly.slop_and_merged.bed.gz
# now create 100bp sliding windows, slide by 1bp
bedtools makewindows -b gencode.v25.annotation.transcriptsOnly.slop_and_merged.bed.gz -w 100 -s 1 | gzip -f > gencode.v25.annotation.transcriptsOnly.slop_and_merged.windowed.bed.gz

####
# 2 create a 100bp window around each variant  
# get ensembl variation data
# I'm using 2016-07-01 release
wget http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz
wget http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz.tbi
zcat /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz | grep -v ^# | awk -v OFS='\t' '{print $1, $2-50, $2+50, $3}' | gzip -f > 20160701_ensembl_homo_sapiens_variation.100bp_window.bed.gz &

######
#1 loj the ensemb variation data back onto the sliding 100bp windows 
bedtools intersect -b /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz -a gencode.v25.annotation.transcriptsOnly.slop_and_merged.windowed.ensembl.bed.gz -loj -sorted | gzip -f > 100bp_gene_windows_1bp_slide.20160701_ensembl_homo_sapiens_variation.loj.dat.gz

#2 loj the ensembl variation data back onto the 100bp window variant coordinates. One line for each overlap.
bedtools intersect -b /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz -a 20160701_ensembl_homo_sapiens_variation.100bp_window.bed.gz -loj -sorted | gzip -f > 20160701_ensembl_homo_sapiens_variation.100bp_window.loj.dat.gz &


