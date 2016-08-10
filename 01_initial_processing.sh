#!/bin/sh
module load bedtools

# get ensembl variation data
# I'm using 2016-07-01 release
cd /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/
wget http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz
wget http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/Homo_sapiens_incl_consequences.vcf.gz.tbi

# create a 100bp window around each variant  
zcat /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz | grep -v ^# | awk -v OFS='\t' '{print $1, $2-50, $2+50, $3}' | gzip -f > /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/20160701_ensembl_homo_sapiens_variation.100bp_window.bed.gz &

# loj the ensembl variation data back onto the variant coordinates created above. One line for each overlap.
bedtools intersect -b /data/mcgaugheyd/genomes/GRCh38/Homo_sapiens_incl_consequences.vcf.gz -a /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/20160701_ensembl_homo_sapiens_variation.100bp_window.bed.gz -loj -sorted | gzip -f > /data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/20160701_ensembl_homo_sapiens_variation.100bp_window.loj.dat.gz &
