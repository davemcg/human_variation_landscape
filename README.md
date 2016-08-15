# human_variation_landscape
Analysis of human variation in defined windows

Question: Does looking at the amount of known variation around a given variant inform whether the given variant is more or less likely to be functional? We would assume that since functional protein domains tend to be more conserved (citation needed), that variants that were in regions sparse of known variation are more likely to be deleterious. 

Two angles are being taken right now for the project:

1. Labeling each base pair of all known transcripts
2. Annotation of known variation

Angle 1:
Using the Gencode v25 annotation (gtf) file (see https://github.com/davemcg/human_variation_landscape/blob/master/01_initial_processing.sh for details), I convert it to a bed file, then keep the transcript coordiantes. Then I increase in each direction by 1000 base pairs, then merge overlapping transcripts together. 100bp windows, stepping by 1bp are created across the transcript coordinates. 

OK, so this new bed file is then left outer joined (loj) with release 85 of Ensembl's Homo_sapiens_incl_consequences.vcf.gz. A custom script is used (https://github.com/davemcg/human_variation_landscape/blob/master/scripts/loj_groupby_gw.py) to collapse the loj file and calculate four things: 

- number of variants in a 100bp window around each base pair
- number of variants with a MAF
- mean AND median MAF (skipping the variants without a MAF)

Four bed files are created for the four statistics. These can be used in two primary ways: first, they can be compressed into a bedGraph format (https://github.com/davemcg/human_variation_landscape/blob/master/scripts/bed_to_bedGraph.py) and next they can be used as is to merge with our VCF files to more quickly determine the variant lanscape stats for our exomes. 
