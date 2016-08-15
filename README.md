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

Four bed files are created for the four statistics. These can be used in two primary ways: first, they can be compressed into a bedGraph format (https://github.com/davemcg/human_variation_landscape/blob/master/scripts/bed_to_bedGraph.py) to view more easily on the UCSC genome browser and next they can be used as is to merge with our VCF files to more quickly determine the variant lanscape stats for our exomes. 

Angle 2:
Annotate every known variant (Ensembl's human variation version 85) with surrounding variation (again, 100bp window. 50 up and 50 down). Like Angle 1, see https://github.com/davemcg/human_variation_landscape/blob/master/01_initial_processing.sh for details. Simpler than above. Basically just use awk to take each coordinate and add 50 bp up and down. Then do the loj, then use a modification of loj_groupby_gw.py (https://github.com/davemcg/human_variation_landscape/blob/master/scripts/loj_groupby_vcf.py) which creates the four stats in Angle 1, and outputs a single new VCF with the new stats.

What is this for? To test whether this theory holds ANY water. First, we can calculate basic stats (summary stats, histograms, break by type of variant (indel, SNP). Then we can check the ClinVar Path/Not Path annotation to see whether the pathogenic annotated variants have a significant difference in the amount of known variation around them. Could also check against CADD, Phred, etc. 
