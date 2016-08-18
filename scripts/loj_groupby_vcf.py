#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
import fileinput
import re
import sqlite3
import statistics
import operator

"""
Coded for the ensembl v85 targeted variations (approach 2). Will take the loj
file as input (through pipe) created with bedtools intersect with -loj and -sorted.

Requires the input to have CHR_POS_REF_ALT in the vcf 'ID' position

Outputs a vcf with four custom tags (see below for details) for each of the ensembl
variants to stdout

"""
# print VCF header
print('##fileformat=VCFv4.1')
print('##INFO=<ID=VLN100,Number=1,Type=Integer,Description="Number of known variants within 100bp of this variant">')
print('##INFO=<ID=VLN100maf,Number=1,Type=Integer,Description="Number of known variants within 100bp of this variant that have a MAF">')
print('##INFO=<ID=VLN100mean_maf,Number=A,Type=Float,Description="Mean MAF of variants within 100bp">')
print('##INFO=<ID=VLN100median_maf,Number=A,Type=Float,Description="Median MAF of variants within 100bp">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

# rolls through file and groups by key (lambda function)
# you can then process all data in the group/chunk
for key, chunk in groupby(fileinput.input(), lambda x: x.split()[3]):
	
	chunk = list(chunk)
	try:
		# this horrifying list comprehension does the following:
		# for each piece of the chunk (grouped by the key) we check if the key,
		# which is CHR_POS_REF_ALT matches the derived key (the itemgetter part)
		# which makes CHR_POS_REF_ALT for the loj file
		# then the rs1242134whatever is grabbed
		rs = [s.split('\t')[6] for s in chunk if key == '_'.join(operator.itemgetter(4,5,7,8)(s.split('\t'))) ][0]
		num_of_var = len(chunk)-1
		# new chunk without overlapping variant
		chunk = [s for s in chunk if key != '_'.join(operator.itemgetter(4,5,7,8)(s.split('\t'))) ]
	except:
		num_of_var = len(chunk)
		rs = '.'
	vcf_info = key.split('_')
	vcf_info.insert(2,rs)
	vcf_info.append('.\t.')
	vcf_info = '\t'.join(vcf_info)

	chunk = ';'.join(chunk)
	regex=re.compile(r'MAF=0\.\d+', re.I) # regex pattern to find the maf
	maf = [float(m.group().split('=')[1]) for section in chunk.split(';') for m in [regex.search(section)] if m]
	number_of_maf = len(maf)
	if len(maf) == 0:
		mean_maf = 0
		median_maf = 0 
	else:
		mean_maf = sum(maf) / len(maf)		
		median_maf = statistics.median(maf)
	info = '\tVLN100=' + str(num_of_var) + ';VLN100maf=' + \
			str(number_of_maf) + ';VLN100mean_maf=' + \
			str(mean_maf) + ';VLN100median_maf=' + str(median_maf)
	print(vcf_info, info)

