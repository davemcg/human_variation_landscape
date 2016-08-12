#!/usr/local/bin/python3
##!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
import fileinput
import re
import sqlite3
import statistics


# print VCF header
print('##INFO=<ID=VLN100,Number=1,Type=Integer,Description="Number of known variants within 100bp of this variant">')
print('##INFO=<ID=VLN100maf,Number=1,Type=Integer,Description="Number of known variants within 100bp of this variant that have a MAF">')
print('##INFO=<ID=VLN100mean_maf,Number=A,Type=Float,Description="Mean MAF of variants within 100bp">')
print('##INFO=<ID=VLN100median_maf,Number=1,Type=Integer,Description="Median MAF of variants within 100bp">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

# rolls through file and groups by key (lambda function)
# you can then process all data in the group/chunk
for key, chunk in groupby(fileinput.input(), lambda x: x.split()[3]):
		
	chunk = list(chunk)
	num_of_var = len(chunk)-1
	chunk_bit = [s for s in chunk if key in s.split('\t')[6]][0]
	chunk_bit = chunk_bit.split('\t')
	chr = chunk_bit[0]
	pos = int( (int(chunk_bit[1]) + int(chunk_bit[2])) / 2)
	id = chunk_bit[3]
	ref = chunk_bit[7]
	alt = chunk_bit[8]
	vcf_info = str(chr) + '\t' + str(pos) + '\t' + id + '\t' + ref + '\t' + alt + '\t.\t.'

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

