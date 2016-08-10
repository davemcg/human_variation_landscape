#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import re
import fileinput

data = {}
# read from stdin
for line in fileinput.input():	

	line = line.split('\t')	

	key = line[3]
	window_key = line[6]
	info = line[11]

	# grab info to build vcf at the end
	chr = line[4]
	pos = line[5]
	id = line[6]
	ref = line[7]
	alt = line[8]
	# lump together
	vcf_info = str(chr) + '\t' + str(pos) + '\t' + id + '\t' + ref + '\t' + alt + '\t.\t.'

	regex=re.compile(r'.*MAF=0\.\d+', re.I) # regex pattern to find the maf

	if key == window_key: # skip cases where the variants are the same 
		continue
	maf = 0 # to cover when no maf is present
	if regex.search(info):
		maf = [float(m.group().split('=')[1]) for section in info.split(';') for m in [regex.search(section)] if m][0]
	if key not in data:
		list_maf = []
		list_maf.append(maf)
		stats = 1,list_maf, vcf_info
		data[key] = stats
	if key in data:
		num_var, the_mafs, vcf_info = data[key]
		num_var += 1
		the_mafs.append(maf)
		data[key] = num_var, the_mafs, vcf_info

print('##fileDate=20160810')
print('##INFO=<ID=SVN100,Number=1,Type=Integer,Description="Number of variants within 100bp of this variant">')
print('##INFO=<ID=SVN100maf,Number=1,Type=Integer,Description="Number of variants with a recorded MAF within 100bp of this variant">')
print('##INFO=<ID=SVN100mean_maf,Number=A,Type=Float,Description="Mean MAF, excluding variants with no MAF in the 100bp window around this variant">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
for k in data.keys():
	num_var, mafs, vcf_info = data[k]
	num_mafs = sum(i > 0 for i in mafs)
	if num_mafs == 0:
		mean_mafs = 0
	else:
		mean_mafs = sum(mafs) / num_mafs
	stats = '\tSVN100=' + str(num_var) + ';SVN100maf=' + str(num_mafs) + ';SVN100mean_maf=' + str(mean_mafs)
	print(vcf_info, stats)
