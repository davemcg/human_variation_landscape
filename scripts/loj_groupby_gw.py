#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
import fileinput
import re
import statistics
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description= \
	"""
	This script is designed to process the loj (bedtools intersect with -loj and -sorted)
	for the windows made against the genocode transcripts. It will output four bed files
	that differ only in the 5th column, the score column. The following four scores
	are calculated for EVERY position in EVERY transcript:
	1. Number of variants in the 100bp window
	2. Number of variants in the 100bp window with a MAF
	3. Mean MAF
	4. Median MAF
	""", formatter_class=RawTextHelpFormatter)

parser.add_argument('-n','--name', required=True, \
        help= \
		"""
		Give filename for the four files to be outputted. VLN100*.bed.gz will be appended. 
		* Not literally used. * will be change for each of the four outputted files.
		""")
parser.add_argument('-l','--loj_file', required=True, type=argparse.FileType('r'), \
        help = \
        'File to be collapsed and stats calculated')

args = parser.parse_args()

name = args.name

# Create files to write to
f1 = open(name + '.VLN100.bed','wb')
f2 = open(name + '.VLN100maf.bed','wb')
f3 = open(name + '.VLN100mean_maf.bed','wb')
f4 = open(name + '.VLN100median_maf.bed','wb')
# rolls through file and groups by key (lambda function)
# you can then process all data in the group/chunk

file = args.loj_file

for key, chunk in groupby(file, lambda x: x.split()[3]):
		
	chunk = list(chunk)
	num_of_var = len(chunk)
	chunk_bit = chunk[0].split('\t')
	if chunk_bit[4]=='.':
		num_of_var = 0 
	chr = chunk_bit[0]
	pos = int( (int(chunk_bit[1]) + int(chunk_bit[2])) / 2)
	id = chunk_bit[3]

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
	#print(chr,pos-1,pos,id, num_of_var, number_of_maf, mean_maf, median_maf)
	f1_data = chr + '\t' + str(pos-1) + '\t' + str(pos) + '\t' + id + '\t' + str(num_of_var) + '\n'
	f2_data = chr + '\t' + str(pos-1) + '\t' + str(pos) + '\t' + id + '\t' + str(number_of_maf) + '\n'
	f3_data = chr + '\t' + str(pos-1) + '\t' + str(pos) + '\t' + id + '\t' + str(mean_maf) + '\n'
	f4_data = chr + '\t' + str(pos-1) + '\t' + str(pos) + '\t' + id + '\t' + str(median_maf) + '\n'

	f1.write(bytes(f1_data, 'UTF-8'))
	f2.write(bytes(f2_data, 'UTF-8'))
	f3.write(bytes(f3_data, 'UTF-8'))
	f4.write(bytes(f4_data, 'UTF-8'))
