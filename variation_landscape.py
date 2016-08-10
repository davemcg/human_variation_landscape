#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
from argparse import RawTextHelpFormatter
import re
import subprocess
from subprocess import Popen, PIPE
import os
import statistics as stat

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description= \
	"""
	Uses tabix to grab overlapping fields (bed format, user defined) against a user-given \
vcf file. The output is then parsed to give:\n \
	1. The number of variants overlapping the window\n \
	2. The number of variants with a recorded MAF\n \
	3. The average MAF of the above variants (highest MAF used)
	""", formatter_class=RawTextHelpFormatter)
parser.add_argument('-v','--vcf', required=True, \
	help = \
	'Give vcf file from ensembl (http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/, \
\"Homo_sapiens_incl_consequences.vcf.gz\") to search for variation across positions')
parser.add_argument('-b','--bed', required=True, \
	help =\
	'Coordinates to search for variation in bed format. Must be formatted the same as vcf')
parser.add_argument('-e','--expand', default=0, type=int, \
	help = \
	'Expansion (in bp) up and down from position. Default is 0. If, say, 50 is given, then 100 bp \
total will be added')

def call_tabix(vcf, bed, expand):
	# example tabix
	# tabix -p vcf the_file.vcf.gz 1:10000-11000
	bed_data = open(bed,'r')
	stats = []
	for line in bed_data:
		sline = line.split()
		region = sline[0] + ':' + str(int(sline[1])-expand) + '-' + str(int(sline[2])+expand)
		tabix_query = 'tabix -p vcf ' + vcf + ' ' + region
		tabix = subprocess.Popen(tabix_query, shell=True, stdout=subprocess.PIPE, bufsize=0)
		tabix = str(tabix.communicate()[0])
	#	tabix = tabix.split('\n')

#		tabix = os.popen(tabix_query)
#		tabix = tabix.readline()

		stats.append(analyze_tabix(tabix))
	return(stats)

def analyze_tabix(tabix_data):
	# rolls through tabix return, which are lines for each variant in the region given
	# calculates three values:
	# 	1. number of variants in region
	#	2. number of variants with a MAF
	#	3. median MAF (just from variants with a MAF)
	num_of_variant = len(tabix_data.split('\n')) 	# number 1
	regex=re.compile(r'.*MAF=0\.\d+', re.I)	# regex pattern to find the maf 
	mafs = [float(m.group().split('=')[1]) for section in tabix_data[:-1].split(';') for m in [regex.search(section)] if m]
	if len(mafs) == 0:
		median = 0
	if len(mafs) > 0:
		median = stat.median(mafs)
	return(num_of_variant, len(mafs), median)


def main():
	args = parser.parse_args()
	vcf = args.vcf
	bed = args.bed
	expand = args.expand
	stats = call_tabix(vcf, bed, expand)
	for data in stats:
		print(data)

# let's go!
main()

