#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
from itertools import dropwhile
import fileinput
from collections import defaultdict
import argparse
from argparse import RawTextHelpFormatter
import gzip

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description= \
	"""
	This script calculates Variant Landscape (VL) statistics across
	all coding exons. 
	
	Two input files are required:
	1. gencode annotation file (bed format)
	- I recommend using bedops gtf2bed to convert gencode's gtf file to a bed
	2. Human variation information
	- tab delimited
	- output from VEP, which gives the crucial information of the transcript \
size and the variant position in the transcript
	- currently the info is grabbed by index, which will only work with hand \
editing of this script 
	""", formatter_class=RawTextHelpFormatter)

parser.add_argument('-g','--gencode_file', required=True)
parser.add_argument('-v','--human_variation', required=True)
parser.add_argument('-o','--output_file_name', required=True)


# recommended input:
# zcat gencode.v25.annotation.bed.gz | \
# awk '$8=="CDS" {print $0}' | grep -i 'tag \"basic\"' | grep -i 'gene_type \"protein_coding\"' | grep -i 'appris_principal_1' | grep -i 'tag \"CCDS\"' | \
# ~/git/human_variation_landscape/scripts/calculate_CDS_coords.py

class FileOperations:
	def open_files(output_name):
		bed_file = open(output_name,'w')
		errors = open(
				

def build_exon_coords(gencode_file):
	tx_gene_coords = {}
	# split on the transcript entry in the 10th column and remove the quotation marks
	# build a dict with the key as the transcript and gene name and exon coords as values
	for key, chunk in groupby(gencode_file, lambda x: x.split('\t')[9].split(';')[1].split('"')[1]):
		chunk = list(chunk)
		exon_coords = [element.split('\t')[0:3] for element in chunk]
		value = []
		gene_name = chunk[0].split('\t')[9].split('"')[9]
		strand = chunk[0].split('\t')[5]
		value.append(gene_name)
		value.extend(strand)
		value.extend(exon_coords)
		dict_key = key.split('.')[0]	
		tx_gene_coords[dict_key] = value
	return(tx_gene_coords)

#errors = open('errors.txt','w')
#variant_data = open('/data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/Homo_sapiens_incl_consequences__codingOnly.VEPnoPick.GRCh38.k55.k22n.tab')
# use dropwhile iterator to skip lines staring with #
def coding_pos_processor(coding_pos):
	# convert coding positions to int and handle coordinate range
	# '67-70' situations
	int_coding = []
	for x in coding_pos:
		# if a range is given (e.g. 67-70, take the mean)
		# convert all to integers, in int_coding
		if '-' in x:
			try:
				x = int( (int(x.split('-')[0]) + int(x.split('-')[1]) ) / 2)
				int_coding.append(x)
			except:
				continue
		else:
			x = int(x)
			int_coding.append(x)
	return(int_coding)

def window_maker(tx_length, output, key, shifted_exon_coords, offsets, int_coding, chromosome, gene_name, strand):
	VL_output = open(output,'w')
	error_output = open('errors.txt', 'w')
	######
	# build 100bp windows, shifted by 1bp and calculate number of variants in the window	
	######
	for i in range(1,int(tx_length)+1):
		windowDown = -50
		windowUp = 50
		# fix if window is outside transcript
		if i < abs(windowDown):
			windowDown = -(i-1)
		if int(tx_length)-i <= 0:
			windowUp = int(tx_length)-i
		# calc number of variants in the window 
		overlapping_num = sum( [i+windowDown <= pos <= i+windowUp for pos in int_coding if pos] ) 
		# scale for interval size
		overlapping_num = str( round( overlapping_num * (100/ (windowUp - windowDown)) ) )
		# calculate actual coordinate by seeing which (shifted) exon i is in
		exon_offset_index = [x[0] for x in enumerate(shifted_exon_coords) if i in range(x[1][0],x[1][1]+1)]
		#print(gene_name, strand, tx_length, exon_offset_index, i, shifted_exon_coords)
		# now can use the index, with the offset info to calc actual genomic position
		try:
			real_genomic_position = i + offsets[int(exon_offset_index[0])]
			out = (chromosome + '\t' + str(real_genomic_position) + '\t' + \
				str(real_genomic_position + 1) + '\t' + \
				gene_name + '_' + key + '\t' + overlapping_num + '\t' + strand + '\n')
			VL_output.write(out)
		except:			
			out = gene_name + ' ' +  strand + ' ' +  str(tx_length) + ' ' + str(i) + ' ' + str(exon_offset_index) + ' ' + str(shifted_exon_coords)
			error_output.write(out)	


def calculator(variant_data, tx_gene_coords, output):
	skipped_output = open('skipped_tx.txt', 'w')

	# group on transcript name. File MUST BE SORTED BY transcript name
	# if not, then things will go TERRIBLY WRONG!!!!!!!!!!!!!!
	for key, chunk in groupby(variant_data, lambda x: x.split()[4]):
		chunk = list(chunk)
		# creates list of coding sequence positions of all variants in the gene
		coding_pos = [x.split()[8].split('/')[0] for x in chunk]
		# convert to int and handle cases where the sequence is a range (e.g 67-70)
		int_coding = coding_pos_processor(coding_pos)


		tx_length = chunk[0].split()[8].split('/')[1]
		# skip should the transcript not pass the requirements set in the initial pipe for
		# gencode (protein coding, canonical, etc.)
		if key not in tx_gene_coords:
			key = key + ' not in filtered gencode\n'
			skipped_output.write(key)
			continue
		# skip transcripts < 100 bp
		if int(tx_length) < 100:
			key = key + ' less than 100bp long\n'
			skipped_output.write(key)
			continue
		exon_coords = tx_gene_coords[key][2:]; gene_name = tx_gene_coords[key][0]
		strand = tx_gene_coords[key][1]; chromosome = tx_gene_coords[key][2][0]
		
		# remove chr from exon coords for calculations
		exon_coords = [x[1:3] for x in exon_coords]
		# convert to integer and ensure they are sorted
		exon_coords = sorted([[int(pos) for pos in x] for x in exon_coords])
		# shift stop positions by 3 to include stop codon (which CDS doesn't not include)
		if strand == '+':	
			exon_coords[-1][1] += 3
		if strand == '-':
			exon_coords[0][0] -= 3
		# prep for offset calcs
		exon_coords.insert(0,[1,1])
		# calculate offsets for each exon coordinate set
		offsets = [0]
		for i in range(1,len(exon_coords)):
			the_offset = exon_coords[i][0] - exon_coords[i-1][1]
			offsets.append(offsets[i-1]+the_offset)
		offsets.pop(0); exon_coords.pop(0)
		# calculate new exon coords from 1 to len(transcript) continuously
		# can use the offets calculated to re-make the positions
		shifted_exon_coords = [[x[1][0]-offsets[x[0]],x[1][1]-offsets[x[0]]] for x in enumerate(exon_coords)] 
		
		window_maker(tx_length, output, key, shifted_exon_coords, offsets, int_coding, chromosome, gene_name, strand)

# files (should create class)



def main():
	args = parser.parse_args()
#	if args.gencode_file[-2:]=='gz':
#		gencode_file = gzip.open(args.gencode_file,'rt', encoding='utf-8')
#	else:
	gencode_file = open(args.gencode_file)
#	if args.human_variation[-2:]=='gz':
#		human_var_file = gzip.open(args.human_variation, 'rt', encoding='utf-8')
#	else:
	human_var_file = open(args.human_variation)

	output = args.output_file_name	
	tx_gene_coords = build_exon_coords(gencode_file) 
	
	calculator(human_var_file, tx_gene_coords, output)	



# go!
main()


