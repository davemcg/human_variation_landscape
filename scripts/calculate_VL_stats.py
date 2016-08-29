#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
from itertools import dropwhile
from itertools import islice
import fileinput
from collections import defaultdict
import argparse
from argparse import RawTextHelpFormatter
import gzip
import sys

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
	def open_files(self, output_name):
		self.VLall = open(output_name+'.VLall.bed','w')
		self.VLmoderate = open(output_name+'.VLmoderate.bed','w')
		self.VLsynonymous = open(output_name+'.VLsynonymous.bed','w')
		self.error_output = open('errors.txt', 'w')
		self.skipped_output = open('skipped_tx.txt', 'w')
				
	def writeErrors(self, info):
		self.error_output.write(info)
	def writeSkipped(self, info):
		self.skipped_output.write(info)
	def writeVLall(self, info):
		self.VLall.write(info)
		self.VLall.write('\n')
	def writeVLmoderate(self, info):
		self.VLmoderate.write(info)
		self.VLmoderate.write('\n')
	def writeVLsynonymous(self, info):
		self.VLsynonymous.write(info)
		self.VLsynonymous.write('\n')
	def close_files(self):
		self.VLall.close()
		self.VLmoderate.close()
		self.VLsynonymous.close()
		self.error_output.close()
		self.skipped_output.close()

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

def coding_pos_processor(coding_pos):
	# convert coding positions to int and handle coordinate range
	# '67-70' situations
	int_coding = []
	index = 0
	bad_pos = []
	for x in coding_pos:
		# if a range is given (e.g. 67-70, take the mean)
		# convert all to integers, in int_coding
		if '-' in x:
			try:
				x = int( (int(x.split('-')[0]) + int(x.split('-')[1]) ) / 2)
				int_coding.append(x)
			except:
				# need to ID the bad position so this can be fully skipped
				bad_pos.append(index)
		else:
			x = int(x)
			int_coding.append(x)
		index = index + 1
	return(int_coding, bad_pos)

def exon_processor(exon_coords, strand):
	# takes in list of list of exon coordinates [['chr1','100','200'],['chr1','300','600'],...]
	# and processes to remove chr and make numeric
	# and also adjusts end position by 3bp 
	# finally, it calculates the offsets and and shifts the exon coords
	# so the above coords would now be [[100,200],[200,500],...]

	# remove the chromosome from exon coords for calculations
	exon_coords = [x[1:3] for x in exon_coords]
	# convert to integer and ensure they are sorted
	exon_coords = sorted([[int(pos) for pos in x] for x in exon_coords])
	# shift stop positions by 3 to include stop codon (which CDS doesn't not include)
	if strand == '+':	
		exon_coords[-1][1] += 3
	if strand == '-':
		exon_coords[0][0] -= 3
	
	# prep for offset calcs by inserting dummy values at beginning
	exon_coords.insert(0,[1,1])
	# calculate offsets for each exon coordinate set
	offsets = [0]
	for i in range(1,len(exon_coords)):
		the_offset = exon_coords[i][0] - exon_coords[i-1][1]
		offsets.append(offsets[i-1]+the_offset)
	offsets.pop(0); exon_coords.pop(0) # remove dummy values
	
	# calculate new exon coords from 1 to len(transcript) continuously
	# will later use the offsets calculated to re-calce the actual genomic positions
	shifted_exon_coords = [[x[1][0]-offsets[x[0]],x[1][1]-offsets[x[0]]] for x in enumerate(exon_coords)]
	return(shifted_exon_coords, offsets)

def window_calc_print(tx_length, \
						fileHandler, \
						key, \
						exon_coords, \
						coding_pos, \
						chromosome, \
						gene_name, \
						strand, \
						index_to_keep):
	# build 100bp windows, shifted by 1bp and calculate number of variants in the window	
	# prints data to file
	
	# remove chr, make numeric, adjust stop coordinate, calculate offsets, shifts the exon coordinates
	shifted_exon_coords, offsets = exon_processor(exon_coords, strand) 
	
	output = []

	# create new position list, removing nonrelevant entries
	try:
		coding_pos = [x[1] for x in enumerate(coding_pos) if x[0] in index_to_keep]
	except:
		print('Can\'t create new coding_pos list')
		out = gene_name + '\n' + str(index_to_keep) + '\n' + str(coding_pos) + '\n' + key
		print(out)
		sys.exit(1)
	for i in range(1,int(tx_length)+1):
		windowDown = -50
		windowUp = 50
		# fix if window is outside transcript
		if i < abs(windowDown):
			windowDown = -(i-1)
		if int(tx_length)-i <= 0:
			windowUp = int(tx_length)-i
		# only keep variants that meet a certain criteria
		# calc number of variants in the window 
		overlapping_num = sum( [i+windowDown <= pos <= i+windowUp for pos in coding_pos if pos] ) 
		# scale for interval size
		overlapping_num = str( round( overlapping_num * (100 / (windowUp - windowDown)) ) )
		# calculate actual coordinate by first seeing which (shifted) exon i is in
		exon_offset_index = [x[0] for x in enumerate(shifted_exon_coords) if i in range(x[1][0],x[1][1]+1)]
		# print(gene_name, strand, tx_length, exon_offset_index, i, shifted_exon_coords)
		# now can use the index, with the offset info to calc actual genomic position
		try:
			real_genomic_position = i + offsets[int(exon_offset_index[0])]
			out = (chromosome + '\t' + str(real_genomic_position) + '\t' + \
				str(real_genomic_position + 1) + '\t' + \
				gene_name + '_' + key + '\t' + overlapping_num + '\t' + strand)
		except:			
			out = gene_name + ' ' +  strand + ' ' +  str(tx_length) + ' ' + str(i) + ' ' + str(exon_offset_index) + ' ' + str(shifted_exon_coords)
			fileHandler.writeErrors(out)	
			continue

		output.append(out)
	return(output)

	
def central_depot(variant_data, tx_gene_coords, output):
	fileHandler = FileOperations()
	fileHandler.open_files(output)
	# grab header
	top_line = []
	for line in islice(variant_data, 200):
		top_line.append(line[:-1])
	# grab header and parse for various column names
	header = ([x for x in top_line if x[0:19]=='#Uploaded_variation'])[0]
	cds_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'CDS_position'])[0]
	gmaf_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'GMAF'])[0]
	impact_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'IMPACT'])[0]
	consequence_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'Consequence'])[0]
	
	# skip header lines
	for line in dropwhile(lambda line: line.startswith('#'), variant_data):	
		# group on transcript name. File MUST BE PRE-SORTED BY transcript name
		# if not, then things will go TERRIBLY WRONG!!!!!!!!!!!!!!
		for key, chunk in groupby(variant_data, lambda x: x.split()[4]):
			chunk = list(chunk)
			# creates list of coding sequence positions of all variants in the gene
			# then coding_pos_processor:
			# convert to int and handle cases where the sequence is a range (e.g 67-70) 
			coding_pos, bad_pos = coding_pos_processor([x.split('\t')[cds_index].split('/')[0] for x in chunk])
			gmaf = [x.split('\t')[gmaf_index] for x in chunk]
			impact = [x.split('\t')[impact_index] for x in chunk]
			consequence = [x.split('\t')[consequence_index] for x in chunk]
			# if we have poorly formatted positions, remove from our index so we don't process it
			if len(bad_pos)>0:
				gmaf = [x[1] for x in enumerate(gmaf) if x[0] not in bad_pos]
				impact = [x[1] for x in enumerate(impact) if x[0] not in bad_pos]
				consequence = [x[1] for x in enumerate(consequence) if x[0] not in bad_pos]
			moderate_impact_index = [x[0] for x in enumerate(impact) if x[1] == 'MODERATE']
			synonymous_index = [x[0] for x in enumerate(consequence) if [x] == 'synonymous_variant']
			
			tx_length = chunk[0].split()[cds_index].split('/')[1]
			# skip should the transcript not pass the requirements set in the initial pipe for
			# gencode (protein coding, canonical, etc.)
			if key not in tx_gene_coords:
				key = key + ' not in filtered gencode\n'
				fileHandler.writeSkipped(key)
				continue
			# skip transcripts < 100 bp
			if int(tx_length) < 100:
				key = key + ' less than 100bp long\n'
				fileHandler.writeSkipped(key)
				continue
	
			# grab exon coordinates and strand from the tx_gene_coords dictionary
			exon_coords = tx_gene_coords[key][2:]; gene_name = tx_gene_coords[key][0]
			strand = tx_gene_coords[key][1]; chromosome = tx_gene_coords[key][2][0]
			# remove chr, make numeric, adjust stop coordinate, calculate offsets, shifts the exon coordinates
			shifted_exon_coords, offsets = exon_processor(exon_coords, strand)	
			
			# this does all the math and prints
			# first all variants
			all_positions = list(range(0,len(coding_pos)))
			VLall=window_calc_print(tx_length, \
								fileHandler, \
								key, \
								exon_coords, \
								coding_pos, \
								chromosome, \
								gene_name, \
								strand, \
								all_positions)
			fileHandler.writeVLall('\n'.join(VLall))
			# now just print out nums of surrounding moderate alleles
			VLmoderate=window_calc_print(tx_length, \
								fileHandler, \
								key, \
								exon_coords, \
								coding_pos, \
								chromosome, \
								gene_name, \
								strand, \
								moderate_impact_index)
			fileHandler.writeVLmoderate('\n'.join(VLmoderate))
			# now nums of surrounding synonymous alleles
			VLsynonymous=window_calc_print(tx_length, \
								fileHandler, \
								key, \
								exon_coords, \
								coding_pos, \
								chromosome, \
								gene_name, \
								strand, \
								synonymous_index)
			fileHandler.writeVLsynonymous('\n'.join(VLsynonymous))
		fileHandler.close_files()


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
	
	central_depot(human_var_file, tx_gene_coords, output)	



# go!
main()


