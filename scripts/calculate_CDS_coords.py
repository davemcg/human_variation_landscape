#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
from itertools import dropwhile
import fileinput
from collections import defaultdict
# recommended input:
# zcat gencode.v25.annotation.bed.gz | \
# awk '$8=="CDS" {print $0}' | grep -i 'tag \"basic\"' | grep -i 'gene_type \"protein_coding\"' | grep -i 'appris_principal_1' | grep -i 'tag \"CCDS\"' | \
# ~/git/human_variation_landscape/scripts/calculate_CDS_coords.py

tx_gene_coords = {}
# split on the transcript entry in the 10th column and remove the quotation marks
# build a dict with the key as the transcript and gene name and exon coords as values
for key, chunk in groupby(fileinput.input(), lambda x: x.split('\t')[9].split(';')[1].split('"')[1]):
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


variant_data = open('/data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/test.tab')
# use dropwhile iterator to skip lines staring with #
for line in dropwhile(lambda line: line.startswith('#'), variant_data):
	# group on transcript name. File MUST BE SORTED BY transcript name
	# if not, then things will go TERRIBLY WRONG!!!!!!!!!!!!!!
	for key, chunk in groupby(variant_data, lambda x: x.split()[4]):
		chunk = list(chunk)
		coding_pos = [x.split()[8].split('/')[0] for x in chunk]
		int_coding = []
		# if a range is given (e.g. 67-70, take the mean)
		# convert all to integers, in int_coding
		for x in coding_pos:
			if '-' in x:
				try:
					x = int( (int(x.split('-')[0]) + int(x.split('-')[1]) ) / 2)
					int_coding.append(x)
				except:
					continue
			else:
				x = int(x)
				int_coding.append(x)
		tx_length = chunk[0].split()[8].split('/')[1]
		# skip should the transcript not pass the requirements set in the initial pipe for
		# gencode (protein coding, canonical, etc.)
		if key not in tx_gene_coords:
			continue
		# skip transcripts < 100 bp
		if int(tx_length) < 100:
			continue
		exon_coords = tx_gene_coords[key][2:]
		gene_name = tx_gene_coords[key][0]
		strand = tx_gene_coords[key][1]
		chromosome = tx_gene_coords[key][2][0]
			
		# remove chr from exon coords for calculations
		exon_coords = [x[1:3] for x in exon_coords]
		# convert to integer
		exon_coords = [[int(pos) for pos in x] for x in exon_coords]
		exon_coords = sorted(exon_coords)
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
		offsets.pop(0)
		exon_coords.pop(0)
		# calculate new exon coords from 1 to len(transcript) continuously
		# can use the offets calculated to re-make the positions
		shifted_exon_coords = [[x[1][0]-offsets[x[0]],x[1][1]-offsets[x[0]]] for x in enumerate(exon_coords)] 
		
		######
		# build 100bp windows, shifted by 1bp and calculate number of variants in the window	
		######
		for i in range(1,int(tx_length)):
			windowDown = -50
			windowUp = 50
			# fix if window is outside transcript
			if i < windowDown:
				windowDown = -i
			if i > windowUp:
				windowUp = i-1
			interval_size = windowUp - windowDown
			# calc number of variants in the window
			overlapping_num = sum( [i+windowDown <= pos <= i+windowUp for pos in int_coding if pos] ) 
			# scale for interval size
			overlapping_num = str( round( overlapping_num * (100/interval_size) ) )
			# calculate actual coordinate by seeing which (shifted) exon i is in
			exon_offset_index = [x[0] for x in enumerate(shifted_exon_coords) if i in range(x[1][0],x[1][1])]
			#print(gene_name, strand, tx_length, exon_offset_index, i, shifted_exon_coords)
			# now can use the index, with the offset info to calc actual genomic position
			real_genomic_position = str(i + offsets[int(exon_offset_index[0])])
			print(chromosome + '\t' + real_genomic_position + '\t' + overlapping_num + '\t' + gene_name + \
					'\t' + key + '\t' + str(exon_coords[0][0]) + '\t' + str(exon_coords[-1][1]) + \
					'\t' + tx_length + '\t')
			
	
