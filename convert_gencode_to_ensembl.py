#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import fileinput

# https://github.com/dpryan79/ChromosomeMappings
# https://github.com/dpryan79/ChromosomeMappings/blob/master/GRCm38_gencode2ensembl.txt
converter=open('/home/mcgaugheyd/git/human_variation_landscape/GRCh38_gencode2ensembl.txt')

converter_dict = {}
for line in converter:
	line = line.split()
	converter_dict[line[0]] = line[1]

for line in fileinput.input():
	line = line.split()	
	chr = line[0]
	line[0] = converter_dict[chr]
	print('\t'.join(line)) 
