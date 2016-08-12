#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import re
import fileinput
import sqlite3

###
"""
Creates a sqlite file with several stats when run on the 'loj' files created from the bedtools 
intersect tool with -loj -sorted
"""

# create sqlite db
db = sqlite3.connect('/data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/var_info.sqlite')
#db = sqlite3.connect(':memory:')
c = db.cursor()
try:
	c.execute('''CREATE TABLE variation_landscape (ID TEXT PRIMARY KEY, \
				Num_Var INTEGER, MAFs TEXT, DISTANCEs TEXT, VCF_INFO TEXT)''')
except:
	print('Table variation_landscape already exists')

# read from stdin
for line in fileinput.input():	
	line = line.split('\t')	

	key = line[3]
	window_key = line[6]
	info = line[11]

	# grab info to build vcf at the end
	chr = line[4]
	pos = int( (int(line[1]) + int(line[2])) / 2)
	id = line[3]
	ref = line[7]
	alt = line[8]
	# lump together
	vcf_info = str(chr) + '\t' + str(pos) + '\t' + id + '\t' + ref + '\t' + alt + '\t.\t.'

	# calculate distance of landscape variant to primary variant
	orig_pos = int( (int(line[1]) + int(line[2])) / 2) # find center to get original position
	distance = orig_pos - int(line[5])

	if key == window_key: # skip cases where the variants are the same 
		continue
	
	maf = 0 # to cover when no maf is present
	regex=re.compile(r'MAF=0\.\d+', re.I) # regex pattern to find the maf
	if regex.search(info):
		maf = [float(m.group().split('=')[1]) for section in info.split(';') for m in [regex.search(section)] if m][0]
	try:
		c.execute('INSERT INTO variation_landscape (ID, Num_Var, MAFs, DISTANCEs, VCF_INFO) \
			VALUES (?,?,?,?,?)', \
			(key, 1, maf, distance, vcf_info))
	except:
		old_data = c.execute('''SELECT * FROM variation_landscape WHERE ID=? ''',(key,)).fetchone()
		new_mafs = old_data[2] + ',' + str(maf)
		new_distances = old_data[3] + ',' + str(distance)
		c.execute('''UPDATE variation_landscape SET MAFs=?,DISTANCEs=? WHERE ID=?''',(new_mafs,new_distances,key))

db.commit()
#	else:
#		num_var, the_mafs, vcf_info = data[key]
#		num_var += 1
#		the_mafs.append(maf)
#		data[key] = num_var, the_mafs, vcf_info

