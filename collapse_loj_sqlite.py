#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import re
import fileinput
import sqlite3

# create sqlite db
db = sqlite3.connect('/data/mcgaugheyd/projects/nei/mcgaughey/human_variation_landscape/var_info.sqlite')
#db = sqlite3.connect(':memory:')
c = db.cursor()
try:
	c.execute('''CREATE TABLE variation_landscape (ID TEXT PRIMARY KEY, Num_Var INTEGER, MAFs TEXT, VCF_INFO TEXT)''')
except:
	print('Table variation_landscape already exists')#db.commit()

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

	regex=re.compile(r'MAF=0\.\d+', re.I) # regex pattern to find the maf

	if key == window_key: # skip cases where the variants are the same 
		continue
	maf = 0 # to cover when no maf is present
	if regex.search(info):
		maf = [float(m.group().split('=')[1]) for section in info.split(';') for m in [regex.search(section)] if m][0]
	try:
		c.execute('INSERT INTO variation_landscape (ID, Num_Var, MAFs, VCF_INFO) VALUES (?,?,?,?)', \
			(key, 1, maf, vcf_info))
	except:
		old_data = c.execute('''SELECT * FROM variation_landscape WHERE ID=? ''',(key,)).fetchone()
		new_mafs = old_data[2] + ',' + str(maf)
		c.execute('''UPDATE variation_landscape SET MAFs=? WHERE ID=?''',(new_mafs,key))

db.commit()
#	else:
#		num_var, the_mafs, vcf_info = data[key]
#		num_var += 1
#		the_mafs.append(maf)
#		data[key] = num_var, the_mafs, vcf_info

