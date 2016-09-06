#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import fileinput
import re

header = ['Chr','Position','ID','Ref','Alt','Variant_Class','GMAF','ExAC',\
			'CADD_phred','CADD_raw','CADD_raw_rankscore',\
			'PolyPhen', 'SIFT', 'Clin_Sig', \
			'VLN100all','VLN100high','VLN100moderate','VLN100low']

out_file=open('Homo_sapiens_incl_consequences__codingOnly.VLNum.VEPnoPick.GRCh38.20160906.my.tab','w')
out_file.write('\t'.join(header))
for line in fileinput.input():
	if line[0]=='#':
		continue
	line = line.split('\t')
	chr = line[0]; position = line[1]; id = line[2]
	ref = line[3]; alt = line[4]
	info = line[7]
	info = info.split(';')
	try:
		vln = re.compile('VLN100')
		vln_pos = [x[0] for x in enumerate(info) if vln.match(x[1])]
		vln100all = info[vln_pos[0]].split('=')[1]; vln100high = info[vln_pos[1]].split('=')[1]
		vln100mod = info[vln_pos[2]].split('=')[1]; vln100low = info[vln_pos[3]].split('=')[1]
	
		clin = re.compile('^CLIN_')
		clin_pos = [x[0] for x in enumerate(info) if clin.match(x[1])]
		clin_sig = info[clin_pos[0]]

		csq = re.compile('^CSQ=')
		vep_info_index = [x[0] for x in enumerate(info) if csq.match(x[1])][0]
		vep_info = info[vep_info_index]
		variant_class = vep_info.split('|')[1]
		try:
			gmaf = vep_info.split('|')[27].split(':')[1].split('&')[0]
		except:
			gmaf = ''
		try:
			exac = vep_info.split('|')[36].split(':')[1].split('&')[0]
		except:
			exac = ''
		cadd_phred = vep_info.split('|')[-3]
		cadd_raw = vep_info.split('|')[-2]
		cadd_raw_rankscore = vep_info.split('|')[-1][:-1]
		sift = vep_info.split('|')[24]
		polyphen = vep_info.split('|')[25]
	except:
		continue
	out = 	chr, position, id, ref, alt, variant_class, gmaf, exac, \
			cadd_phred,cadd_raw, cadd_raw_rankscore, polyphen, sift, \
			clin_sig, vln100all,vln100high, vln100mod, vln100low
	out_file.write('\t'.join(out))

out_file.close()
