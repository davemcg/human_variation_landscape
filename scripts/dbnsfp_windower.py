#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
from itertools import dropwhile
from itertools import islice
import fileinput
import sys

def central_depot(dbnSFP_file):
	# grab header
	top_line = []
	for line in islice(dbnSFP_file, 10):
		top_line.append(line[:-1])
	# grab header and parse for various column names
	header = ([x for x in top_line if x[0:4]=='#chr'])[0]
	codon_pos = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'aapos'])[0]
	CADD_raw_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'CADD_raw_rankscore'])[0]
	gerp_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'GERP++_RS_rankscore'])[0]
	gmaf_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == '1000Gp3_AF'])[0]
	exac_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'ExAC_AF'])[0]
	clinsig_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'clinvar_clnsig'])[0]
	tx_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'Ensembl_transcriptid'])[0]
	gene_index = ([x[0] for x in enumerate(header.split('\t')) if x[1] == 'Ensembl_geneid'])[0]
	# skip header lines
	for line in dropwhile(lambda line: line.startswith('#'), dbnSFP_file):	
		# group on gene name. dbnSFP MUST BE PRE-SORTED on the gene name (20th column)
		# if not, then things will go TERRIBLY WRONG!!!!!!!!!!!!!!
		for key, chunk in groupby(dbnSFP_file, lambda x: x.split()[gene_index]):
			chunk = list(chunk)
			# get list of all transcripts
			try:
				transcripts = [x.split('\t')[tx_index].split(';') for x in chunk if ';' in x]
				unique_tx = list(set([item for sublist in transcripts for item in sublist]))
			except:
				unique_tx = list(set([x.split('\t')[tx_index] for x in chunk]))
			for tx in unique_tx:
				print(key,tx)
				# grab tx specific positions
				tx_indices = [x[0] for x in enumerate(chunk) if tx in x[1].split('\t')[tx_index]]
				tx_data = [x[1] for x in enumerate(chunk) if x[0] in tx_indices]
				# grab relevant aa position
				# the transcripts (tx) are collapsed on each position, split by ;
				# as we roll through the genome/gene, the transcripts change
				# so we have to check each position and see which position the 
				# transcript is in so we get the correct aa position, which is also split by ;
				transcripts_list_of_list = [x.split('\t')[tx_index].split(';') for x in tx_data] 
				pos_of_tx = [] 
				for sub_list in transcripts_list_of_list:
					pos_of_tx.append([x[0] for x in enumerate(sub_list) if tx in x[1]][0])
				tx_aa_pos_list_of_list = [x.split('\t')[codon_pos].split(';') for x in tx_data]
				actual_aa_pos = []
				for sub_list in tx_aa_pos_list_of_list:
					if len(sub_list)==1:
						actual_aa_pos.append(int(sub_list[0]))
					else:
						actual_aa_pos.append([int(x[1]) for x in enumerate(sub_list) if pos_of_tx[x[0]]==x[0]])
				try:
					len(pos_of_tx) == len(actual_aa_pos)
				except:
					print('Warning, mismatch of transcript index and aa for ' + key, chunk)
					sys.exit(1)
				tx_cadd_rankscores = [x.split('\t')[CADD_raw_index] for x in tx_data]
				tx_gerp_rankscores = [x.split('\t')[gerp_index] for x in tx_data]
				print(tx_aa_pos_list_of_list)
				print(transcripts_list_of_list)
				print(pos_of_tx)
				print(actual_aa_pos)
				print(tx_cadd_rankscores)
				print(tx_gerp_rankscores)
				print(str(len(pos_of_tx))+' '+str(len(actual_aa_pos))+' '+str(len(tx_cadd_rankscores))+' '+str(len(tx_gerp_rankscores)))
				#[print(x.strip()) for x in tx_data]	

def main():
	central_depot(fileinput.input())

main()
