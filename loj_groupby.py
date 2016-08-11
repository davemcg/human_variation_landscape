#!/usr/local/Anaconda/envs/py3.4.3/bin/python

from itertools import groupby
import fileinput
from operator import itemgetter

for key, line in groupby(fileinput.input(), lambda x: x.split()[3]):
	print(key)
	for sublines in line:
		print(sublines)
	print('-'*20)	
