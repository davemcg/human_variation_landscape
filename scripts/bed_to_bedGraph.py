#!/usr/local/Anaconda/envs/py3.4.3/bin/python


"""
Groups sorted bed lines based on the value column (4th one). Will merge start and stop 
from. Example:
1 1000 1100 43.2
1 1100 1200 43.2
1 1200 1300 23.2
Will be convered to:
1 1000 1200 43.2
1 1200 1300 23.2
"""

from itertools import groupby
import fileinput

for key, chunk in groupby(fileinput.input(), lambda x: x.split()[0] + x.split()[4]):
	chunk = list(chunk)
	chr = chunk[0].split()[0]
	start = chunk[0].split()[1]
	stop = chunk[-1].split()[2]
	new_key = chr + '_' + start + '_' + stop
	value = chunk[0].split()[4]
	print(chr, start, stop, new_key, value)
