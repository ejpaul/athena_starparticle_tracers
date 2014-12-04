#!python
import numpy as np
"""
  handling history data files from Athena
"""
def read(filename,sortable_key=False,silent=False,headerfile=None):

	import re
	import os

	if not silent: print "Reading a history file:" + filename
	base=os.path.basename(filename)
	split=re.split('\.',base)
	ext=split[-1]
	file=open(filename,'r')
	lines=file.readlines()
	file.close()

	data={}

	header1 = lines[0]
	header = lines[1]
	if headerfile != None:
		ftmp = open(headerfile,'r')
		header1 = ftmp.readline()	
		header = ftmp.readline()	
		file.close()

	if ext == 'hst':
		match=re.search('-?\d(\.\d+|)[Ee][+\-]\d\d',header1)
		data['vol']=eval(match.group(0))

	if header[0] != "#":
		split=re.findall('[+\-]?\d\.?\d*[Ee][+\-]\d\d?',header)
		nvar=len(split)
		varlist=[]
		for i in range(nvar):
			varlist.append('var')
	else:
		varlist=re.split("\[\d+]\=|\n",header)
		for i in range(len(varlist)): varlist[i]=re.sub("\s|\W","",varlist[i])
		varlist=varlist[1:-1]
		nvar=len(varlist)

	if sortable_key:
		for i in range(nvar): 
			head= '%02d' % (i+1)
			varlist[i] = head+varlist[i] 

	for var in varlist:
		data[var] = []

	for line in lines:
		split=re.findall('[+\-]?\d\.?\d*[Ee][+\-]\d\d?',line)
	 	if nvar == len(split):
			for var, value  in zip(varlist, split):
				data[var].append(eval(value))

	for var in varlist:
		data[var] = np.array(data[var])

	return data
