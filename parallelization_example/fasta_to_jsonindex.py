import json
import sys
import re

fname=sys.argv[1]

if not fname.endswith('.fasta'):
	print("INPUT NEEDS TO BE A FASTA FILE")
	exit()

d=open(fname,'r')
lines=d.readlines()
d.close()

fastadict={}

headerline=None
sequenceline=None

c=0
for line in lines:
	if line.startswith('>'):
		isheaderline=True
		sequenceline=None
		headerline=line
	if re.match("[C|A|T|G]",line):
		isheaderline=False
		if headerline is not None:
			sequenceline=line
	if headerline is not None and sequenceline is not None:
		fastadict[c]={'sequence':sequenceline,'headerdata':headerline}
		c+=1

fnamestripped=re.sub("\.fasta","",fname)

d=open(fnamestripped+'.json','w')
d.write(json.dumps(fastadict,indent=2))
d.close()