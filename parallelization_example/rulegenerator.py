import json
import sys
import random

testfilename=sys.argv[1]
comparisonfilename=sys.argv[2]

d=open(testfilename,'r')
t=d.read()
d.close()
test_data=json.loads(t)

d=open(comparisonfilename,'r')
t=d.read()
d.close()
comparison_data=json.loads(t)

def shouldthesebeprocessedtogether(this_test_data,this_comparison_data):
	
	## logic
	
	return random.choice([True,False])
	


ruleset={}

for test_idx in test_data:
	this_test_data=test_data[test_idx]
	if test_idx not in ruleset:
		ruleset[test_idx]={}
	for comparison_idx in comparison_data:
		this_comparison_data=comparison_data[comparison_idx]
		
		decision=shouldthesebeprocessedtogether(this_test_data,this_comparison_data)
		
		ruleset[test_idx][comparison_idx]=decision

d=open('rules.json','w')
d.write(json.dumps(ruleset,indent=2))
d.close()