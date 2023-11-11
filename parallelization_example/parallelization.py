import sys
import json
import os
import random
import numpy as np

def pairwise_compare(test_data,comparison_data,worker):
	
	##yourlogichere
	
	output=random.choice([True,False])
	
	return output


def main(testfilename,comparisonfilename,worker_number,rules_file,number_of_workers):
	
	d=open(testfilename,'r')
	t=d.read()
	d.close()
	testdata=json.loads(t)
	
	d=open(comparisonfilename,'r')
	t=d.read()
	d.close()
	comparisondata=json.loads(t)
	
	d=open(rules_file,'r')
	t=d.read()
	d.close()
	rules=json.loads(t)
	
	#load the completed work, which we'll store as a tsv rather than as a json file
	checkpointfiles=[f for f in os.listdir('checkpoints/') if f.startswith('checkpoint')]
	completed_work=[]
	for checkpointfile in checkpointfile:
		d=open(comparisonfilename,'r')
		t=d.read()
		d.close()
		lines=[line for line in t.split('\n') if line!='']
		for line in lines:
			this_completed_work=[int(i.strip()) for i in line.split('\t')]
			completed_work.append(this_completed_work)
	
	fullworklist=[]
	
	for test_idx in test_data:
		for comparison_idx in comparison_data:
			if [test_idx,comparison_idx] in rules:
				test_compare_idx_tuple=[test_idx,comparison_idx]
				fullworklist.append(test_compare_idx_tuple)
	
	print("full work list is %d items long" %len(fullworklist))
	
	finishedworklist=[item for item in fullworklist if item not in completed_work]
	
	print("remaining work list is %d items long" $len(finishedworklist))
	
	#now divvy up the work amongst the number of workers and have each worker pull its slice
	numpyarrayworklist=np.array(finishedworklist)
	thisworklist=array.split(numpyarrayworklist,worker_number,number_of_workers)[worker_number]
	
	for tup in thisworklist:
		test_idx,comparison_idx=tup
		this_comparison_data=comparisondata[comparison_idx]
		this_test_data=testdata[test_idx]
		
		#get and then record the result -- in an albeit messy format, as lines of json dumps
		#don't checkpoint until you've ensured that your work successfully completed
		pairwise_result=pairwise_compare(this_test_data,this_comparison_data)
		
		#but immediately write that work out and checkpoint it as closely together as possible
		#so that you don't accidentally record finished work but fail to checkpoint it, and thereby duplicate results
		d=open('outputs/%d_pairwise.txt' %worker_number ,'a')
		d.write(
			json.dumps({
			'test_idx':test_idx,
			'comparison_idx':comparison_idx,
			'result':pairwise_result
			} + '\n'
		))
		d.close()
		
		d=open('checkpoints/%d_checkpoints.txt' %worker_number,'a')
		d.write('\t'.join([str(i) for i in [test_idx,comparison_idx]])+'\n')
		d.close()

if __name__=="__main__":
	testfilename=sys.argv[1]
	comparisonfilename=sys.argv[2]
	worker_number=int(sys.argv[2])
	rules_file=int(sys.argv[3])
	number_of_workers=int(sys.argv[4])
	main(testfilename,comparisonfilename,worker_number,rules_file,number_of_workers)