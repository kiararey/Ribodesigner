## Bacillus_halotolerans.fasta

first 50 lines of your example file

## fasta_to_jsonindex

converts the fasta file to a json index file

it assumes that the file alternates between header lines and actg lines

in the future you may want to do away with these obvious indices and just have tab-delimited files --they'll be more efficient to read in one line at a time

but right now they're not too long -- and you can see how you'd go about generating semantic metadata on these for your json output file

run it like ```python fasta_to_jsonindex.py Bacillus_halotolerans.fasta```

## the json files

### the sequence dictionaries

these represent your sample test and comparison datasets -- i don't know precisely what they look like.

### the rules dictionary (rules.json)

this is generated by running rulegenerator.py -- which you'll have to put your own logic into, obviously.

so for instance you would run this one like ```python rulegenerator.py generated.json Bacillus_halotolerans.json```

## parallelize.py

this takes a test file, a comparison file, a ruleslist file, a worker index, and a number of workers, to parallelize as a job array.

given all of the above setup steps, this would look like ```python parallelization.py generated.json Bacillus_halotolerans.json rules.json 2 50``` to run worker #3 on a 50-worker job array.

## parallelize.slurm

finally, you can now run a slurm script that parallelizes that python wrapper script with a job array, and passes in the count of workers and the individual worker index as arguments

an example is provided that gives you 50 workers for 4 hours (so, 200 CPU hours) on the scavenge queue

but you'll want to work out how to install your requisite software other than python.