# Step by step to use the cluster  
#### Last edited 2024.02.01 KRG

**What is the cluster and why to use it:** NOTS is Rice's supercomputing cluster. If your process takes a ridiculously 
long time to run and uses up a lot of resources on your own local computer (eg. say you can't use your computer for 
hours or days on end if you run your dataset) it may be time to run it on the cluster! You can set up your program to
basically divide up the work into smaller chunks that will be assigned to separate resources on NOTS. Depending on how
beefy your computer is each resource that was allocated might not be that much faster, but if you can divide up the
work properly your task could finish a lot faster by running it in parallel on the cluster.  

Alternatively, say you've vectorized and/or parallelized your code properly to run really fast on your local computer
but it's using all of the memory and processing power your computer has, making it unusable for the amount of time it 
needs to run your job. Even if you don't get a good speed gain, running your job on the cluster is still advantageous
as you can run your job remotely while maintaining your computer's functionality.

An intro to NOTS can be found in https://researchcomputing.rice.edu/rice-supercomputing-nots  
More details can be found in https://kb.rice.edu/page.php?id=108237

## Before you start!
**Go to OIT before you begin!!!** They're super nice and have folks that can help you set up your remote files properly
and might even help you brainstorm a way to make your code parallel and maybe write a script to run it on the cluster 
properly. You'll want to be ready to explain your project and show your code, especially the parts that take longer to 
run. You'll also want to make sure to have a git repository set up with your code, ideally with an easy to copy file 
structure so you can copy files in and out of the cluster easily.
Make sure to ask them to help you set up the following:
- Key to enter NOTS (that's what the `eval` and `ssh-add` commands will do later)
- Permissions to clone the git repository where you have your code
- Make a slurm script (I made two for my own project: a test one `parallelize_test.slurm` 
and a real one`parallelize.slurm`)
- Help with setting up your first run!  

I also recommend setting up Globus to transfer files back and forth from the cluster - more details in the 
"Upload files to test" section

## Preparing to test
if off campus, can enter NOTS through the backdoor by going to:  
`ssh <netID>@gw.crc.rice.edu`  
Then use your computer password to allow access, and then continue as if you were testing from on-campus.
Log into NOTS:  
`ssh <netID>@nots.rice.edu`  
Add key to use:  
`eval $(ssh-agent)`  
`ssh-add NOTS`
Pull the folder where you have your repo  
`cd <repo folder name>/`  
`git pull`  
Check you’re on the right branch  
`git branch`  
Check you’re on the right commit  
`git log` (quit by pressing shift + Q)  
switch to scratch to simulate running – should change  
`cd /scratch/<netID>/<repo folder name>`  
Copy all files that are not in a folder  
`cp -R /home/<netID>/<repo folder name>/* .`  
Check that slurm script is good  
`cat parallelize_test.slurm`
The script should start something look like:  
```pycon
module purge
module load GCCcore/12.3.0  
module load Python/3.11.3  
pip install --user -r requirements.txt 
```
Test now  
`sbatch parallelize_test.slurm`  

#### AFTER TESTING DELETE ANY OUTPUT FILES!!!    
exit session  
`exit`

Here is a one line slurm script to test a script in interactive mode. More info can be found
https://kb.rice.edu/108437  
`srun --partition=interactive --pty --export=ALL --ntasks=1 --time=00:30:00 /bin/bash`

## Upload files to test
The easiest way to do this is to set up Globus. Here's a page on how to do this:
https://kb.rice.edu/page.php?id=108242  
Otherwise, follow these instructions: 
First, log onto Rice VPN  
Run this command to copy files over  
`scp -r <local folder where you have the files you want to copy over>* <netID>@nots.rice.edu:<path to folder where you want to copy the files to>`

## Notes
cmd + K (connect to server)  
smb://smb.rdf.rice.edu/research

## Running a job
#### REMEMBER TO RUN ALL JOBS ON SCRATCH!!!
`cd /scratch/<netID>/<repo folder name>`  
Check that slurm script is good  
`cat parallelize.slurm`  
Submit batch  
`sbatch parallelize.slurm`  
Check queue  
`squeue -u <netID>`  
Check the output from a run  
`cat <.out file name>`  
Move results file into your results folder (or use Globus ideally)  
`cp -R /scratch/<netID>/<path to results folder>`  
Compress file  
`tar -czvf results.tgz results`  

## Getting results back
In separate terminal, now download results  
Move to folder where you want results to be  
`cd <path to folder>`  
Download compressed results folder  
`scp -r <netID>@nots.rice.edu:./results.tgz .`  
You can now exit your ssh session and analyze your files locally!
