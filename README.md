# RiboDesigner

### Setting up

I've finally uploaded sample data to use, I've added this in the Ribodesigner_test_files folder.

To use, unzip this folder and extract the ``Datasets_used`` and ``test_output_files`` directories into the same folder 
where you'll run ribodesigner on.

Make sure you have Muscle5 installed in your path! More info found here: https://www.drive5.com/muscle/

This program makes use of Biopython, numpy, pandas, seaborn, and matplotlib. I think the rest of the dependencies are
default libraries but please let me know if this is not the case.

To run, I've honestly been opening main.py on Pycharm. if you're running this for the first time I highly suggest
replacing all instances of ``test_data_folders`` with ``test_data_folders_test``. I also strongly suggest to set all 
``percent_of_target_seqs_used``, ``percent_of_target_seqs_used`` and ``percent_of_background_seqs_used`` to something 
like 0.01 (except for the reference sequence and control sequences, those only have one target).

### Troubleshooting

IF MUSCLE STOPS WORKING!!!!!

- Install [homebrew]( https://brew.sh/) 

``/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"``  

- Install [gcc@11](https://gcc.gnu.org/) 

``brew install gcc@11``  

**In finder** (or in terminal if you are more computer literate or have more patience than me):  
- Find ``/opt/homebrew/opt/gcc/lib/gcc``  
- Open a separate tab in finder.  
- In this separate tab, find ``/opt/homebrew/opt/gcc@11/lib/gcc``  
- Duplicate the file called 11 in this separate tab and move that to ``/opt/homebrew/opt/gcc/lib/gcc``  

- Run muscle5 and check that everything works properly!