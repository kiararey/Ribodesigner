# RiboDesigner

### Setting up

I've finally uploaded sample data to use, I've added this in the Ribodesigner_test_files folder.

To use, unzip this folder and extract the ``test_output_files`` directory into the same folder 
where you'll run ribodesigner on.

Make sure you have Muscle5 installed in your path! More info found here: https://www.drive5.com/muscle/

This program makes use of Biopython, numpy, pandas, seaborn, and matplotlib, playsound, date_util, icecream, alive_progress. I think the rest of the dependencies are
default libraries but please let me know if this is not the case. I'll make an environment file later too

### Running sample data
1. Edit main.py to run ``ribodesigner_routine`` with test data
2. Once test data is done, run ``run_local`` to analyze data
3. Exctract data to later analyze with ``import_data_to_df``

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