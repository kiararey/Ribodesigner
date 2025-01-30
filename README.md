# RiboDesigner

### Setting up
notes from Kiara:
- Install [muscle5](https://www.drive5.com/muscle/) - make sure to follow their installation instructions!
- Download python 3.11
- Download all files, unzip
- open the folder in your favorite IDE or code editor (I like opening mine as a project in PyCharm)
- Extract ``Ribodesigner_test_files`` and put ``test_output_files`` in the same folder
- Open ``main.py`` and follow running sample data


I've finally uploaded sample data to use, I've added this in the Ribodesigner_test_files folder.

To use, unzip this folder and extract the ``test_output_files`` directory into the same folder 
where you'll run ribodesigner on.

Make sure you have Muscle5 installed in your path! More info found here: https://www.drive5.com/muscle/

This program makes use of Biopython, numpy, pandas, seaborn, and matplotlib, playsound, date_util, icecream, alive_progress. I think the rest of the dependencies are
default libraries but please let me know if this is not the case. I'll make an environment file later too

### Running sample data
1. Edit ``main.py`` to run ``ribodesigner_routine`` with test data
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

### What do the columns headers on my data mean?
- **id** (str) : A unique identifier for each cat-RNA design including the IGS sequence and the location of the catalytic U (U-IGS) it targets using the reference sequence indexing. Each U-IGS at a specific location only has one design.
- **igs** (str) : The IGS sequence used by the cat-RNA design.
- **reference_idx** (int) : The location of the catalytic U the cat-RNA design is targeting using the reference sequence numbering. For example, 16 would mean that that U occurs at position 16 on the reference sequence used.
- **optimized_to_targets** (bool) : Whether the optimization step was performed on the guide sequence for the cat-RNA design. Check out our publication for more details.
- **tested** (bool) : Whether the design was tested against a set of test sequences.
- **guide** (str) : The guide sequence used by the cat-RNA design.
- **full_guide_g_igs** (str) : The guide, g, and IGS in the correct order used by this design.
- **num_of_targets** (int) : The number of target sequences that had a U-pentanucleotide matching the IGS of the cat-RNA design at the location targeted.
- **score_type** (str) : The type or scoring method used to score the guide. For more information on scoring method options check out our publication.
- **score** (float) : The guide score for the cat-RNA design. A 0 means there is only ambiguity in the design, while a 1 means that there is no ambiguity (only A, T, G, C).
- **U_IGS_coverage** (float) : The fraction of target sequences that have a U-pentanucleotide complementary to the IGS at any location out of all target sequences. (number of target sequences with a U-pentanucleotide / total number of target sequences)
- **U_IGS_on_target** (float) : The fraction of target sequences with a U-pentanucleotide complementary to the IGS where that U is at the correct targeted location out of sequences that have a matching U-pentanucleotide at any position. (number of target sequences with a matching U-pentanucleotide at the correct location / number of target sequences with a U-pentanucleotide)
- **true_U_IGS_coverage** (float) : The fraction of target sequences with a U-pentanucleotide complementary to the IGS at the correct targeted location out of all target sequences. It's the multiplication of U-IGS-on-target and U-IGS-coverage. (number of target sequences with a matching U-pentanucleotide at the correct location / total number of target sequences)
- **composite_score** (float) : The guide score * true U-IGS coverage score (score * true_U_IGS_coverage)
- **name_of_test_dataset** (str) : The name of the fasta file used for the test dataset.
- **num_of_targets_test** (int) : The number of test sequences that had a U-pentanucleotide matching the IGS of the cat-RNA design at the location targeted.
- **u_conservation_test** (float) : The fraction of test sequences that have a U at the correct position targetted by the cat-RNA. 
- **test_score** (float) : The average guide score for the cat-RNA design when evaluated against the test sequences. A 0 means there is no conservation on average between test sequences and the cat-RNA guide sequence at the location targetted and 1 means that there is perfect conservation with all of them. *(caveat below)*
- **tm_nn_vs_test** (float) : The average melting temperature of the guide with the possible binding sequences in the test sequences. Uses GC content to evaluate. It's calculated but I would not recommend using this as we did not go further with this type of data. *(caveat below)*
- **U_IGS_coverage_test** (float) : The fraction of test sequences that have a U-pentanucleotide complementary to the IGS at any location out of all test sequences. (number of test sequences with a matching U-pentanucleotide / total number of test sequences)
- **U_IGS_on_target_test** (float) : The fraction of test sequences with a U-pentanucleotide complementary to the IGS where that U is at the correct targeted location out of sequences that have a matching U-pentanucleotide at any position. (number of test sequences with a matching U-pentanucleotide at the correct location / number of test sequences with a matching U-pentanucleotide)
- **true_U_IGS_coverage_test** (float) : The fraction of test sequences with a U-pentanucleotide complementary to the IGS at the correct targeted location out of all test sequences. It's the multiplication of U_IGS_on_target_test and U_IGS_coverage_test. (number of test sequences with a matching U-pentanucleotide at the correct location / total number of test sequences)
- **composite_test_score** (float) : The test guide score * test true U-IGS coverage score (test_score * true_U_IGS_coverage_test)
- **delta_guide_vs_test** (float) : true_U_IGS_coverage - true_U_IGS_coverage_test
- **delta_vs_test** (float) : test_score - score

The last few columns are divided by order of taxonomy and are only relevant if you are using SILVA database formatted headers. They contain dictionaries with the taxonomic group targetted (either target sequences (e.g. target_Domain) or test sequences (e.g. target_Domain_test) with the number of sequences within that group. For example, if target_Domain_test has {'Archaea':53}, that means that the cat-RNA found 53 sequences in test dataset with a matching U-IGS for the domain Archaea.  *(caveat below)*

caveat: *By default all sequences in the test dataset with a U (not necessarily a U-pentanucleotide) at the correct position are evaluated against the designed guide sequence, but this can be changed by setting ``flexible_igs`` in ``ribodesigner_routine`` to False. This will make it so that only sequences with an exact U-pentanucleotide match in the correct location will be evaluated. This was done because we were exploring the parameter space in the paper but it might not be relevant for you.*
