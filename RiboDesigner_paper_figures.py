# Figure 2: synthetic community data violin plots
import Figure_plotting_functions as fig_plt
import os
import numpy as np

# Basically, look at the datasets used data and then run RiboDesigner on them. Generate a violin plot using the results.
# Run RiboDesigner on all datasets we are looking at
# First, set up base files and parameters
m = 5
n = 50
minlen = 50

# Barcode sequence is split sfGFP just cuz. This does not affect guide sequence design.
barcode_seq_file = 'Common_sequences/sfGFP_2_seq_barcode.txt'

# We'll be using the ribozyme published in the RAM paper
ribobody_file = 'Common_sequences/ribozyme_body.txt'

# Prepare the datasets
datasets_path = 'Datasets_used/zymo_files/'

# Output folder
output_path = 'Figure_2_output_files/'

# Reference sequence - will have to be E. coli to graph the variable regions
ref_path = 'Common_sequences/e-coli-16s-mg1655.fasta'

########################################################
# Score vs. True Coverage graphs
# Generate data here. If already generated we'll just import it with pandas
# For each folder in datasets used run RiboDesigner and keep the data for plotting later
# get rid of that pesky .DS_Store file
datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
datasets.sort()
ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True]
fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, None)
fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)


########################################################
# SILVA squished datasets
datasets_path = 'Datasets_used/SILVA_squished_datasets'
output_path = 'SILVA_figure_output_files/'

datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
datasets.sort()
ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True]
fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)
fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)
