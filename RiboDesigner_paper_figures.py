# Figure 2: synthetic community data violin plots
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import pandas as pd
from RiboDesigner import RiboDesigner, read_fasta_folder
import os
import math
import cmcrameri.cm as cmc
import numpy as np

# Basically, look at the datasets used data and then run RiboDesigner on them. Generate a violin plot using the results.
# Run RiboDesigner on all datasets we are looking at
# First, set up base files and parameters
m = 5
n = 50
minlen = 50

# Barcode sequence is split sfGFP just cuz. This does not affect guide sequence design.
barcode_seq_file = 'Barcode_and_ribozyme_body_sequences/sfGFP_2_seq_barcode.txt'

# We'll be using the ribozyme published in the RAM paper
ribobody_file = 'Barcode_and_ribozyme_body_sequences/ribozyme_body.txt'

# Prepare the datasets
datasets_path = 'Datasets_used/'

# Output folder
output_path = 'Figure_2_output_files/'

# Generate data here. If already generated we'll just import it with pandas
# For each folder in datasets used run RiboDesigner and keep the data for plotting later
# get rid of that pesky .DS_Store file
datasets = np.array([set for set in os.listdir(datasets_path) if set != '.DS_Store'])
datasets.sort()
data = [None] * len(datasets)
# Set up subplots

plt.figure()

custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30, 16)}
cmap = np.append([cmc.batlow.colors[0]], [cmc.batlow.colors[-1]], axis=0)
sns.set_palette(palette=cmap)
sns.set_theme(context='talk', style="ticks", rc=custom_params, palette=cmap)

fig, axs = plt.subplots(math.ceil(len(datasets) / 3), 3, sharex=True, sharey=True, layout="constrained")

for i in range(0, len(datasets)):

    dataset_name = datasets[i].replace('_', ' ')

    output_path_folder = output_path + datasets[i] + '_results'
    target_sequences_folder = datasets_path + datasets[i]
    org_nums = len(read_fasta_folder(target_sequences_folder))

    if os.path.isdir(output_path_folder) is False:
        os.mkdir(output_path_folder)

        out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, target_sequences_folder, min_true_cov=0,
                                identity_thresh=0.7, fileout=True, folder_to_save=output_path_folder)
        out_data_df = pd.DataFrame(data=out_data, index=None, columns=['IGS', 'Reference index', 'Score', '% cov',
                                                                       '% on target', 'True % cov',
                                                                       '(Target name, Target idx, Other occurrences of'
                                                                       ' IGS in target sequence)', 'Optimized guide',
                                                                       'Optimized guide + G + IGS',
                                                                       'Full ribozyme design'],
                                   dtype=object).sort_values(by=['True % cov', 'Score'], ascending=[False, False])

    else:

        output_path_file = output_path + datasets[i] + \
                           '_results/Ranked Ribozyme Designs with Optimized Guide Sequence Designs quantitative.csv'
        out_data_df = pd.read_csv(output_path_file)

    colors = [None] * len(out_data_df.loc[:, ['True % cov']])

    for j in range(0, len(colors)):
        row = out_data_df.loc[j, ['True % cov', 'Score']]
        if row[0] > 0.7 and row[1] > 0.7:
            colors[j] = True
        else:
            colors[j] = False

    perc = sum(colors) / len(colors) * 100

    # relevant data here: dataset name, percent above threshold, species number, ribodesigner output
    data[i] = (dataset_name, perc, org_nums, out_data_df, colors)

    ax_r_coord = math.floor(i / 3)
    if math.floor(i / 3) == math.ceil(i / 3):
        ax_c_coord = 0
    elif math.floor(i / 3) == round(i / 3):
        ax_c_coord = 1
    else:
        ax_c_coord = 2

    axs[ax_r_coord, ax_c_coord].set_title(dataset_name + '\nn = ' + str(org_nums) + ', ' + str(round(perc, 2)) +
                                          '% designs above threshold.')
    sns.scatterplot(data=out_data_df, x='True % cov', y='Score', alpha=0.8, hue=colors, legend=False,
                    ax=axs[ax_r_coord, ax_c_coord])
sns.despine()
plt.savefig(output_path + '9 panel figure.png', transparent=False)
plt.show()

# # Now make a final plot comparing number of species vs. percentage above threshold
custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (12, 8)}
sns.set_theme(context='talk', style="ticks", rc=custom_params, palette=cmap)
plt.figure()

plt.title('% above threshold vs. number of species')

x_to_plot = [data[i][1] for i in range(len(data) - 1)]
y_to_plot = [data[i][2] for i in range(len(data) - 1)]
sns.scatterplot(x=x_to_plot, y=y_to_plot, alpha=0.8)
sns.lineplot(x=x_to_plot, y=y_to_plot, alpha=0.8)
plt.xlabel('Percent of designs above threshold')
plt.ylabel('Number of species in dataset')
sns.despine()

plt.savefig(output_path + 'species vs percent above threshold.png', transparent=False)
plt.show()
