# Figure 2: synthetic community data violin plots
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from RiboDesigner import RiboDesigner
import os

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

# For each folder in datasets used run RiboDesigner and keep the data for plotting later
data = [None] * len(os.listdir(datasets_path))
col = 0

for dataset in os.listdir(datasets_path):
    output_path_folder = output_path + dataset + '_results'

    if os.path.isdir(output_path_folder) is False:
        os.mkdir(output_path_folder)

    target_sequences_folder = datasets_path + dataset
    data[0] = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, target_sequences_folder, min_true_cov=0.7,
                            identity_thresh=0.7, fileout=True, folder_to_save=output_path_folder)
    col += 1


# Feed the data into plotting function
make_split_plot(data, output_path + 'Figure 2 plot.png')

# Violin plot should look like this: x-axis is the dataset name, y-axis is the score, split between 70% minimum
# true coverage and 100% true coverage.



def make_split_plot(target_data, title, filename, names, var_regs, colr, splt=False, left_data=[], right_data=[]):



    fig, ax = plt.subplots(2, 1, sharey='col', constrained_layout=True)

    names = [tup[0].replace('_', ' ') for tup in target_data]
    refidx = [tup[1] for tup in target_data]
    IGS = [tup[2] for tup in target_data]
    score = [tup[3] for tup in target_data]
    repeat_cols = [tup[4] for tup in target_data]
    off_target = [tup[5] for tup in target_data]

    plot_var_reg(ax[0], var_regs)

    ax[0].set(xlim=(-0.1, 1599), ylim=(0, 1.1))

    # plot score vs. site of U (ref--> E.coli)
    ax[0].scatter(refidx, score, alpha=0.5, color=repeat_cols, edgecolors=off_target)
    ax[0].set_ylabel('Score')
    ax[0].set_xlabel('Location of U on E. coli 16s rRNA (base pair location)')
    ax[0].set_title(title)
    sns.despine(ax=ax[0], offset=5)

    # prepare names for plotting
    new_names = [None] * len(names)
    col = 0

    for name in names:
        index = name.find(' ')
        new_name = name[0] + '.' + name[index:]
        new_names[col] = new_name

    # plot score vs. target name
    sns.set_palette('viridis')

    plot_var_reg(ax[1], var_regs)

    ax[1].set_ylabel('Score')
    ax[1].set_xlabel('Target name')
    ax[1].tick_params(bottom=False)
    if splt:
        scoreOn = [tup[3] for tup in left_data]
        namesOn = [fixNames(tup[0].replace('_', ' ')) for tup in left_data]
        scoreOff = [tup[3] for tup in right_data]
        namesOff = [fixNames(tup[0].replace('_', ' ')) for tup in right_data]

        colms = ['on off', 'Score', 'Target name']
        offRows = pd.DataFrame(data={colms[0]: ['Off target']*len(scoreOff), colms[1]: scoreOff, colms[2]: namesOff})
        onRows = pd.DataFrame(data={colms[0]: ['On target']*len(scoreOn), colms[1]: scoreOn, colms[2]: namesOn})
        allRows = pd.concat([onRows,offRows], axis=0, join='outer')
        ax[1] = sns.violinplot(x='Target name', y='Score', hue='on off', palette=colr, data=allRows, split=splt,
                               inner=None, cut=0)
        ax[1].legend([], [], frameon=False)
    else:
        ax[1] = sns.violinplot(x=new_names, y=score, color=colr, inner=None, cut=0)
        # ax[1] = sns.swarmplot(x=newNames, y=score, color='white', edgecolor='gray', s=0.8)
    sns.despine(ax=ax[1], offset=5, bottom=True)

    plt.savefig(filename, transparent=True)
    return

def plot_var_reg(ax, var_regs):
    for V in var_regs:
        ax.axvspan(V[0], V[1], facecolor='g', alpha=0.2)
    return

