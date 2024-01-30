import matplotlib.pyplot as plt
import pickle
import numpy as np
from collections import defaultdict
import pandas as pd
import seaborn as sns
from alive_progress import alive_bar
import logomaker as lm
import random
import icecream as ic
import ribodesigner as rd
import scipy
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
import os


def import_data_to_df(designs_path, name, check_integrity: bool = True):
    lengths = [int(name.split('/')[-1].split('_')[0]) for name in designs_path]
    dfs = []
    all_dupes_check = defaultdict(lambda: [])
    all_without_dupes = defaultdict(lambda: [])

    for i in range(0, len(designs_path)):
        temp = extract_info(designs_path[i], name)
        if check_integrity:
            dupes_check = len(temp.index)
            all_dupes_check[lengths[i]].append(dupes_check)
            without_dupes = len(temp.drop_duplicates().index)
            all_without_dupes[lengths[i]].append(without_dupes)
            if dupes_check != without_dupes:
                print(f'\n{designs_path} has duplicate outputs!')
            dfs.append(temp)

    for key, val in all_without_dupes.items():
        all_vals = sum(val)
        if key == all_vals:
            print(f'\n{name} {key} designs file has correct number of outputs')
        else:
            print(f'\n{name} {key} designs file is not done cooking! Missing {key - all_vals} designs. '
                  f'Run again to complete analysis.')

    designs_df = pd.concat(dfs)
    return designs_df


def make_graphs(var_regs: list[tuple[int, int]], control_designs_path: list[str],
                universal_designs_path: list[str] = '', selective_designs_path: list[str] = '',
                ref_seq_designs_path: list[str] = '', random_seq_designs_path: list[str] = '', save_fig: bool = False,
                file_type: str = 'png', save_file_loc: str = ''):
    print('Loading data...')
    fn = lambda x: int(x.split('/')[-1].split('_')[0])
    percentages = []
    for paths in [control_designs_path, universal_designs_path, selective_designs_path, ref_seq_designs_path,
                  random_seq_designs_path]:
        if not paths:
            percentages.append(0)
            continue
        total_seqs = sum([fn(path) for path in paths])
        percentages.append(total_seqs)

    total_designs = sum(percentages)
    # percentages = [perc/total_designs for perc in percentages]
    with alive_bar(total=total_designs, spinner='fishes') as bar:
        control_designs_df = import_data_to_df(control_designs_path, 'control')
        bar(percentages[0])
        print('Control data loaded!')
        if universal_designs_path:
            print('Now loading universal data...')
            universal_designs_df = import_data_to_df(universal_designs_path, 'universal')
            bar(percentages[1])
            print('Universal data loaded!')
            all_data_df = pd.concat([control_designs_df, universal_designs_df])
        else:
            all_data_df = control_designs_df
        if selective_designs_path:
            print('Now loading selective data...')
            selective_designs_df = import_data_to_df(selective_designs_path, 'selective')
            bar(percentages[2])
            print('Selective data loaded!')
            all_data_df = pd.concat([all_data_df, selective_designs_df])
        if ref_seq_designs_path:
            print('Now loading reference data...')
            ref_seq_designs_df = import_data_to_df(ref_seq_designs_path, 'ref_seq')
            bar(percentages[3])
            print('Reference sequence data loaded!')
            all_data_df = pd.concat([all_data_df, ref_seq_designs_df])
        if random_seq_designs_path:
            print('Now loading random sequence data...')
            random_seq_designs_df = import_data_to_df(random_seq_designs_path, 'random_seq')
            bar(percentages[4])
            print('Reference sequence data loaded!')
            all_data_df = pd.concat([all_data_df, random_seq_designs_df])

    control_designs_df_1376 = control_designs_df[control_designs_df['reference_idx'] == 1376]
    # Get new column on all_data_df that has the target dataset
    labels = {'designs_e-coli': 'Reference', 'designs_Lacto': 'Random', 'designs_Bacteria': 'Bacteria',
              'designs_Archaea': 'Archaea', 'designs_Eukaryota': 'Eukaryota', 'designs_All': 'All kingdoms',
              'designs_TTCAC1376': 'Control'}
    new_target_labels = []
    for line in all_data_df['name_of_test_dataset']:
        txt = line.split('/')[-1].split('vs_test_sequences')[0]
        for key, name in labels.items():
            if key in txt:
                new_target_labels.append(name)
                break
    all_data_df.insert(0, 'target_name', new_target_labels)

    print('Generating graphs...')

    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 16 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # colors = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']  # viridis
    colors = ['#0D365E', '#3F6D54', '#9B892D', '#F9A281', '#FACDFB']  # batlow

    max_vals = all_data_df.max(numeric_only=True)

    # Fig 1: MG1655
    test_file_names = ['vs_test_sequences_Bac', 'vs_test_sequences_Arc', 'vs_test_sequences_Euk',
                       'vs_test_sequences_All']

    all_subsets = defaultdict(lambda: {})
    keys_for_subsets = {'control': ('control', None), 'reference': ('ref_seq', None), 'random': ('random_seq', None),
                        'universal bacterial': ('universal', '_designs_Bac'),
                        'universal archaeal': ('universal', 'designs_Arc'),
                        'universal eukaryotic': ('universal', 'designs_Euk'),
                        'universal all': ('universal', 'designs_All')}
    test_file_names = {'vs Bacteria': 'vs_test_sequences_Bac', 'vs Archaea': 'vs_test_sequences_Arc',
                       'vs Eukaryota': 'vs_test_sequences_Euk', 'vs all': 'vs_test_sequences_All'}
    # paired_names = [('reference', 'vs Bacteria'), ('random', 'vs Bacteria'),
    #                 ('universal bacterial', 'vs Bacteria'), ('universal archaeal', 'vs Archaea'),
    #                 ('universal eukaryotic', 'vs Eukaryota'), ('universal all', 'vs all')]
    paired_names = [('universal bacterial', 'vs Bacteria'), ('universal archaeal', 'vs Archaea'),
                    ('universal eukaryotic', 'vs Eukaryota'), ('universal all', 'vs all')]

    for subset_key, (design_type, target_name) in keys_for_subsets.items():
        for test_key, test_name in test_file_names.items():
            # First select out the design type
            target_subset = all_data_df.loc[all_data_df.index == design_type]
            # Next, get everyone in that design type that has the correct subset, if needed
            if target_name:
                target_subset = target_subset[target_subset['name_of_test_dataset'].str.contains(target_name)]
            # Finally, get everyone in correct design type and subset that has been tested against the correct dataset
            target_subset = target_subset[target_subset['name_of_test_dataset'].str.contains(test_name)]
            all_subsets[subset_key][test_key] = target_subset

    # Figures 2a and 2b
    bins_wanted = 100
    tm_max = max(all_data_df['tm_nn_vs_test'])
    tm_min = min(all_data_df['tm_nn_vs_test'])

    vars_to_plot = [('u_conservation_test', 'true_%_cov_test', 'figure_2a',
                     'U site conservation', 'IGS true percent coverage'),
                    ('test_score', 'tm_nn_vs_test', 'figure_2b', 'Guide score vs test', 'Tm GC of guide vs test')]

    # Get common values to later normalize
    min_max = {'u_conservation_test_vs_true_%_cov_test': [], 'test_score_vs_tm_nn_vs_test': []}
    for target_name, test_name in paired_names:
        for test_name in test_file_names.keys():
            target_subset = all_subsets[target_name][test_name]
            for xvar, yvar, _, _, _ in vars_to_plot:
                hist_data = np.histogram2d(target_subset[xvar], target_subset[yvar], bins=bins_wanted)
                min_max[f'{xvar}_vs_{yvar}'].append(hist_data[0])
    for key, item in min_max.items():
        vmin = min(subdata.min() for subdata in item)
        vmax = max(subdata.max() for subdata in item)
        min_max[key] = (vmin, vmax)

    # for target_name, test_name in paired_names:
    #     graph = all_subsets[target_name][test_name]
    #     top_control = all_subsets['control'][test_name]
    #     for xvar, yvar, figname, xlabel, ylabel in vars_to_plot:
    #         jointplot_fig = plt.figure()
    #         gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
    #         joint_ax = {
    #             0: jointplot_fig.add_subplot(gridspec[0:6, 7:14]),
    #             1: jointplot_fig.add_subplot(gridspec[0:3, 0:6]),
    #             2: jointplot_fig.add_subplot(gridspec[3:6, 0:6]),
    #         }
    #         # Plot scatter and kde plots
    #         # # below is for individual log norms
    #         # hist_data = jointplot_fig.axes[0].hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
    #         #                                          norm='log')
    #         hist_data = jointplot_fig.axes[0].hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
    #                                                  norm=LogNorm(vmin=1, vmax=min_max[f'{xvar}_vs_{yvar}'][-1]))
    #         jointplot_fig.colorbar(hist_data[-1], ax=jointplot_fig.axes[0])
    #         jointplot_fig.axes[0].scatter(x=top_control[xvar], marker='^', y=top_control[yvar], alpha=1, c='#000000',
    #                                       label='Control', edgecolors='black', linewidth=0.8,
    #                                       sizes=[150] * len(top_control[yvar]))
    #         regstats = scipy.stats.linregress(x=graph[xvar], y=graph[yvar])
    #         jointplot_fig.axes[0].annotate(f'$r^2$={round(regstats.rvalue ** 2, 3)}', xy=(0.1, 0.9),
    #                                        xycoords='axes fraction')
    #         jointplot_fig.axes[0].set_xlabel(xlabel)
    #         jointplot_fig.axes[0].set_ylabel(ylabel)
    #         # Set variable regions for location plots
    #         plot_variable_regions(joint_ax[1], var_regs)
    #         plot_variable_regions(joint_ax[2], var_regs)
    #         # Plot test data for each testing condition
    #         sns.scatterplot(x='reference_idx', y=xvar, data=graph, ax=joint_ax[1], alpha=0.3, legend=False, linewidth=0,
    #                         size=0.5)
    #         sns.scatterplot(x='reference_idx', y=yvar, data=graph, ax=joint_ax[2], alpha=0.3, legend=False, linewidth=0,
    #                         size=0.5)
    #         jointplot_fig.axes[2].set_xlabel('16s rRNA sequence position on reference sequence')
    #         # Plot control data
    #         jointplot_fig.axes[1].scatter(x=top_control['reference_idx'], y=top_control[xvar], alpha=1, c='#000000',
    #                                       label='control', marker='^', edgecolors='black', linewidth=0.8,
    #                                       sizes=[150] * len(top_control[yvar]))
    #         jointplot_fig.axes[1].set_ylabel(xlabel)
    #         jointplot_fig.axes[2].scatter(x=top_control['reference_idx'], y=top_control[yvar], alpha=1, c='#000000',
    #                                       label='control', marker='^', edgecolors='black', linewidth=0.8,
    #                                       sizes=[150] * len(top_control[yvar]))
    #         jointplot_fig.axes[2].set_ylabel(ylabel)
    #         if yvar == 'tm_nn_vs_test':
    #             jointplot_fig.axes[0].set(ylim=[tm_min, tm_max + 20], xlim=[0, 1])
    #         else:
    #             jointplot_fig.axes[0].set(ylim=[0, 1], xlim=[0, 1])
    #         jointplot_fig.axes[1].set(xlabel=None, ylim=[0, 1], xlim=[-0.1, max_vals['reference_idx'] + 20])
    #         jointplot_fig.axes[2].set(xlabel='Reference 16s rRNA index')
    #         jointplot_fig.axes[2].sharex(jointplot_fig.axes[1])
    #         jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
    #         jointplot_fig.axes[1].tick_params(labelbottom=False)
    #         title = f'Designs from {target_name} target sequences against {test_name} test dataset'
    #         jointplot_fig.suptitle(title)
    #         plt.tight_layout()
    #         save_file_name = f'{save_file_loc}/{figname}_{target_name}_{test_name}.{file_type}'
    # Figures 4, 5
    for target_name, test_name in paired_names:
        graph = all_subsets[target_name][test_name]
        top_control = all_subsets['control'][test_name]
        for xvar, yvar, figname, xlabel, ylabel in vars_to_plot:
            # Set plot parameters
            sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
            jointplot_fig = plt.figure()
            if xvar == 'u_conservation_test':
                gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=7)
                joint_ax = {
                    0: jointplot_fig.add_subplot(gridspec[0:6, 0:6]),
                }
                do_hist_plot(jointplot_fig, joint_ax[0], graph, xvar, yvar, bins_wanted, min_max, top_control, xlabel,
                             ylabel, tm_min, tm_max)
                jointplot_fig.axes[0].set(ylim=[0, 1], xlim=[0, 1])
                title = f'Designs from {target_name} target sequences against {test_name} test dataset'
                jointplot_fig.suptitle(title)
                plt.tight_layout()

                if save_fig:
                    save_file_name = f'{save_file_loc}/{figname}_{target_name}_{test_name}_a.{file_type}'
                    plt.savefig(fname=save_file_name, format=file_type)
                plt.show()
                save_file_name = f'{save_file_loc}/{figname}_{target_name}_{test_name}_b.{file_type}'
                plot_three_panel_graph(var_regs, graph['reference_idx'], graph[xvar], 0.3, graph[yvar], xlabel,
                                       ylabel, len(graph['reference_idx']), title, save_file_name, file_type)
            else:
                jointplot_fig = plt.figure()
                gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
                joint_ax = {
                    0: jointplot_fig.add_subplot(gridspec[0:6, 7:14]),
                    1: jointplot_fig.add_subplot(gridspec[0:3, 0:6]),
                    2: jointplot_fig.add_subplot(gridspec[3:6, 0:6]),
                }
                # Plot scatter and kde plots
                # # below is for individual log norms
                # hist_data = jointplot_fig.axes[0].hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
                #                                          norm='log')
                hist_data = jointplot_fig.axes[0].hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
                                                         norm=LogNorm(vmin=1, vmax=min_max[f'{xvar}_vs_{yvar}'][-1]))
                jointplot_fig.colorbar(hist_data[-1], ax=jointplot_fig.axes[0])
                jointplot_fig.axes[0].scatter(x=top_control[xvar], marker='^', y=top_control[yvar], alpha=1,
                                              c='#FCC2DC', label='Control', edgecolors='black', linewidth=0.8,
                                              sizes=[150] * len(top_control[yvar]))
                regstats = scipy.stats.linregress(x=graph[xvar], y=graph[yvar])
                jointplot_fig.axes[0].annotate(f'$r^2$={round(regstats.rvalue ** 2, 3)}', xy=(0.1, 0.9),
                                               xycoords='axes fraction')
                jointplot_fig.axes[0].set_xlabel(xlabel)
                jointplot_fig.axes[0].set_ylabel(ylabel)
                # Set variable regions for location plots
                plot_variable_regions(joint_ax[1], var_regs)
                plot_variable_regions(joint_ax[2], var_regs)
                # Plot test data for each testing condition
                sns.scatterplot(x='reference_idx', y=xvar, data=graph, ax=joint_ax[1], alpha=0.3, legend=False,
                                linewidth=0, size=0.5, color='#000000')
                sns.scatterplot(x='reference_idx', y=yvar, data=graph, ax=joint_ax[2], alpha=0.3, legend=False,
                                linewidth=0, size=0.5, color='#000000')
                jointplot_fig.axes[2].set_xlabel('16s rRNA sequence position on reference sequence')
                # Plot control data
                jointplot_fig.axes[1].scatter(x=top_control['reference_idx'], y=top_control[xvar], alpha=1, c='#FCC2DC',
                                              label='control', marker='^', edgecolors='black', linewidth=0.8,
                                              sizes=[150] * len(top_control[yvar]))
                jointplot_fig.axes[1].set_ylabel(xlabel)
                jointplot_fig.axes[2].scatter(x=top_control['reference_idx'], y=top_control[yvar], alpha=1, c='#FCC2DC',
                                              label='control', marker='^', edgecolors='black', linewidth=0.8,
                                              sizes=[150] * len(top_control[yvar]))
                jointplot_fig.axes[2].set_ylabel(ylabel)
                jointplot_fig.axes[0].set(ylim=[tm_min, tm_max + 20], xlim=[0, 1])
                jointplot_fig.axes[1].set(xlabel=None, ylim=[0, 1], xlim=[-0.1, max_vals['reference_idx'] + 20])
                jointplot_fig.axes[2].set(xlabel='Reference 16s rRNA index')
                jointplot_fig.axes[2].sharex(jointplot_fig.axes[1])
                jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
                jointplot_fig.axes[1].tick_params(labelbottom=False)
                title = f'Designs from {target_name} target sequences against {test_name} test dataset'
                jointplot_fig.suptitle(title)
                plt.tight_layout()

    # Fig 6
    # Set plot parameters
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=8, ncols=27)
    joint_ax = {
        0: jointplot_fig.add_subplot(gridspec[0:2, 0:6]),
        1: jointplot_fig.add_subplot(gridspec[0:2, 7:13]),
        2: jointplot_fig.add_subplot(gridspec[0:2, 14:20]),
        3: jointplot_fig.add_subplot(gridspec[0:2, 21:27]),

        4: jointplot_fig.add_subplot(gridspec[2:4, 0:6]),
        5: jointplot_fig.add_subplot(gridspec[2:4, 7:13]),
        6: jointplot_fig.add_subplot(gridspec[2:4, 14:20]),
        7: jointplot_fig.add_subplot(gridspec[2:4, 21:27]),

        8: jointplot_fig.add_subplot(gridspec[4:6, 0:6]),
        9: jointplot_fig.add_subplot(gridspec[4:6, 7:13]),
        10: jointplot_fig.add_subplot(gridspec[4:6, 14:20]),
        11: jointplot_fig.add_subplot(gridspec[4:6, 21:27]),

        12: jointplot_fig.add_subplot(gridspec[6:8, 0:6]),
        13: jointplot_fig.add_subplot(gridspec[6:8, 7:13]),
        14: jointplot_fig.add_subplot(gridspec[6:8, 14:20]),
        15: jointplot_fig.add_subplot(gridspec[6:8, 21:27]),
    }
    for i in range(0, len(joint_ax)):
        plot_variable_regions(joint_ax[i], var_regs)

    # Each row is a group, each column is a metric to graph
    vars_and_titles = [('U conservation', 'u_conservation_test', [0, 1]),
                       ('IGS conservation', 'true_%_cov_test', [0, 1]),
                       ('Guide score', 'test_score', [0, 1]),
                       ('Tm GC', 'tm_nn_vs_test', [tm_min, tm_max])]
    current_graph = 0
    for target_name, test_name in paired_names:
        graph = all_subsets[target_name][test_name]
        top_control = all_subsets['control'][test_name]

        # Plot test data for each testing condition
        for i, (title, yvar, lims) in enumerate(vars_and_titles):
            sns.scatterplot(x='reference_idx', y=yvar, linewidth=0, size=0.5, data=graph,
                            ax=joint_ax[current_graph + i], alpha=0.3, legend=False, color='#000000')
            jointplot_fig.axes[current_graph + i].scatter(x=top_control['reference_idx'], y=top_control[yvar],
                                                          alpha=1, c='#FCC2DC', edgecolors='black', linewidth=0.8,
                                                          marker='^', sizes=[100] * len(top_control[yvar]))
            jointplot_fig.axes[current_graph + i].set(xlabel=None, ylabel=None, ylim=lims,
                                                      xlim=[-0.1, max_vals['reference_idx'] + 20])
            jointplot_fig.axes[current_graph + i].set_title(f'{title} {target_name} {test_name}', loc='left')
            jointplot_fig.axes[current_graph + i].tick_params(labelbottom=False)
        current_graph += 4

    # Now make the axes conserved
    multipliers = [4, 8, 12]
    for i in range(4):
        for j in multipliers:
            jointplot_fig.axes[j].sharey(jointplot_fig.axes[i])
            jointplot_fig.axes[j].sharex(jointplot_fig.axes[i])
        jointplot_fig.axes[multipliers[-1]].set(xlabel='Reference 16s rRNA index')
        multipliers = [val + 1 for val in multipliers]
    plt.tight_layout()
    if save_fig:
        plt.savefig(fname=f'{save_file_loc}/figure_3a.{file_type}', format=file_type)
    plt.show()

    # Fig 7: : Assessing universal design quality. 2a) IGS true percent coverage vs. guide score of bacterial designs
    # evaluated against datasets of different kingdoms. 2b) 16s rRNA location of all designs along the reference
    # E. coli MG1655 16s rRNA sequence.
    # Prepare axes
    xvar = 'true_%_cov_test'
    yvar = 'test_score'
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (20 * 0.8, 20 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # Prepare histograms:
    all_vals = []
    for target_name, _ in paired_names:
        for test_name in test_file_names.keys():
            target_subset = all_subsets[target_name][test_name]
            data = np.histogram2d(target_subset[xvar], target_subset[yvar], bins=bins_wanted)
            all_vals.append(data[0])
    # vmin = min(subdata.min() for subdata in all_vals)
    vmax = max(subdata.max() for subdata in all_vals)
    # plot each test dataset in a different quadrant
    for target_name, target_data in paired_names:
        jointplot_fig = plt.figure()
        gridspec = jointplot_fig.add_gridspec(nrows=15, ncols=15)
        joint_ax = {
            0: jointplot_fig.add_subplot(gridspec[0:7, 0:7]),
            1: jointplot_fig.add_subplot(gridspec[0:7, 8:15]),
            2: jointplot_fig.add_subplot(gridspec[8:15, 0:7]),
            3: jointplot_fig.add_subplot(gridspec[8:15, 8:15])
        }
        for i, test_name in enumerate(test_file_names.keys()):
            # Plot a testing set that's relevant
            graph = all_subsets[target_name][test_name]
            jointplot_fig.axes[i].axhline(y=0.7, color='orange', linestyle='-')
            jointplot_fig.axes[i].axvline(x=0.7, color='orange', linestyle='-')
            hist_data = jointplot_fig.axes[i].hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
                                                     norm=LogNorm(vmin=1, vmax=vmax))
            # Plot control datapoints
            top_control = all_subsets['control'][test_name]
            jointplot_fig.axes[i].scatter(x=top_control[xvar], y=top_control[yvar], alpha=1, c='#FCC2DC',
                                          label='Control', marker='^', edgecolors='black', linewidth=0.8,
                                          sizes=[150] * len(top_control[yvar]))
            jointplot_fig.axes[i].set(xlabel=None, ylabel=None, xlim=[0, 1], ylim=[0, 1])
            jointplot_fig.axes[i].set_title(test_name, loc='left')
            num_df = graph[((graph[xvar] >= 0.7) & (graph[yvar] >= 0.7))]
            jointplot_fig.axes[i].annotate(f'n = {len(num_df.index)}\ngood guides', xy=(0.75, 0.15),
                                           xycoords='axes fraction')
            jointplot_fig.axes[i].label_outer()
        jointplot_fig.colorbar(hist_data[-1], ax=jointplot_fig.axes)
        jointplot_fig.suptitle(f'Designs from {target_name} evaluated against different kingdoms')
        jointplot_fig.axes[2].set_xlabel('IGS true percent coverage')
        jointplot_fig.axes[3].set_xlabel('IGS true percent coverage')
        jointplot_fig.axes[0].set_ylabel('Guide score')
        jointplot_fig.axes[2].set_ylabel('Guide score')
        plt.tight_layout()
        if save_fig:
            plt.savefig(fname=f'{save_file_loc}/figure_3b_{target_name}.{file_type}', format=file_type)
        plt.show()

    # # Fig 3b: : Assessing universal design quality. 2a) IGS true percent coverage vs. guide score of bacterial designs
    # # evaluated against datasets of different kingdoms. 2b) 16s rRNA location of all designs along the reference
    # # E. coli MG1655 16s rRNA sequence.
    # for target_name, inner_dict in all_subsets.items():
    #     if target_name == 'control':
    #         continue
    #     subset_to_graph = pd.DataFrame()
    #     for df in inner_dict.values():
    #         subset_to_graph = pd.concat([subset_to_graph, df])
    #
    #     # Prepare axes
    #     jointplot_fig = plt.figure()
    #     gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
    #     joint_ax = {
    #         0: jointplot_fig.add_subplot(gridspec[1:7, 7:13]),
    #         1: jointplot_fig.add_subplot(gridspec[0:1, 7:13]),
    #         2: jointplot_fig.add_subplot(gridspec[1:7, 13:14]),
    #         3: jointplot_fig.add_subplot(gridspec[0:2, 0:6]),
    #         4: jointplot_fig.add_subplot(gridspec[2:4, 0:6]),
    #         5: jointplot_fig.add_subplot(gridspec[4:6, 0:6])
    #     }
    #     xvar = 'true_%_cov_test'
    #     yvar = 'test_score'
    #     # Plot scatter and kde plots
    #     sns.scatterplot(x=xvar, y=yvar, hue='name_of_test_dataset', data=subset_to_graph, ax=joint_ax[0], alpha=0.2,
    #                     legend=False)
    #     # sns.kdeplot(x=xvar, y=yvar, hue='name_of_test_dataset', data=subset_to_graph, ax=joint_ax[0], fill=True,
    #     #             cmap='viridis', legend=False)
    #     markers = ['v', '^', '<', '>']
    #     for i, marker in enumerate(markers):
    #         jointplot_fig.axes[0].scatter(x=control_designs_df_1376[xvar].iloc[i],
    #                                       y=control_designs_df_1376[yvar].iloc[i],
    #                                       alpha=1, c='#fde725', marker=marker, label='Control')
    #
    #     # sns.scatterplot(x=xvar, y=yvar, style='name_of_test_dataset', data=control_designs_df_1376,
    #     #                 ax=jointplot_fig.axes[0], alpha=1, legend=False, palette=['#fde725'])
    #     # jointplot_fig.axes[0].scatter(x=control_designs_df_1376[xvar], y=control_designs_df_1376[yvar],
    #     #                               alpha=1, c='#fde725', marker=markers.pop(), label='Control')
    #     sns.kdeplot(x=xvar, hue='name_of_test_dataset', data=subset_to_graph, ax=joint_ax[1], fill=True,
    #                 common_norm=True, alpha=.3, legend=False)
    #     sns.kdeplot(y=yvar, hue='name_of_test_dataset', data=subset_to_graph, ax=joint_ax[2], fill=True,
    #                 common_norm=True, alpha=.3, legend=False)
    #     jointplot_fig.axes[0].set_xlabel('IGS true percent coverage')
    #     jointplot_fig.axes[0].set_ylabel('Guide score')
    #     # Set variable regions for location plots
    #     plot_variable_regions(joint_ax[3], var_regs)
    #     plot_variable_regions(joint_ax[4], var_regs)
    #     plot_variable_regions(joint_ax[5], var_regs)
    #     # Plot test data for each testing condition
    #     sns.scatterplot(x='reference_idx', y='u_conservation_test', hue='name_of_test_dataset', data=subset_to_graph,
    #                     ax=joint_ax[3], alpha=0.2, legend=False)
    #     sns.scatterplot(x='reference_idx', y=xvar, hue='name_of_test_dataset', data=subset_to_graph, ax=joint_ax[4],
    #                     alpha=0.2, legend=False)
    #     sns.scatterplot(x='reference_idx', y=yvar, hue='name_of_test_dataset', data=subset_to_graph, ax=joint_ax[5],
    #                     alpha=0.2)
    #     jointplot_fig.axes[4].set_xlabel('16s rRNA sequence position on reference sequence')
    #     # Plot control data
    #     for i, marker in enumerate(markers):
    #         jointplot_fig.axes[3].scatter(x=control_designs_df_1376['reference_idx'].iloc[i],
    #                                       y=control_designs_df_1376['u_conservation_test'].iloc[i],
    #                                       alpha=1, c='#fde725', marker=marker, label='Control')
    #
    #         jointplot_fig.axes[4].scatter(x=control_designs_df_1376['reference_idx'].iloc[i],
    #                                       y=control_designs_df_1376[xvar].iloc[i],
    #                                       alpha=1, c='#fde725', marker=marker, label='control')
    #         jointplot_fig.axes[5].scatter(x=control_designs_df_1376['reference_idx'].iloc[i],
    #                                       y=control_designs_df_1376[yvar].iloc[i], alpha=1,
    #                                       c='#fde725', marker=marker, label='control')
    #     # sns.scatterplot(x='reference_idx', y='u_conservation_test', style='name_of_test_dataset',
    #     #                 data=control_designs_df_1376, ax=joint_ax[3], palette=['#fde725'], alpha=1, legend=False)
    #     # sns.scatterplot(x='reference_idx', y=xvar, style='name_of_test_dataset',
    #     #                 data=control_designs_df_1376, ax=joint_ax[4], palette=['#fde725'], alpha=1, legend=False)
    #     # sns.scatterplot(x='reference_idx', y=yvar, style='name_of_test_dataset',
    #     #                 data=control_designs_df_1376, ax=joint_ax[5], palette=['#fde725'], alpha=1, legend=False)
    #
    #     # jointplot_fig.axes[3].scatter(x=top_score_control['reference_idx'], y=top_score_control['u_conservation_test'],
    #     #                               alpha=1, c='#fde725', label='control')
    #     jointplot_fig.axes[3].set_ylabel('U site conservation')
    #     # jointplot_fig.axes[4].scatter(x=top_score_control['reference_idx'], y=top_score_control[xvar],
    #     #                               alpha=1, c='#fde725', label='control')
    #     jointplot_fig.axes[4].set_ylabel('IGS true percent coverage')
    #     # jointplot_fig.axes[5].scatter(x=top_score_control['reference_idx'], y=top_score_control[yvar], alpha=1,
    #     #                               c='#fde725', label='control')
    #     jointplot_fig.axes[5].set_ylabel('Guide score')
    #     # Set graph settings for pretti graphing
    #     jointplot_fig.axes[0].set(xlim=[-0.1, 1.1], ylim=[-0.1, 1.1])
    #     jointplot_fig.axes[0].set(ylim=[-0.1, 1.1])
    #     jointplot_fig.axes[1].set(xlabel=None)
    #     # jointplot_fig.axes[1].set_title('D', loc='left', fontsize=30)
    #     jointplot_fig.axes[2].set(ylabel=None)
    #     # jointplot_fig.axes[3].set_title('A', loc='left', fontsize=30)
    #     jointplot_fig.axes[3].set(xlabel=None, xlim=[-0.1, max_vals['reference_idx'] + 20], ylim=[-0.1, 1.1])
    #     jointplot_fig.axes[4].sharex(jointplot_fig.axes[3])
    #     jointplot_fig.axes[4].sharey(jointplot_fig.axes[3])
    #     jointplot_fig.axes[5].sharex(jointplot_fig.axes[3])
    #     jointplot_fig.axes[5].sharey(jointplot_fig.axes[3])
    #     jointplot_fig.axes[4].set(xlabel=None)
    #     # jointplot_fig.axes[4].set_title('B', loc='left', fontsize=30)
    #     # jointplot_fig.axes[5].set_title('C', loc='left', fontsize=30)
    #     jointplot_fig.axes[5].set(xlabel='Reference 16s rRNA index')
    #     jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    #     jointplot_fig.axes[1].tick_params(labelbottom=False, labelleft=False, left=False)
    #     jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
    #     jointplot_fig.axes[2].tick_params(labelbottom=False, labelleft=False, bottom=False)
    #     jointplot_fig.axes[3].tick_params(labelbottom=False)
    #     jointplot_fig.axes[4].tick_params(labelbottom=False)
    #
    #     legend = plt.legend(bbox_to_anchor=(1.4, -0.15), ncols=5, title='Test dataset')
    #     for label in legend.get_texts():
    #         txt = label.get_text().split('/')[-1].split('vs_test_sequences')[-1].replace('_', ' ')
    #         for name in inner_dict.keys():
    #             if name[3:].casefold() in txt.casefold():
    #                 label.set_text(name)
    #                 break
    #
    #     jointplot_fig.suptitle(f'Designs from {target_name} target sequences against test datasets')
    #     plt.tight_layout()
    #
    #     if save_fig:
    #         plt.savefig(fname=f'{save_file_loc}/figure_2a_{target_name}.{file_type}', format=file_type)
    #     plt.show()

    # # Graph 1a:
    # # Find max y label:
    # y_max = all_data_df.max(numeric_only=True)
    # fig1a, ax1a = plt.subplots(nrows=1, ncols=2, layout='constrained',
    #                            sharex='all')
    # fig1a.suptitle('Number of targets vs. delta composite scores')
    # sns.scatterplot(x='delta_vs_test', y='num_of_targets_test', data=universal_designs_df,
    #                 ax=ax1a[0], alpha=0.7, hue=universal_designs_df.index.name)
    # ax1a[0].set_title('Universal vs. Test Dataset')
    # ax1a[0].set_ylim(bottom=0, top=max(y_max['num_of_targets_test'], y_max['num_of_targets_background']))
    # ax1a[0].set_xlabel('Delta composite score vs. test dataset')
    # ax1a[0].set_ylabel('Number of targets in test dataset')
    # sns.scatterplot(x='delta_vs_background', y='num_of_targets_background', data=selective_designs_df,
    #                 ax=ax1a[1], alpha=0.7, hue=selective_designs_df.index.name)
    # ax1a[1].set_title('Selective vs. Background Dataset')
    # ax1a[1].set_ylim(bottom=0, top=max(y_max['num_of_targets_test'], y_max['num_of_targets_background']))
    # ax1a[1].set_xlabel('Delta composite score vs. backpround dataset')
    # ax1a[1].set_ylabel('Number of targets in backpround dataset')
    #
    # plt.tight_layout()
    # if save_fig:
    #     plt.savefig(fname=f'{save_file_loc}/fig1a.{file_type}', format=file_type)
    # plt.show()
    #
    # def make_jointplot(x_var, y_var, name):
    #     jointplot_fig = plt.figure()
    #     gridspec = jointplot_fig.add_gridspec(nrows=4, ncols=9)
    #     joint_ax = {
    #         0: jointplot_fig.add_subplot(gridspec[1:4, 0:3]),
    #         1: jointplot_fig.add_subplot(gridspec[0:1, 0:3]),
    #         2: jointplot_fig.add_subplot(gridspec[1:4, 3:4]),
    #         3: jointplot_fig.add_subplot(gridspec[1:4, 5:8]),
    #         4: jointplot_fig.add_subplot(gridspec[0:1, 5:8]),
    #         5: jointplot_fig.add_subplot(gridspec[1:4, 8:9])
    #     }
    #     for i, dset, labels in zip([0, 3], [universal_designs_df, selective_designs_df],
    #                                [universal_labels, selective_labels]):
    #         slope, intercept, r, p, sterr = scipy.stats.linregress(x=dset[x_var],
    #                                                                y=dset[y_var])
    #         sns.scatterplot(x=x_var, y=y_var, hue=dset.index.name, data=dset, ax=joint_ax[i],
    #                         alpha=0.7)
    #         sns.kdeplot(x=x_var, hue=dset.index.name, data=dset, ax=joint_ax[i + 1],
    #                     fill=True, common_norm=True, alpha=.3, legend=False)
    #         sns.kdeplot(y=y_var, hue=dset.index.name, data=dset, ax=joint_ax[i + 2],
    #                     fill=True, common_norm=True, alpha=.3, legend=False)
    #
    #         joint_ax[i].annotate(f'$r^2$={round(r, 3)}', xy=(0.1, 0.9), xycoords='axes fraction')
    #
    #     jointplot_fig.axes[0].set(xlim=[-0.1, 1.1], ylim=[-0.1, 1.1])
    #     jointplot_fig.axes[3].sharex(jointplot_fig.axes[0])
    #     jointplot_fig.axes[3].sharey(jointplot_fig.axes[0])
    #     jointplot_fig.axes[1].set(xlabel=None)
    #     jointplot_fig.axes[2].set(ylabel=None)
    #     jointplot_fig.axes[4].set(xlabel=None)
    #     jointplot_fig.axes[5].set(ylabel=None)
    #     jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    #     jointplot_fig.axes[1].tick_params(labelbottom=False)
    #     jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
    #     jointplot_fig.axes[2].tick_params(labelleft=False)
    #     jointplot_fig.axes[4].sharex(jointplot_fig.axes[3])
    #     jointplot_fig.axes[4].tick_params(labelbottom=False)
    #     jointplot_fig.axes[5].sharey(jointplot_fig.axes[3])
    #     jointplot_fig.axes[5].tick_params(labelleft=False)
    #
    #     if save_fig:
    #         plt.savefig(fname=f'{save_file_loc}/{name}.{file_type}', format=file_type)
    #     plt.show()
    #     return
    #
    # make_jointplot('true_%_cov_test', 'test_score', name='/fig1b')
    # make_jointplot('delta_igs_vs_test', 'delta_guide_vs_test', name='/fig1c')
    #
    # # Graph 2: y-axis is the composite score, x-axis is the 16s rRNA gene, plot the universal, control, and each
    # # selective for all designs in different panels (same data as above but order along gene)
    # fig2a, ax2a = plt.subplots(nrows=max(num_universal, num_selective), ncols=2, layout='constrained', sharey='all',
    #                            sharex='all')
    # fig2a.suptitle('Scores along 16s rRNA sequences in test datasets: universal/ selective')
    # fig_2a_titles = ['Average conservation of catalytic site', '% of reads with IGS', 'Guide score']
    # fig_2a_columns = ['u_conservation_test', 'true_%_cov_test', 'test_score']
    # ax2a[0, 0].set_title('Universal designs')
    # for i, name, col in zip(range(3), fig_2a_titles, fig_2a_columns):
    #     plot_variable_regions(ax2a[i, 0], var_regs)
    #     ax2a[i, 0].scatter(x=top_score_control['reference_idx'], y=top_score_control[col], alpha=0.7, c='#fde725',
    #                        label='control')
    #     for j, dataset in enumerate(universal_labels):
    #         ax2a[i, 0].scatter(x=universal_designs_df.loc[dataset, 'reference_idx'],
    #                            y=universal_designs_df.loc[dataset, col], alpha=0.5, c=colors[j], label=dataset)
    #     ax2a[i, 0].set_ylabel(name)
    #     ax2a[i, 0].legend()
    #     ax2a[i, 0].set_xlabel('Position in reference sequence')
    # ax2a[0, 1].set_title('Selective designs')
    # for i, name, col in zip(range(3), fig_2a_titles, fig_2a_columns):
    #     plot_variable_regions(ax2a[i, 1], var_regs)
    #     ax2a[i, 1].scatter(x=top_score_control['reference_idx'], y=top_score_control[col], alpha=0.7, c='#fde725',
    #                        label='control')
    #     for j, dataset in enumerate(selective_labels):
    #         ax2a[i, 1].scatter(x=selective_designs_df.loc[dataset, 'reference_idx'],
    #                            y=selective_designs_df.loc[dataset, col], alpha=0.5, c=colors[j], label=dataset)
    #     ax2a[i, 1].legend()
    # ax2a[i, 0].set_xlabel('Position in reference sequence')
    # plt.tight_layout()
    # if save_fig:
    #     plt.savefig(fname=f'{save_file_loc}/fig2a.{file_type}', format=file_type)
    # plt.show()
    #
    # fig2b, ax2b = plt.subplots(nrows=2, ncols=1, sharex='all')
    # fig2b.suptitle('Composite score along test 16s sequences')
    #
    # for i, (dset, name) in enumerate(zip([universal_designs_df, selective_designs_df], ['Universal', 'Selective'])):
    #     plot_variable_regions(ax2b[i], var_regs)
    #     ax2b[i].scatter(x=top_score_control['reference_idx'], y=top_score_control['composite_test_score'],
    #                     alpha=0.7, c='#fde725', label='control')
    #     for j, dataset in enumerate(list(set(dset.index))):
    #         ax2b[i].scatter(x=dset.loc[dataset, 'reference_idx'],
    #                         y=dset.loc[dataset, 'composite_test_score'], alpha=0.5, c=colors[j], label=dataset)
    #     ax2b[i].set_title(name)
    #     ax2b[i].set_ylabel('Composite test score')
    # ax2b[1].set_xlabel('16s rRNA sequence position on reference sequence')
    # plt.tight_layout()
    # if save_fig:
    #     plt.savefig(fname=f'{save_file_loc}/fig2b.{file_type}', format=file_type)
    # plt.show()
    #
    # # Graph 3: violin plot, showing composite score distribution in all designs of a given category
    # back_score_data = all_data_df.loc[:, ['id', 'test_score']]
    # back_score_data['score_type'] = 'Guide score on test data'
    # back_score_data.rename(columns={'test_score': 'Score'}, inplace=True)
    # back_true_cov_data = all_data_df.loc[:, ['id', 'true_%_cov_test']]
    # back_true_cov_data['score_type'] = 'True % coverage on test data'
    # back_true_cov_data.rename(columns={'true_%_cov_test': 'Score'}, inplace=True)
    # fig_3_data = pd.concat([back_score_data, back_true_cov_data])
    # plt.title('Guide and IGS scores distribution on test dataset')
    # # sns.violinplot(x=fig_3_data.index, y='score', hue='score_type', data=fig_3_data, inner='point', cut=0, split=True)
    # sns.boxplot(x=fig_3_data.index, y='Score', hue='score_type', data=fig_3_data, notch=True, boxprops={'alpha': 0.7})
    # sns.stripplot(x=fig_3_data.index, y='Score', hue='score_type', data=fig_3_data, dodge=True)
    #
    # plt.tight_layout()
    # if save_fig:
    #     plt.savefig(fname=f'{save_file_loc}/fig3.{file_type}', format=file_type)
    # plt.show()
    #
    # # Graph 4: y-axis is the delta score, x-axis is the 16s rRNA gene, plot all selective designs in different panels
    # fig4, ax4 = plt.subplots(nrows=2, ncols=1, sharex='all')
    # for i, dset in enumerate([universal_designs_df, selective_designs_df]):
    #     plot_variable_regions(ax4[i], var_regs)
    #     ax4[i].scatter(x=top_score_control['reference_idx'], y=top_score_control['delta_vs_test'],
    #                    alpha=0.7, c='#fde725', label='control')
    #     for j, dataset in enumerate(list(set(dset.index))):
    #         ax4[i].scatter(x=dset.loc[dataset, 'reference_idx'], y=dset.loc[dataset, 'delta_vs_test'],
    #                        alpha=0.5, c=colors[j], label=dataset)
    #     ax4[i].set_ylabel('Delta composite score vs. test dataset')
    # ax4[0].set_title('Universal')
    # ax4[1].set_title('Selective')
    # ax4[1].set_xlabel('16s rRNA sequence position on reference sequence')
    #
    # if save_fig:
    #     plt.savefig(fname=f'{save_file_loc}/fig4.{file_type}', format=file_type)
    # plt.show()
    #
    # # Extract taxonomy levels from test dataset too
    # for i, test_folder in enumerate(test_folders):
    #     test_names_and_seqs = read_silva_fasta(in_file=test_folder)
    #     fn_align_to_ref = np.vectorize(give_taxonomy, otypes=[TargetSeq], excluded=['level'])
    #     test_dataset_taxonomies = fn_align_to_ref(test_names_and_seqs, level=taxonomy)
    #
    #     test_counts_and_names = dict(np.asarray(np.unique(test_dataset_taxonomies, return_counts=True)).T)
    #     # Extract possible taxonomy labels
    #     test_target_col_name = f'target_{taxonomy}_test'
    #     all_targets_set = set()
    #     all_targets_raw = dict(top_test_scores[test_target_col_name])
    #
    #     all_targets_dict = {f'Test data {i}': test_counts_and_names}
    #     for key, val in all_targets_raw.items():
    #         val = val.replace(';', ',')
    #         all_targets_dict[key] = eval(val)
    #         all_targets_set.update(set(all_targets_dict[key].keys()))
    #
    # fig_5_data = pd.DataFrame(all_targets_dict).T
    #
    # filter_nans = fig_5_data.groupby(fig_5_data.index).sum()
    #
    # fig5, ax5 = plt.subplots(ncols=2, nrows=1, sharey='all', layout='constrained')
    # # Figure 5a:
    # # Get 100%
    # totals = filter_nans.sum(axis=1)
    # # Calculate fractions
    # fractions = filter_nans.div(totals, axis=0)
    # fractions.plot.barh(stacked=True, ax=ax5[0], legend=0)
    # ax5[0].set_xlabel('Fraction')
    # # Figure 5b
    # filter_nans.plot.barh(stacked=True, ax=ax5[1], legend=0)
    # ax5[1].set_xlabel('Counts')
    # fig5.suptitle(f'{taxonomy} targeted of dataset')
    # leg_handles, leg_labels = ax5[0].get_legend_handles_labels()
    # fig5.legend(leg_handles, leg_labels, ncols=math.ceil(len(leg_labels) / 14), loc='lower center', fancybox=True,
    #             fontsize='x-small')
    # # show the graph
    # plt.tight_layout()
    # plt.subplots_adjust(bottom=0.27)
    # if save_fig:
    #     plt.savefig(fname=f'{save_file_loc}/fig5.{file_type}', format=file_type)
    # plt.show()
    print('All graphs generated!')
    return


def make_test_seqs_graph(title: str, x_data: list, xlabel: str, y_data: list, ylabel: str,
                         loc_data: list[int], var_regs: list, save_file_name: str, alpha=0.5,
                         file_type: str = 'png', dataset_len: int = None):
    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 18 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # Prepare axes for first figure
    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=7, ncols=14)
    joint_ax = {
        0: jointplot_fig.add_subplot(gridspec[1:7, 0:7]),
        1: jointplot_fig.add_subplot(gridspec[0:1, 0:7]),
        2: jointplot_fig.add_subplot(gridspec[1:7, 7:14]),
        3: jointplot_fig.add_subplot(gridspec[0:1, 7:14]),
    }
    # Prepare hue data: color based on yes no var regs
    yes_no_var_reg = ['Conserved region'] * len(loc_data)
    for i, loc in enumerate(loc_data):
        for range_min, range_max in var_regs:
            if range_min <= loc <= range_max:
                yes_no_var_reg[i] = 'Variable region'
                break
    # Plot scatter and kde plots
    # color based on location:
    yes_no_palette = {'Conserved region': '#000000', 'Variable region': 'g'}
    sns.scatterplot(x=x_data, y=y_data, linewidth=0, alpha=alpha, hue=loc_data, ax=joint_ax[0], palette='viridis')
    norm = plt.Normalize(min(loc_data), max(loc_data))
    colorbar_data = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    colorbar_data.set_array([])
    joint_ax[0].get_legend().remove()
    jointplot_fig.colorbar(colorbar_data, ax=joint_ax[0], pad=0.1, orientation='horizontal',
                           label='Location along E. coli 16s')
    sns.scatterplot(x=x_data, y=y_data, linewidth=0, alpha=alpha, hue=yes_no_var_reg, ax=joint_ax[2],
                    palette=yes_no_palette)
    cmap = ListedColormap(colors=['#000000', 'g'])
    bounds = [0, 0.5, 1]
    norm = BoundaryNorm(bounds, cmap.N)
    colorbar_data = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    colorbar_data.set_array([])
    joint_ax[2].get_legend().remove()

    formatter = plt.FixedFormatter(['Conserved region', 'Variable region'])
    plt.colorbar(colorbar_data, ticks=plt.FixedLocator([0.25, 0.75]), format=formatter, ax=joint_ax[2], pad=0.1,
                 orientation='horizontal', label='Variable region location')

    sns.kdeplot(x=x_data, ax=joint_ax[1], fill=True, common_norm=True, alpha=.3, legend=False, color='#000000')
    # sns.kdeplot(y=y_data, ax=joint_ax[2], fill=True, common_norm=True, alpha=.3, legend=False, color='#000000')
    jointplot_fig.axes[1].tick_params(axis='both', labelleft=False)
    sns.kdeplot(x=x_data, ax=joint_ax[3], hue=yes_no_var_reg, fill=True, common_norm=True, alpha=.3, legend=False,
                palette=yes_no_palette)
    # sns.kdeplot(y=y_data, ax=joint_ax[5], hue=yes_no_var_reg, fill=True, common_norm=True, alpha=.3, legend=False,
    #             palette=yes_no_palette)
    jointplot_fig.axes[0].set(xlabel=xlabel, ylabel=ylabel, xlim=[0, 1], ylim=[0, 1])
    jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[1].tick_params(axis='both', labelleft=False, labelbottom=False, left=False)
    jointplot_fig.axes[2].set(xlabel=xlabel, xlim=[0, 1], ylim=[0, 1])
    jointplot_fig.axes[3].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[3].tick_params(axis='both', labelleft=False, labelbottom=False, left=False)
    jointplot_fig.axes[2].tick_params(axis='y', labelleft=False)
    jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} distributions '
                                                     f'n = {dataset_len}')
    plt.tight_layout()
    plt.savefig(fname=f'{save_file_name}_var_regs.{file_type}', format=file_type)
    plt.show()

    # Prepare axes for second figure
    plot_three_panel_graph(var_regs, loc_data, x_data, alpha, y_data, xlabel, ylabel, dataset_len, title,
                           save_file_name, file_type)
    return


def plot_three_panel_graph(var_regs, loc_data, x_data, alpha, y_data, xlabel, ylabel, dataset_len, title,
                           save_file_name, file_type):
    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 18 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # Prepare axes for figure
    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=9, ncols=7)
    joint_ax = {
        0: jointplot_fig.add_subplot(gridspec[0:3, 0:7]),
        1: jointplot_fig.add_subplot(gridspec[3:6, 0:7]),
        2: jointplot_fig.add_subplot(gridspec[6:9, 0:7])
    }
    # Set variable regions for location plots
    plot_variable_regions(joint_ax[0], var_regs)
    plot_variable_regions(joint_ax[1], var_regs)
    plot_variable_regions(joint_ax[2], var_regs)
    # Plot test data for each testing condition
    unique_locs = {loc: [set(), 0] for loc in set(loc_data)}
    for loc, xval in zip(loc_data, x_data):
        unique_locs[loc][0].add(xval)
        unique_locs[loc][1] += 1
    to_plot_x = []
    to_plot_loc = []
    to_plot_nums_of_vals = []
    for loc, (val, num_of_vals) in unique_locs.items():
        to_plot_x.append(*val)
        to_plot_loc.append(loc)
        to_plot_nums_of_vals.append(num_of_vals)
    if len(to_plot_loc) != len(to_plot_x):
        print('x data contains non-unique values and cannot be plotted by reducing duplicates')
        sns.scatterplot(x=loc_data, y=x_data, linewidth=0, size=0.5, color='#000000',
                        ax=jointplot_fig.axes[0], alpha=alpha / 2, legend=False)
    else:
        sns.scatterplot(x=to_plot_loc, y=to_plot_x, linewidth=0, size=0.5, color='#000000',
                        ax=jointplot_fig.axes[0], alpha=1, legend=False)
        jointplot_fig.axes[0].set_title(f'{len(to_plot_loc)} unique U sites', loc='left')
        jointplot_fig.axes[0].label_outer()
    jointplot_fig.axes[1].bar(to_plot_loc, to_plot_nums_of_vals, color='#000000', edgecolor='#000000')
    # jointplot_fig.axes[1].set_title(f'{len(loc_data)} unique U-IGS sites', loc='left')
    jointplot_fig.axes[1].label_outer()
    sns.scatterplot(x=loc_data, y=y_data, linewidth=0, size=0.5, color='#000000',
                    ax=jointplot_fig.axes[2], alpha=alpha / 2, legend=False)
    jointplot_fig.axes[2].set_title(f'{len(loc_data)} unique U-IGS sites', loc='left')
    # Set graph settings for pretti graphing
    jointplot_fig.axes[0].set_ylabel(xlabel)
    jointplot_fig.axes[1].set_ylabel('Number of U-IGS')
    jointplot_fig.axes[2].set_ylabel(ylabel)
    max_loc = max(loc_data)
    jointplot_fig.axes[0].set(xlabel=None, xlim=[-0.1, max_loc + 20], ylim=[-0.1, 1.1])
    jointplot_fig.axes[2].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
    jointplot_fig.axes[2].set(xlabel=None)
    jointplot_fig.axes[2].set(xlabel='Reference 16s rRNA index')
    jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[0].tick_params(labelbottom=False)
    jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} distributions '
                                                     f'n = {dataset_len}')
    plt.tight_layout()
    plt.savefig(fname=save_file_name + '.' + file_type, format=file_type)
    plt.show()
    return


def plot_variable_regions(ax, var_regs):
    for V in var_regs:
        ax.axvspan(V[0], V[1], facecolor='g', alpha=0.2)
    return


def extract_info(results_file_path: str, dataset: str):
    with open(results_file_path, 'r') as read_file:
        list_of_designs = read_file.readlines()
    # Now extact data here:
    column_types = {'id': str, 'igs': str, 'reference_idx': int, 'optimized_to_targets': bool,
                    'optimized_to_background': bool, 'tested': bool, 'tested_design': bool, 'guide': str,
                    'num_of_targets': int, 'score_type': str, 'score': float, '%_coverage': float, '%_on target': float,
                    'true_%_cov': float, 'composite_score': float, 'num_of_targets_background': int,
                    'u_conservation_background': float, 'background_score': float, '%_coverage_background': float,
                    '%_on target_background': float, 'true_%_cov_background': float,
                    'composite_background_score': float,
                    'delta_igs_vs_background': float, 'delta_guide_vs_background': float, 'delta_vs_background': float,
                    'name_of_test_dataset': str, 'num_of_targets_test': int, 'u_conservation_test': float,
                    'test_score': float, 'tm_nn_vs_test': float, '%_coverage_test': float, '%_on target_test': float,
                    'true_%_cov_test': float, 'composite_test_score': float, 'delta_igs_vs_test': float,
                    'delta_guide_vs_test': float, 'delta_vs_test': float, 'target_Domain': str,
                    'target_Domain_background': str, 'target_Domain_test': str, 'target_Phylum': str,
                    'target_Phylum_background': str, 'target_Phylum_test': str, 'target_Class': str,
                    'target_Class_background': str, 'target_Class_test': str, 'target_Order': str,
                    'target_Order_background': str, 'target_Order_test': str, 'target_Family': str,
                    'target_Family_background': str, 'target_Family_test': str, 'target_Genus': str,
                    'target_Genus_background': str, 'target_Genus_test': str, 'target_Species': str,
                    'target_Species_background': str, 'target_Species_test': str, 'target_Taxon': str,
                    'target_Taxon_background': str, 'target_Taxon_test': str}

    designs = []
    for design in list_of_designs:
        design_dict = {}
        individual_attributes = design.strip(' {').split(', "')
        for key_val in individual_attributes:
            key, val = key_val.split('": ')
            key = key.strip('"')
            try:
                val_typed = column_types[key](val.strip('"\''))
            except ValueError:
                if column_types[key] == int or column_types[key] == float:
                    val_typed = column_types[key](0)
                else:
                    print(column_types[key])
            design_dict[key] = val_typed
        # extracted_design = pd.DataFrame(design_dict, index=[dataset])
        # design_df = pd.concat([design_df, extracted_design])
        designs.append(design_dict)
    design_df = pd.DataFrame.from_records(designs, index=[dataset] * len(designs))
    return design_df


def make_sequence_logo_graph(test_data_path: str, design_data_path: list[str], ref_data_path: list[str],
                             save_fig: bool = False,
                             file_type: str = 'png', save_file_loc: str = '', guide_len: int = 50, igs_len: int = 5):
    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 16 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')

    # Extract data
    with open(test_data_path, 'rb') as handle:
        test_seqs = pickle.load(handle)

    design_data_df = import_data_to_df(design_data_path, 'universal')
    ref_data_df = import_data_to_df(ref_data_path, 'ref_data')

    # The test data is formatted for analysis so it will be a little clunky to extract necessary data here
    # Extract all positions all that one of the requirements below
    # 'high, high' = at least 0.9 in each
    # 'low, high' = at least 0.9 in u, at most 0.1 in igs
    # 'low, low' = at most 0.1 in each
    # 'mid, mid' = between 0.45 and 0.55 in each
    key_names = ['high, high', 'high, low', 'low, low', 'mid, mid']

    # matching_test_ref_idx = {key: {} for key in key_names}
    # igs_data = {}
    # igs_ids_in_pos = defaultdict(lambda: set())

    def truth_table(to_test_u, to_test_igs):
        truth_dict = {
            'high, high': lambda u, igs: (u >= 0.9) and (igs >= 0.9),
            'high, low': lambda u, igs: (u >= 0.9) and (igs <= 0.1),
            'low, low': lambda u, igs: (u <= 0.1) and (igs <= 0.1),
            'mid, mid': lambda u, igs: (0.40 <= u <= 0.60) and (0.40 <= igs <= 0.60),
        }
        for i in truth_dict.keys():
            if truth_dict[i](to_test_u, to_test_igs):
                return i
        return

    test_data_dicts = []
    for ref_idx, u_vals in test_seqs.items():
        u_conservation_test = u_vals[1]

        # Recall each inner dictionary is of form IGS_id: [guides, targets, perc cov, perc on target, true perc cov]
        for design_id, igs_vals in u_vals[0].items():
            igs = design_id[0:5]
            igs_true_perc_cov = igs_vals[4]
            condition = truth_table(u_conservation_test, igs_true_perc_cov)
            if not condition:
                continue
            # if igs.contains('N'):
            #     continue
            test_data_dicts.extend([{'condition': condition, 'id': design_id, 'reference_idx': ref_idx, 'igs': igs,
                                     'guide': str(single_guide), 'u_conservation': u_conservation_test,
                                     'igs_true_perc_cov': igs_true_perc_cov} for single_guide in igs_vals[0]])

            # matching_test_ref_idx[condition][ref_idx] = (u_conservation_test, igs_vals[0])
            # igs_ids_in_pos[ref_idx].update([design_id])
            # igs_data[design_id] = igs_true_perc_cov
            # matching_test_igs_id[condition][design_id] = (u_conservation_test, igs_true_perc_cov, igs_vals[0])
    test_data_df = pd.DataFrame.from_records(test_data_dicts, index=['test seqs'] * len(test_data_dicts))

    # Find designs and ref_sequence guides at those positions
    data_to_graph_u = {key: [] for key in key_names}
    data_to_graph_igs = {key: [] for key in key_names}
    u_sites_not_in_both = []
    igs_sites_not_in_both = []  # This will make it a bit easier to later extract information
    for key in key_names:
        vals = test_data_df[test_data_df['condition'] == key]
        matching_designs = design_data_df[design_data_df['reference_idx'].isin(vals['reference_idx'])]
        matching_ref_seqs = ref_data_df[ref_data_df['reference_idx'].isin(vals['reference_idx'])]
        # doing a set because pandas is being temperamental with inner joins and so am I
        matching_idxes = set(matching_designs['reference_idx']) & set(matching_ref_seqs['reference_idx'])

        if not matching_idxes:
            if not set(matching_designs['reference_idx']) and set(matching_ref_seqs['reference_idx']):
                u_sites_not_in_both.append(('only_in_ref', [set(matching_ref_seqs['reference_idx'])]))
            elif set(matching_designs['reference_idx']) and not set(matching_ref_seqs['reference_idx']):
                u_sites_not_in_both.append(('only_in_design', [set(matching_designs['reference_idx'])]))
            else:  # should be empty
                u_sites_not_in_both.append(('neither', [set(matching_designs['reference_idx']),
                                                        set(matching_ref_seqs['reference_idx'])]))
            continue
        # Ok now extract the data for only the indexes that match:
        ref_data = matching_ref_seqs[matching_ref_seqs['reference_idx'].isin(matching_idxes)][
            ['id', 'reference_idx', 'igs', 'guide', 'composite_test_score']]
        design_data = matching_designs[matching_designs['reference_idx'].isin(matching_idxes)][
            ['id', 'reference_idx', 'igs', 'guide', 'composite_test_score']]
        matching_igses_test = set([igs for idx in matching_idxes for igs in vals[vals['reference_idx'] == idx]['id']])
        matching_igses = set(design_data['id']) & set(ref_data['id']) & matching_igses_test

        # We only want to graph igses that we can analyze is both designs and ref_seqs
        if not matching_igses:
            if not set(design_data['igs']) and set(ref_data['igs']):
                igs_sites_not_in_both.append(('only_in_ref', [set(ref_data['igs'])]))
            elif set(design_data['igs']) and not set(ref_data['igs']):
                igs_sites_not_in_both.append(('only_in_design', [set(design_data['igs'])]))
            else:  # should be empty
                igs_sites_not_in_both.append(('neither', [set(ref_data['igs']), set(design_data['igs'])]))
            continue

        idx_with_igs = set(int(igs_id[5:]) for igs_id in matching_igses)
        with alive_bar(len(idx_with_igs), spinner='fishes') as bar:
            for idx in idx_with_igs:
                # Extract data to graph for each dataset
                # separate into two steps so we can use the dataframe in the next loop
                ref_data_temp = ref_data[ref_data['reference_idx'] == idx][['id', 'igs', 'guide']]
                ref_u_counted = lm.alignment_to_matrix(ref_data_temp['igs'])
                design_data_temp = design_data[design_data['reference_idx'] == idx][
                    ['id', 'igs', 'guide', 'composite_test_score']]
                design_u_counted = lm.alignment_to_matrix(design_data_temp['igs'])
                test_data_temp = test_data_df[test_data_df['reference_idx'] == idx][['id', 'igs', 'guide']]
                test_u_counted = lm.alignment_to_matrix(test_data_temp['igs'])
                best_score_design = design_data_temp.nlargest(1, 'composite_test_score')
                best_u_counted = lm.alignment_to_matrix(best_score_design['igs'])

                # # get consensus sequence for test dataset
                # test_consensus = Bio.motifs.create(test_data_temp['igs'],
                #                                    alphabet='GATCRYWSMKHBVDN-').degenerate_consensus.strip('-')
                mat = rd.words2countmatrix(test_data_temp['igs'])
                test_consensus = pairwise_comparison_consensus = rd.consensus(mat)
                # And now get score for best U and ref U with the consensus
                design_vs_test_u = rd.pairwise_comparison(seq_a=test_consensus,
                                                          seq_b=best_score_design['igs'].values[0], opti_len=igs_len,
                                                          score_type='weighted', only_consensus=True)
                ref_vs_test_u = rd.pairwise_comparison(seq_a=test_consensus, seq_b=ref_data_temp['igs'].values[0],
                                                       score_type='weighted', only_consensus=True, opti_len=igs_len)

                # Save in a tuple with:
                # (idx, [extracted_data_test, best_design, extracted_data_design, extracted_data_ref], u_cons)
                data_to_graph_u[key].append((idx, [test_u_counted, best_u_counted, design_u_counted, ref_u_counted],
                                             vals['u_conservation'].iloc[0], design_vs_test_u, ref_vs_test_u))
                to_scan = set([igs for igs in vals[vals['reference_idx'] == idx]['id']])
                for igs_id in to_scan:
                    if igs_id not in matching_igses:
                        continue
                    # Extract data to graph for each dataset
                    igs_true_perc_cov = vals[vals['id'] == igs_id]['igs_true_perc_cov']
                    if len(igs_true_perc_cov.unique()) != 1:
                        print('uh oh')
                    igs_probs = lm.alignment_to_matrix(['G' + igs_id[0:5]])

                    design_data_igs_id = design_data[design_data['id'] == igs_id]
                    ref_data_igs_id = ref_data_temp[ref_data_temp['id'] == igs_id]
                    test_data_igs_id = test_data_temp[test_data_temp['id'] == igs_id]

                    # get consensus sequence for test dataset
                    try:
                        # test_consensus = Bio.motifs.create(test_data_temp['guide'],
                        #                                    alphabet='GATCRYWSMKHBVDN-').degenerate_consensus.strip('-')
                        mat = rd.words2countmatrix(test_data_temp['guide'])
                        test_consensus = pairwise_comparison_consensus = rd.consensus(mat)
                    except ValueError:
                        print(f'{key}: {igs_id} does not have the same length in all guides in test sequences.')
                        continue

                    design_igs_counted = pd.concat([lm.alignment_to_matrix(design_data_igs_id['guide']),
                                                    igs_probs], ignore_index=True).fillna(0)
                    ref_igs_counted = pd.concat([lm.alignment_to_matrix(ref_data_igs_id['guide']),
                                                 igs_probs], ignore_index=True).fillna(0)
                    test_size = test_data_igs_id.shape[0]
                    igs_probs = igs_probs.apply(lambda x: x * test_size)
                    test_igs_counted = pd.concat([lm.alignment_to_matrix(test_data_igs_id['guide']),
                                                  igs_probs], ignore_index=True).fillna(0)

                    # And now get score for best U and ref U with the consensus
                    # there should only be one design per igs, so take the first one
                    design_vs_test_igs = rd.pairwise_comparison(seq_a=test_consensus, opti_len=guide_len,
                                                                seq_b=design_data_igs_id['guide'].values[0],
                                                                score_type='weighted', only_consensus=True)
                    ref_vs_test_igs = rd.pairwise_comparison(seq_a=test_consensus, opti_len=guide_len,
                                                             seq_b=ref_data_igs_id['guide'].values[0],
                                                             score_type='weighted', only_consensus=True)

                    # best_igs_counted = lm.alignment_to_matrix(for_best_finding.nlargest(1, 'composite_test_score')
                    #                                           ['guide'])

                    # Save in a tuple with:
                    # (idx, [extracted_data_test, extracted_data_design, extracted_data_ref], igs_true_%_cov)
                    data_to_graph_igs[key].append((igs_id, [test_igs_counted, design_igs_counted, ref_igs_counted],
                                                   igs_true_perc_cov.iloc[0], design_vs_test_igs, ref_vs_test_igs))
                    # except lm.src.error_handling.LogomakerError:
                    #     print(f'{key}: {igs_id} does not have the same length in all guides in test sequences.')
                    #     # basically if not all guides are the same length. This will happen if minlen was given
                    #     # or if the guide is at the end of the test sequence
                    #     continue
                bar()

    # Now plot!
    titles = {key_names[0]: 'Hign U, high IGS', key_names[1]: 'High U, low IGS', key_names[2]: 'Low U, low IGS',
              key_names[3]: 'Medium U, medium IGS'}

    for key in key_names:
        try:
            random_num = random.sample(range(len(data_to_graph_igs[key])), 1)[0]
        except ValueError:
            print(f'{key} did not have any sequences in one of the datasets given')
            continue
        # for igs_id, data, u_con in data_to_graph_u[key][random_num]:
        igs_id, data, u_con, design_v_test_igs, ref_v_test_igs = data_to_graph_u[key][random_num]
        fig, ax = plt.subplots(4, 1)
        lm.Logo(data[0], ax=ax[0], stack_order='small_on_top', color_scheme='colorblind_safe')
        lm.Logo(data[1], ax=ax[1], stack_order='small_on_top', color_scheme='colorblind_safe')
        lm.Logo(data[2], ax=ax[2], stack_order='small_on_top', color_scheme='colorblind_safe')
        lm.Logo(data[3], ax=ax[3], stack_order='small_on_top', color_scheme='colorblind_safe')

        ax[0].set_title('Testing sequences', loc='left')
        ax[1].set_title(f'Best composite score design: conservation with consensus test: '
                        f'{round(design_v_test_igs, 4)}', loc='left')
        ax[2].set_title('All design sequences', loc='left')
        ax[3].set_title(f'Reference sequence: conservation with consensus test: '
                        f'{round(ref_v_test_igs, 4)}', loc='left')

        fig.suptitle(f'{titles[key]} igs data at position {igs_id}\nU conservation = {round(u_con, 4)}')
        plt.tight_layout()

        if save_fig:
            plt.savefig(fname=f'{save_file_loc}/igs_alignments_at_U_{key}.{file_type}', format=file_type)

        plt.show()

        igs_id, data, igs_con, design_v_test_guide, ref_v_test_guide = data_to_graph_igs[key][random_num]
        fig, ax = plt.subplots(3, 1)
        logo_1 = lm.Logo(data[0], ax=ax[0], stack_order='small_on_top', color_scheme='colorblind_safe')
        logo_2 = lm.Logo(data[1], ax=ax[1], stack_order='small_on_top', color_scheme='colorblind_safe')
        logo_3 = lm.Logo(data[2], ax=ax[2], stack_order='small_on_top', color_scheme='colorblind_safe')

        for i, char in enumerate('G' + igs_id[0:5]):
            logo_1.style_single_glyph(c=char, p=guide_len + i, color=[0.9, 0.9, 0.9])
            logo_2.style_single_glyph(c=char, p=guide_len + i, color=[0.9, 0.9, 0.9])
            logo_3.style_single_glyph(c=char, p=guide_len + i, color=[0.9, 0.9, 0.9])

        ax[0].set_title('Testing sequences', loc='left')
        ax[1].set_title(f'Design sequence at position: conservation with consensus test: '
                        f'{round(design_v_test_guide, 4)}',
                        loc='left')
        ax[2].set_title(f'Reference sequence: conservation with consensus test: '
                        f'{round(ref_v_test_guide, 4)}', loc='left')

        fig.suptitle(f'{titles[key]} guide data for igs {igs_id}\nIGS true % coverage = {round(igs_con, 4)}')
        plt.tight_layout()

        if save_fig:
            plt.savefig(fname=f'{save_file_loc}/guide_alignments_at_igs_{key}.{file_type}', format=file_type)

        plt.show()

    # Graph them in a sequence logo with probability on the y axis
    # And add the consensus motif underneath, and the design at that location, and the ref_sequence at that location
    # Also add the consensus score between the test sequences, between the design and consensus, and between the refseq
    # and the consensus.
    return



def do_hist_plot(fig, ax, graph, xvar, yvar, bins_wanted, min_max, top_control, xlabel, ylabel, tm_min, tm_max):

    # Plot scatter and kde plots
    # # below is for individual log norms
    # hist_data = jointplot_fig.axes[0].hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
    #                                          norm='log')
    hist_data = ax.hist2d(graph[xvar], graph[yvar], cmap='viridis', bins=bins_wanted,
                                             norm=LogNorm(vmin=1, vmax=min_max[f'{xvar}_vs_{yvar}'][-1]))
    fig.colorbar(hist_data[-1], ax=ax)
    ax.scatter(x=top_control[xvar], marker='^', y=top_control[yvar], alpha=1, c='#000000',
                                  label='Control', edgecolors='black', linewidth=0.8,
                                  sizes=[150] * len(top_control[yvar]))
    regstats = scipy.stats.linregress(x=graph[xvar], y=graph[yvar])
    ax.annotate(f'$r^2$={round(regstats.rvalue ** 2, 3)}', xy=(0.1, 0.9),
                                   xycoords='axes fraction')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if yvar == 'tm_nn_vs_test':
        ax.set(ylim=[tm_min, tm_max + 20], xlim=[0, 1])
    else:
        ax.set(ylim=[0, 1], xlim=[0, 1])
    return

def make_violin_plots(folder: str, vars_to_plot: list[str], file_type: str = 'png'):
    print('Now loading data...')
    to_analyze = {}
    files = [(f, f'{folder}/{f}/coupled/results') for f in os.listdir(folder) if f != '.DS_Store']
    sorted_files = sorted(files, key=lambda a: int(a[0][:-3]))
    for name, file in sorted_files:
        try:
            items = [file + '/' + data_file for data_file in os.listdir(file) if data_file.endswith('.txt')]
            to_analyze[name.replace('_', ' ')] = items
        except FileNotFoundError:
            continue

    with alive_bar(total=len(to_analyze), spinner='fishes') as bar:
        all_data_df = None
        for key, file in to_analyze.items():
            next_design_df = import_data_to_df(file, key, True)
            if all_data_df is None:
                all_data_df = next_design_df
            else:
                all_data_df = pd.concat([all_data_df, next_design_df])
            bar()
    print('Data loaded!')
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 20 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')

    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=3 * len(vars_to_plot), ncols=7)
    joint_ax = {}
    for i, var in enumerate(vars_to_plot):
        joint_ax[i] = jointplot_fig.add_subplot(gridspec[i * 3:i * 3 + 3, 0:7])
        # sns.boxplot(x=all_data_df.index, y=var, data=all_data_df, notch=True, boxprops={'alpha': 0.7}, whis=(0, 100),
        #             fill=False, ax=joint_ax[i], color='#000000')
        sns.violinplot(x=all_data_df.index, y=var, data=all_data_df, cut=0, split=True, ax=joint_ax[i], inner=None,
                       color='#000000')
        # sns.stripplot(x=all_data_df.index, y=var, linewidth=0, data=all_data_df, dodge=True, alpha=0.01, ax=joint_ax[i],
        #               color='gray')
        joint_ax[i].set(xlabel='Guide length', ylabel=var.replace('_', ' '))
        joint_ax[i].label_outer()
    jointplot_fig.suptitle('Scores vs. guide length')
    text = '_and_'.join(vars_to_plot)
    plt.tight_layout()
    plt.savefig(fname=f'{folder}/{text}.{file_type}', format=file_type)
    plt.show()

    return

# Delta difference plot between phyla?? for selective??
