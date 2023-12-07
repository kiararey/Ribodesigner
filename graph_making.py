import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from alive_progress import alive_bar
import logomaker


def make_graphs(var_regs: list[tuple[int, int]], control_designs_path: str, universal_designs_path: str = '',
                selective_designs_path: str = '',
                ref_seq_designs_path: str = '', save_fig: bool = False,
                file_type: str = 'png', save_file_loc: str = ''):
    print('Loading data...')
    fn = lambda x: int(x.split('/')[-1].split('_')[0])
    percentages = []
    for path in [control_designs_path, universal_designs_path, selective_designs_path, ref_seq_designs_path]:
        if path:
            percentages.append(fn(path))
        else:
            percentages.append(0)
    total_designs = sum(percentages)
    # percentages = [perc/total_designs for perc in percentages]
    with alive_bar(total=total_designs, spinner='fishes') as bar:
        control_designs_df = extract_info(control_designs_path, 'control')
        bar(percentages[0])
        print('Control data loaded!')
        if universal_designs_path:
            print('Now loading universal data...')
            universal_designs_df = extract_info(universal_designs_path, 'universal')
            bar(percentages[1])
            print('Universal data loaded!')
            all_data_df = pd.concat([control_designs_df, universal_designs_df])
        else:
            all_data_df = control_designs_df
        if selective_designs_path:
            print('Now loading selective data...')
            selective_designs_df = extract_info(selective_designs_path, 'selective')
            bar(percentages[2])
            print('Selective data loaded!')
            all_data_df = pd.concat([all_data_df, selective_designs_df])
        if ref_seq_designs_path:
            print('Now loading reference data...')
            ref_seq_designs_df = extract_info(ref_seq_designs_path, 'ref_seq')
            bar(percentages[3])
            print('Reference sequence data loaded!')
            all_data_df = pd.concat([all_data_df, ref_seq_designs_df])

    control_designs_df_1376 = control_designs_df[control_designs_df['reference_idx'] == 1376]

    # get top scores in each: 1 control, one of each group of selective and universal
    top_score_control = control_designs_df_1376[control_designs_df_1376[
        'name_of_test_dataset'].str.contains('Bac')].nlargest(1, 'composite_test_score')
    file_name_checks = ['Bac', 'Euk', 'Arc', 'All']
    for name in file_name_checks[1:]:
        top_score_temp = control_designs_df_1376[control_designs_df_1376[
            'name_of_test_dataset'].str.contains(name)].nlargest(1, 'composite_test_score')
        top_score_control = pd.concat([top_score_control, top_score_temp])

    top_scores_universal = universal_designs_df.loc[universal_designs_df.index ==
                                                    'universal'].nlargest(1, 'composite_test_score')
    print('Generating graphs...')
    # for label in universal_labels[1:]:
    #     score_temp = universal_designs_df[universal_designs_df['tested_targets'] ==
    #                                       label].nlargest(1, 'composite_test_score')
    #     top_scores_universal = pd.concat([top_scores_universal, score_temp])
    #
    # top_scores_selective = selective_designs_df.loc[selective_designs_df['tested_targets'] ==
    #                                                 selective_labels[0]].nlargest(1, 'composite_test_score')
    # for label in selective_labels[1:]:
    #     score_temp = selective_designs_df.loc[selective_designs_df['tested_targets'] ==
    #                                           label].nlargest(1, 'composite_test_score')
    #     top_scores_selective = pd.concat([top_scores_selective, score_temp])

    # top_test_scores = pd.concat([top_score_control, top_scores_universal, top_scores_selective])
    top_test_scores = pd.concat([top_score_control, top_scores_universal])

    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 16 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # colors = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']  # viridis
    colors = ['#0D365E', '#3F6D54', '#9B892D', '#F9A281', '#FACDFB']  # batlow

    # target_names = ['reference', 'bacterial', 'archaeal', 'eukaryotic', 'all kingdoms']
    target_names = ['reference', 'bacterial', 'archaeal', 'eukaryotic', 'all kingdoms']
    to_analyze = ['ref_seq', 'universal', 'control']

    max_vals = all_data_df.max(numeric_only=True)

    # Fig 1: MG1655
    test_name = 'vs_test_sequences_Bac'
    reference_subset = ref_seq_designs_df[ref_seq_designs_df['name_of_test_dataset'].str.contains(test_name)]
    bacterial_subset = universal_designs_df[universal_designs_df['name_of_test_dataset'].str.contains(test_name)]
    bacterial_subset = bacterial_subset.loc[bacterial_subset.index == 'universal']
    top_score_control_subset = top_score_control[top_score_control['name_of_test_dataset'].str.contains(test_name)]

    # testing_subset = 'All'
    # reference_subset = ref_seq_designs_df[ref_seq_designs_df['name_of_test_dataset'].str.contains(testing_subset)]
    # bacterial_subset = universal_designs_df[universal_designs_df['name_of_test_dataset'].str.contains(testing_subset)]
    # bacterial_subset = bacterial_subset.loc[bacterial_subset.index == 'universal']
    # top_score_control_subset = top_score_control[top_score_control['name_of_test_dataset'].str.contains(testing_subset)]
    # Only going to plot bacteria test dataset
    # First, plot x = u conservation, y = igs conservation, hue = guide score

    # Prepare axes
    for graph, name in zip([reference_subset, bacterial_subset], ['reference seq', 'bacterial']):
        # Prepare axes
        jointplot_fig = plt.figure()
        gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
        joint_ax = {
            0: jointplot_fig.add_subplot(gridspec[1:7, 7:13]),
            1: jointplot_fig.add_subplot(gridspec[0:1, 7:13]),
            2: jointplot_fig.add_subplot(gridspec[1:7, 13:14]),
            3: jointplot_fig.add_subplot(gridspec[0:2, 0:6]),
            4: jointplot_fig.add_subplot(gridspec[2:4, 0:6]),
            5: jointplot_fig.add_subplot(gridspec[4:6, 0:6])
        }
        xvar = 'u_conservation_test'
        yvar = 'true_%_cov_test'
        hue_var = 'test_score'
        # Plot scatter and kde plots
        # sns.scatterplot(x=xvar, y=yvar, hue=hue_var, data=graph, ax=joint_ax[0], alpha=0.5, legend=False)
        jointplot_fig.axes[0].scatter(x=graph[xvar], y=graph[yvar], alpha=0.5, c=graph[hue_var], vmin=0, vmax=1,
                                      cmap='viridis')
        jointplot_fig.axes[0].scatter(x=top_score_control_subset[xvar], vmin=0, vmax=1, marker='^',
                                      y=top_score_control_subset[yvar], alpha=1, c=top_score_control_subset[hue_var],
                                      edgecolors='#9B892D', label='Control')
        sns.kdeplot(x=xvar, data=graph, ax=joint_ax[1], fill=True, common_norm=True, alpha=.3, legend=False)
        sns.kdeplot(y=yvar, data=graph, ax=joint_ax[2], fill=True, common_norm=True, alpha=.3, legend=False)
        jointplot_fig.axes[0].legend(bbox_to_anchor=(1.6, -0.15), title='Guide score')
        jointplot_fig.axes[0].set_xlabel('U conservation')
        jointplot_fig.axes[0].set_ylabel('IGS true percent coverage')
        # Set variable regions for location plots
        plot_variable_regions(joint_ax[3], var_regs)
        plot_variable_regions(joint_ax[4], var_regs)
        plot_variable_regions(joint_ax[5], var_regs)
        # Plot test data for each testing condition
        sns.scatterplot(x='reference_idx', y='u_conservation_test', data=graph, ax=joint_ax[3], alpha=0.5, legend=False)
        sns.scatterplot(x='reference_idx', y='true_%_cov_test', data=graph, ax=joint_ax[4], alpha=0.5, legend=False)
        sns.scatterplot(x='reference_idx', y='test_score', data=graph, ax=joint_ax[5], alpha=0.5)
        jointplot_fig.axes[4].set_xlabel('16s rRNA sequence position on reference sequence')
        # Plot control data
        jointplot_fig.axes[3].scatter(x=top_score_control_subset['reference_idx'],
                                      y=top_score_control_subset['u_conservation_test'],
                                      alpha=1, c='#fde725', label='control')
        jointplot_fig.axes[3].set_ylabel('U site conservation')
        jointplot_fig.axes[4].scatter(x=top_score_control_subset['reference_idx'],
                                      y=top_score_control_subset['u_conservation_test'],
                                      alpha=1, c='#fde725', label='control')
        jointplot_fig.axes[4].set_ylabel('IGS true percent coverage')
        jointplot_fig.axes[5].scatter(x=top_score_control_subset['reference_idx'],
                                      y=top_score_control_subset[yvar], alpha=1,
                                      c='#fde725', label='control')
        jointplot_fig.axes[5].set_ylabel('Guide score')
        # Set graph settings for pretti graphing
        jointplot_fig.axes[0].set(xlim=[-0.1, 1.1], ylim=[-0.1, 1.1])
        jointplot_fig.axes[0].set(ylim=[-0.1, 1.1])
        jointplot_fig.axes[1].set(xlabel=None)
        # jointplot_fig.axes[1].set_title('D', loc='left', fontsize=30)
        jointplot_fig.axes[2].set(ylabel=None)
        # jointplot_fig.axes[3].set_title('A', loc='left', fontsize=30)
        jointplot_fig.axes[3].set(xlabel=None, xlim=[-0.1, max_vals['reference_idx'] + 20], ylim=[-0.1, 1.1])
        jointplot_fig.axes[4].sharex(jointplot_fig.axes[3])
        jointplot_fig.axes[4].sharey(jointplot_fig.axes[3])
        jointplot_fig.axes[5].sharex(jointplot_fig.axes[3])
        jointplot_fig.axes[5].sharey(jointplot_fig.axes[3])
        jointplot_fig.axes[4].set(xlabel=None)
        # jointplot_fig.axes[4].set_title('B', loc='left', fontsize=30)
        # jointplot_fig.axes[5].set_title('C', loc='left', fontsize=30)
        jointplot_fig.axes[5].set(xlabel='Reference 16s rRNA index')
        jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
        jointplot_fig.axes[1].tick_params(labelbottom=False, labelleft=False, left=False)
        jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
        jointplot_fig.axes[2].tick_params(labelbottom=False, labelleft=False, bottom=False)
        jointplot_fig.axes[3].tick_params(labelbottom=False)
        jointplot_fig.axes[4].tick_params(labelbottom=False)
        jointplot_fig.suptitle(f'Designs from {name} target sequences against bacterial test dataset')

        if save_fig:
            plt.savefig(fname=f'{save_file_loc}/figure_1_{name}.{file_type}', format=file_type)
        plt.show()

    # Fig 2: : Assessing universal design quality. 2a) IGS true percent coverage vs. guide score of bacterial designs
    # evaluated against datasets of different kingdoms. 2b) 16s rRNA location of all designs along the reference
    # E. coli MG1655 16s rRNA sequence.
    for i, target_name in enumerate(target_names):
        universal_subset = all_data_df.loc[all_data_df.index == to_analyze[i]]

        # Prepare axes
        jointplot_fig = plt.figure()
        gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
        joint_ax = {
            0: jointplot_fig.add_subplot(gridspec[1:7, 7:13]),
            1: jointplot_fig.add_subplot(gridspec[0:1, 7:13]),
            2: jointplot_fig.add_subplot(gridspec[1:7, 13:14]),
            3: jointplot_fig.add_subplot(gridspec[0:2, 0:6]),
            4: jointplot_fig.add_subplot(gridspec[2:4, 0:6]),
            5: jointplot_fig.add_subplot(gridspec[4:6, 0:6])
        }
        xvar = 'true_%_cov_test'
        yvar = 'test_score'
        legend_names = ['Eukaryota only', 'Archaea only', 'Bacteria only', 'All kingdoms', 'Control']
        # Plot scatter and kde plots
        sns.scatterplot(x=xvar, y=yvar, hue='name_of_test_dataset', data=universal_subset, ax=joint_ax[0], alpha=0.5,
                        legend=False)
        jointplot_fig.axes[0].scatter(x=top_score_control[xvar], y=top_score_control[yvar],
                                      alpha=1, c='#fde725', label='Control')
        sns.kdeplot(x=xvar, hue='name_of_test_dataset', data=universal_subset, ax=joint_ax[1], fill=True,
                    common_norm=True, alpha=.3, legend=False)
        sns.kdeplot(y=yvar, hue='name_of_test_dataset', data=universal_subset, ax=joint_ax[2], fill=True,
                    common_norm=True, alpha=.3, legend=False)
        jointplot_fig.axes[0].set_xlabel('IGS true percent coverage')
        jointplot_fig.axes[0].set_ylabel('Guide score')
        # Set variable regions for location plots
        plot_variable_regions(joint_ax[3], var_regs)
        plot_variable_regions(joint_ax[4], var_regs)
        plot_variable_regions(joint_ax[5], var_regs)
        # Plot test data for each testing condition
        sns.scatterplot(x='reference_idx', y='u_conservation_test', hue='name_of_test_dataset', data=universal_subset,
                        ax=joint_ax[3],
                        alpha=0.5, legend=False)
        sns.scatterplot(x='reference_idx', y=xvar, hue='name_of_test_dataset', data=universal_subset, ax=joint_ax[4],
                        alpha=0.5, legend=False)
        sns.scatterplot(x='reference_idx', y=yvar, hue='name_of_test_dataset', data=universal_subset, ax=joint_ax[5],
                        alpha=0.5)
        jointplot_fig.axes[4].set_xlabel('16s rRNA sequence position on reference sequence')
        # Plot control data
        jointplot_fig.axes[3].scatter(x=top_score_control['reference_idx'], y=top_score_control['u_conservation_test'],
                                      alpha=1, c='#fde725', label='control')
        jointplot_fig.axes[3].set_ylabel('U site conservation')
        jointplot_fig.axes[4].scatter(x=top_score_control['reference_idx'], y=top_score_control[xvar],
                                      alpha=1, c='#fde725', label='control')
        jointplot_fig.axes[4].set_ylabel('IGS true percent coverage')
        jointplot_fig.axes[5].scatter(x=top_score_control['reference_idx'], y=top_score_control[yvar], alpha=1,
                                      c='#fde725', label='control')
        jointplot_fig.axes[5].set_ylabel('Guide score')
        # Set graph settings for pretti graphing
        jointplot_fig.axes[0].set(xlim=[-0.1, 1.1], ylim=[-0.1, 1.1])
        jointplot_fig.axes[0].set(ylim=[-0.1, 1.1])
        jointplot_fig.axes[1].set(xlabel=None)
        # jointplot_fig.axes[1].set_title('D', loc='left', fontsize=30)
        jointplot_fig.axes[2].set(ylabel=None)
        # jointplot_fig.axes[3].set_title('A', loc='left', fontsize=30)
        jointplot_fig.axes[3].set(xlabel=None, xlim=[-0.1, max_vals['reference_idx'] + 20], ylim=[-0.1, 1.1])
        jointplot_fig.axes[4].sharex(jointplot_fig.axes[3])
        jointplot_fig.axes[4].sharey(jointplot_fig.axes[3])
        jointplot_fig.axes[5].sharex(jointplot_fig.axes[3])
        jointplot_fig.axes[5].sharey(jointplot_fig.axes[3])
        jointplot_fig.axes[4].set(xlabel=None)
        # jointplot_fig.axes[4].set_title('B', loc='left', fontsize=30)
        # jointplot_fig.axes[5].set_title('C', loc='left', fontsize=30)
        jointplot_fig.axes[5].set(xlabel='Reference 16s rRNA index')
        jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
        jointplot_fig.axes[1].tick_params(labelbottom=False, labelleft=False, left=False)
        jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
        jointplot_fig.axes[2].tick_params(labelbottom=False, labelleft=False, bottom=False)
        jointplot_fig.axes[3].tick_params(labelbottom=False)
        jointplot_fig.axes[4].tick_params(labelbottom=False)
        legend = plt.legend(bbox_to_anchor=(1.4, -0.15), ncols=5, title='Test dataset')
        for label in legend.get_texts():
            txt = label.get_text().split('/')[-1].split('vs_test_sequences')[-1].replace('_', ' ')
            new_label = [name for name in legend_names if name.casefold() in txt.casefold()]
            label.set_text(new_label[0])
        jointplot_fig.suptitle(f'Designs from {target_name} target sequences against test datasets')

        if save_fig:
            plt.savefig(fname=f'{save_file_loc}/figure_2_{target_name}.{file_type}', format=file_type)
        plt.show()

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


def make_split_graph(title: str, x_data: list, xlabel: str, y_data: list, ylabel: str, hue_data: list, huelabel: str,
                     loc_data: list[int], var_regs: list, save_file_name: str, file_type: str = 'png',
                     add_control: bool = False, control_x_data: list = None, control_y_data: list = None,
                     control_hue_data: list = None, control_loc_data: list[int] = None, alpha=0.5, control_alpha=1,
                     vmin: float = 0.0, vmax: float = 1.0):
    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 16 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')

    # Prepare axes
    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
    joint_ax = {
        0: jointplot_fig.add_subplot(gridspec[1:7, 7:13]),
        1: jointplot_fig.add_subplot(gridspec[0:1, 7:13]),
        2: jointplot_fig.add_subplot(gridspec[1:7, 13:14]),
        3: jointplot_fig.add_subplot(gridspec[0:3, 0:6]),
        4: jointplot_fig.add_subplot(gridspec[3:6, 0:6])
    }
    # Plot scatter and kde plots
    # sns.scatterplot(x=xvar, y=yvar, hue=hue_var, data=graph, ax=joint_ax[0], alpha=0.5, legend=False)
    if hue_data:
        joint_ax[0].scatter(x=x_data, y=y_data, alpha=alpha, c=hue_data, vmin=vmin, vmax=vmax, cmap='viridis')
    else:
        joint_ax[0].scatter(x=x_data, y=y_data, alpha=alpha, vmin=vmin, vmax=vmax, cmap='viridis')
    if add_control:
        if control_hue_data:
            jointplot_fig.axes[0].scatter(x=control_x_data, c=control_hue_data, marker='^', y=control_y_data,
                                          alpha=control_alpha, edgecolors='#9B892D', label='Control', vmin=0, vmax=1)
        else:
            jointplot_fig.axes[0].scatter(x=control_x_data, marker='^', y=control_y_data,
                                          alpha=control_alpha, edgecolors='#9B892D', label='Control', vmin=0, vmax=1)
    sns.kdeplot(x=x_data, ax=joint_ax[1], fill=True, common_norm=True, alpha=.3, legend=False)
    sns.kdeplot(y=y_data, ax=joint_ax[2], fill=True, common_norm=True, alpha=.3, legend=False)
    jointplot_fig.axes[0].legend(bbox_to_anchor=(1.6, -0.15), title=huelabel)
    jointplot_fig.axes[0].set_xlabel(xlabel)
    jointplot_fig.axes[0].set_ylabel(ylabel)
    # Set variable regions for location plots
    plot_variable_regions(joint_ax[3], var_regs)
    plot_variable_regions(joint_ax[4], var_regs)

    # Plot test data for each testing condition
    # jointplot_fig.axes[3].scatter(x=loc_data, y=x_data, alpha=alpha, c=y_data, vmin=0, vmax=1, cmap='viridis')
    # jointplot_fig.axes[4].scatter(x=loc_data, y=y_data, alpha=alpha, c=x_data, vmin=0, vmax=1, cmap='viridis')
    jointplot_fig.axes[3].scatter(x=loc_data, y=x_data, alpha=alpha / 2, vmin=0, vmax=1, cmap='viridis')
    jointplot_fig.axes[4].scatter(x=loc_data, y=y_data, alpha=alpha, vmin=0, vmax=1, cmap='viridis')

    # Plot control data
    if add_control:
        jointplot_fig.axes[3].scatter(x=control_loc_data, y=control_x_data, alpha=control_alpha, c='#fde725',
                                      label='control')
        jointplot_fig.axes[4].scatter(x=control_loc_data, y=control_y_data, alpha=control_alpha, c='#fde725',
                                      label='control')

    # Set graph settings for pretti graphing
    jointplot_fig.axes[3].set_ylabel(xlabel)
    jointplot_fig.axes[4].set_ylabel(ylabel)
    jointplot_fig.axes[0].set(xlim=[-0.1, 1.1], ylim=[-0.1, 1.1])
    jointplot_fig.axes[0].set(ylim=[-0.1, 1.1])
    jointplot_fig.axes[1].set(xlabel=None)
    # jointplot_fig.axes[1].set_title('D', loc='left', fontsize=30)
    jointplot_fig.axes[2].set(ylabel=None)
    # jointplot_fig.axes[3].set_title('A', loc='left', fontsize=30)
    if add_control:
        temp_loc = control_loc_data.extend(loc_data)
        max_loc = max(temp_loc)
    else:
        max_loc = max(loc_data)
    jointplot_fig.axes[3].set(xlabel=None, xlim=[-0.1, max_loc + 20], ylim=[-0.1, 1.1])
    jointplot_fig.axes[4].sharex(jointplot_fig.axes[3])
    jointplot_fig.axes[4].sharey(jointplot_fig.axes[3])
    jointplot_fig.axes[4].set(xlabel=None)
    # jointplot_fig.axes[4].set_title('B', loc='left', fontsize=30)
    # jointplot_fig.axes[5].set_title('C', loc='left', fontsize=30)
    jointplot_fig.axes[4].set(xlabel='Reference 16s rRNA index')
    jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[1].tick_params(labelbottom=False, labelleft=False, left=False)
    jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
    jointplot_fig.axes[2].tick_params(labelbottom=False, labelleft=False, bottom=False)
    jointplot_fig.axes[3].tick_params(labelbottom=False)
    # jointplot_fig.axes[4].tick_params(labelbottom=False)
    jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} '
                                                     f'distributions n = {len(loc_data)}')

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
                    'test_score': float, '%_coverage_test': float, '%_on target_test': float, 'true_%_cov_test': float,
                    'composite_test_score': float, 'delta_igs_vs_test': float, 'delta_guide_vs_test': float,
                    'delta_vs_test': float, 'target_Domain': str, 'target_Domain_background': str,
                    'target_Domain_test': str, 'target_Phylum': str, 'target_Phylum_background': str,
                    'target_Phylum_test': str, 'target_Class': str, 'target_Class_background': str,
                    'target_Class_test': str, 'target_Order': str, 'target_Order_background': str,
                    'target_Order_test': str, 'target_Family': str, 'target_Family_background': str,
                    'target_Family_test': str, 'target_Genus': str, 'target_Genus_background': str,
                    'target_Genus_test': str, 'target_Species': str, 'target_Species_background': str,
                    'target_Species_test': str, 'target_Taxon': str, 'target_Taxon_background': str,
                    'target_Taxon_test': str}
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


def make_sequence_logo_graph(data):
    # Pick out 3 positions at edge cases and one in the middle

    # Find designs and ref_sequence guides at those positions

    # Ok, now find all IGSes in the dataset at that position

    # Graph them in a sequence logo with probability on the y axis
    # And add the consensus motif underneath, and the design at that location, and the ref_sequence at that location
    # Also add the consensus score between the test sequences, between the design and consensus, and between the refseq
    # and the consensus.

    # Find all guide sequences in the dataset at that position

    # Pull out corresponding design and ref sequence

    # Graph them in a sequence logo with probability on the y axis
    # And add the consensus motif underneath, and the design at that location, and the ref_sequence at that location
    # Also add the consensus score between the test sequences, between the design and consensus, and between the refseq
    # and the consensus.
    return