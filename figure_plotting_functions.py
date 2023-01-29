from RiboDesigner import RiboDesigner, read_fasta_folder, read_fasta, find_cat_sites, find
from Bio.Seq import Seq
from Bio import pairwise2
import re
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import pandas as pd
import os
import math
import cmcrameri.cm as cmc
import numpy as np

pd.options.mode.chained_assignment = None


def insert_into_title(title, key_word, to_insert):
    index = title.find(key_word)
    newtitle = title[:index] + to_insert + title[index:]
    return newtitle


def fix_names(names):
    index = names.find(' ')
    new_names = names[0] + '.' + names[index:]
    return new_names


def plot_variable_regions(ax, var_regs):
    for V in var_regs:
        ax.axvspan(V[0], V[1], facecolor='g', alpha=0.2)
    return


def make_plot(target_data, title, filename, var_regs, colr, splt=False, on_target=[], off_target=[]):
    fig, ax = plt.subplots(2, 1, sharey='col', constrained_layout=True)

    names = [tup[0].replace('_', ' ') for tup in target_data]
    refidx = [tup[1] for tup in target_data]
    IGS = [tup[2] for tup in target_data]
    score = [tup[3] for tup in target_data]
    repeat_cols = [tup[4] for tup in target_data]
    off_target = [tup[5] for tup in target_data]

    plot_variable_regions(ax[0], var_regs)

    ax[0].set(xlim=(-0.1, 1599), ylim=(0, 1.1))

    # plot score vs. site of U (ref--> E.coli)
    ax[0].scatter(refidx, score, alpha=0.5, color=repeat_cols, edgecolors=off_target)
    ax[0].set_ylabel('Score')
    ax[0].set_xlabel('Location of U on E. coli 16s rRNA (base pair location)')
    ax[0].set_title(title)
    sns.despine(ax=ax[0], offset=5)

    # prepare names for plotting
    new_names = []
    for name in names:
        new_names.append(fix_names(name))

    # plot score vs. target name
    # sns.set_palette('viridis')

    plot_variable_regions(ax[1], var_regs)

    ax[1].set_ylabel('Score')
    ax[1].set_xlabel('Target name')
    ax[1].tick_params(bottom=False)
    if splt:
        score_on = [tup[3] for tup in on_target]
        names_on = [fix_names(tup[0].replace('_', ' ')) for tup in on_target]
        score_off = [tup[3] for tup in off_target]
        names_off = [fix_names(tup[0].replace('_', ' ')) for tup in off_target]

        colms = ['on off', 'Score', 'Target name']
        off_rows = pd.DataFrame(
            data={colms[0]: ['Off target'] * len(score_off), colms[1]: score_off, colms[2]: names_off})
        on_rows = pd.DataFrame(data={colms[0]: ['On target'] * len(score_on), colms[1]: score_on, colms[2]: names_on})
        all_rows = pd.concat([on_rows, off_rows], axis=0, join='outer')
        ax[1] = sns.violinplot(x='Target name', y='Score', hue='on off', palette=colr, data=all_rows, split=splt,
                               inner=None, cut=0)
        ax[1].legend([], [], frameon=False)
    else:
        ax[1] = sns.violinplot(x=new_names, y=score, color=colr, inner=None, cut=0)
        # ax[1] = sns.swarmplot(x=new_names, y=score, color='white', edgecolor='gray', s=0.8)
    sns.despine(ax=ax[1], offset=5, bottom=True)

    plt.savefig(filename, transparent=True)
    return


def align_to_ref(n, file, barcode_seq_file, ribobody_file, target_path, ref_path, riboscore_test=False):
    # n is the guide length
    # file is where the data output will be saved
    # barcode_seq_file is the FASTA or .txt file with the barcode sequence
    # ribobody_file is the FASTA or .txt file with the Ribozyme body sequence
    # target_path is the folder where the data to check is
    # ref_path is the folder where the reference genome is
    # riboscore_test is a Boolean that tells if we also want to dcompare NGS vs. RiboDesigner scores

    # variable regions V1-V9 (start, end) 1-based indexing:
    var_regs = [(68, 100), (137, 226), (440, 496), (590, 650), (829, 856), (999, 1037), (1119, 1156),
                (1243, 1295), (1435, 1465)]

    # set up plots to look pretty
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params, palette='viridis')
    sns.set_context("paper")
    # use this website to get this color map https://waldyrious.net/viridis-palette-generator/
    color_dict = {0: '#440154', 1: '#3b528b', 2: '#21918c', 3: '#5ec962', 4: '#fde725'}

    # this is U64
    ribozyme_to_test = Seq("CANCCCACTCCCATGGTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGTGTTCAC").upper().transcribe()
    ribozyme_to_test_index = 1376

    # # these are U8, U34, and U73
    # otherToTest = [Seq("ATATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGNCCGTGTCTCAGTNCCGGTGTG").upper().transcribe(),
    #                Seq("GCGTTGTGCGGGCCCCCGTCAATTCATTTGAGTTTTANNCTTGCGGCCGTGCTCCC").upper().transcribe(),
    #                Seq("AGGAGGTGATCCAGCCNCAGGTTCCCCTANGGCTACCTTGTTACGACTTCGCCCCA").upper().transcribe()]

    # prepare barcode sequence
    with open(barcode_seq_file) as f:
        for i in f:
            bar_seq = Seq(i).upper().transcribe()

    # prepare Ribozyme template sequence
    with open(ribobody_file) as f:
        for i in f:
            ribo_body = Seq(i).upper().transcribe()

    # set variables
    ribo_seq = ribo_body + bar_seq
    m = 5  # IGS length
    minlen = n

    # import original sequences
    target_names, target_seqs = read_fasta_folder(target_path)

    # run ribodesigner but get the final raw list
    ############## UPDATE THIS!!! ##############
    ranked_and_sorted_igs, sorted_opti_seqs = RiboDesigner(target_seqs, target_names, [], [], ribo_body, bar_seq, m, n,
                                                           minlen, file, fileout=False, thresh=0.8, idxOut=True)

    new_data = ranked_and_sorted_igs[['IGS sequence', '%' + ' coverage', 'Target', 'Occurrences in Target Sequence',
                                      'Index of Splice Site']]
    print('\n')

    # find idx of each U
    good_ribozyme_seqs, good_ig_sand_guide_seqs, good_ig_sseqs, \
    good_guide_seqs, goodidx_list, new_names = find_cat_sites(target_seqs, target_names, ribo_seq, m, n, minlen)

    # Get e coli sequence
    ref_name, ref_seq = read_fasta(ref_path)

    # prepare patterns to look for to extract individual sequences from pairwise alignment
    pattern_a = 'seqA=\'(.*?)\''
    pattern_b = 'seqB=\'(.*?)\''

    # initialize data structures
    new_aligns = []
    new_sequence_list = []
    igs_for_later = []
    color_to_use = color_dict[1]
    conversion_dicts = {name: {} for name in new_names}

    for i in range(len(target_seqs)):
        current_dict = conversion_dicts[new_names[i]]
        print('Now re-indexing ' + new_names[i].replace('_', ' ') + ' to reference ' + ref_name[0].replace('_', ' '))
        # will have to keep in mind the potential lengths of the sequences and add length m to our final E.coli index
        alignments = pairwise2.align.globalxx(target_seqs[i][m:-minlen], ref_seq[0][m:-minlen])
        new_aligns.append(''.join(str(alignments[0])))

        seq_a = re.search(pattern_a, new_aligns[i]).group(1)
        seq_b = re.search(pattern_b, new_aligns[i]).group(1)

        # obtain index of new Us
        idx_seq_a = find(seq_a, 'U')

        # get sub dataframe with only the analyzed sequence
        sub_data = new_data[new_data['Target'] == new_names[i]]

        temp_idx_storage = []
        cnt = 0
        for idx in idx_seq_a:
            # find what index that is based on the reference sequence
            ref_string = seq_b[:idx]
            refidx = len(ref_string.replace('-', '')) + m + 1  # turns zero based indexing to refseq numbering

            # Use location of index to pair with old U idx
            og_idx = goodidx_list[i][cnt]
            coverage = sub_data.loc[sub_data['Index of Splice Site'] == og_idx, ['%' + ' coverage',
                                                                                 'Occurrences in Target Sequence']]
            current_dict[og_idx] = refidx  # fill a dictionary containing the conversion
            # obtain % coverage for each old U and assign it to new U
            if len(coverage.index) == 0:  # since RiboDesigner doesn't bother with unique sequences unless asked to
                temp_idx_storage.append((refidx, 1 / len(new_names), 'lightgrey', 'lightgrey'))
            elif len(coverage.index) > 1:
                print('WARNING: Two indexes found per redidx')
            else:
                cov = coverage.iloc[0, 0]
                # Check if the IGS is degenerate (occurs in multiple spots in the same target sequence) or unique
                check = coverage.iloc[0, 1] > 1
                if check:  # only fill circle if degenerate IGS
                    if cov == 1:  # if the percent coverage is 100% use color to line and fill
                        color_cov = color_to_use
                        occur = color_to_use
                    else:  # else line and fill with grey
                        color_cov = 'grey'
                        occur = 'grey'

                else:  # if unique IGS, do not fill circle
                    if cov == 1:  # if the percent coverage is 100% use color to line
                        color_cov = color_to_use
                        occur = 'white'
                    else:  # else line with grey
                        color_cov = 'grey'
                        occur = 'white'

                # store info as: (refseq idx, % coverage, inner color, line color)
                temp_idx_storage.append((refidx, cov, occur, color_cov))
                if cov == 1:
                    igs_series = sub_data.loc[sub_data['Index of Splice Site'] == og_idx, 'IGS sequence']
                    for igs in igs_series:
                        igs_for_later.append((new_names[i], refidx, og_idx, igs))

            cnt += 1
        new_sequence_list.append(temp_idx_storage)

    # plot %coverage vs. site of U (ref--> E.coli)
    fig, ax = plt.subplots()

    # plot variable regions
    plot_variable_regions(ax, var_regs)

    for i in range(len(new_names)):
        x = [tup[0] for tup in new_sequence_list[i]]
        y = [tup[1] for tup in new_sequence_list[i]]
        col = [tup[2] for tup in new_sequence_list[i]]
        edg = [tup[3] for tup in new_sequence_list[i]]
        ax.scatter(x, y, color=col, alpha=0.5, edgecolors=edg)

    ax.set_ylabel('%' + ' coverage')
    ax.set_xlabel('Location of U on E. coli 16s rRNA (base pair location)')
    ax.set_title('%' + ' coverage of generated designs on selected organisms')

    sns.despine(offset=10, trim=False)
    fig.savefig(file + '/coverage.pdf', transparent=True)

    off_target_data = []
    on_target_data = []
    all_data = []

    # these will be used for comparing scores
    reference_scores = []
    reference_indexes = []

    for i in range(len(igs_for_later)):
        name = igs_for_later[i][0]
        refidx = igs_for_later[i][1]
        og_idx = igs_for_later[i][2]
        igs = igs_for_later[i][3]  # extract IGS

        # extract scores matching this IGS
        scores_df = (sorted_opti_seqs.loc[sorted_opti_seqs['IGS'] == igs,
                                          ['Score', 'Repeats', 'Indexes', 'Guide + G + IGS']])
        scores = scores_df.values.tolist()

        for row in scores:
            # check if there are repeats
            if row[1] == 'No repeats' and row[3] != ribozyme_to_test:
                repeat = 'nan'
            elif row[3] == ribozyme_to_test and refidx == ribozyme_to_test_index:  # if it's U64
                repeat = 'r'
            # elif row[3] in otherToTest:  # if it's any of the other sequences
            #     repeat = 'r'
            else:
                repeat = 'w'
            # check whether indexes match
            index_info = row[2]

            degen_test = []
            for key in index_info.keys():
                idx_to_test = index_info[key]
                idx_ref = conversion_dicts[key][idx_to_test]
                degen_test.append(idx_ref)

                if key == 'Escherichia_coli':
                    reference_scores.append(row[0])
                    reference_indexes.append(idx_ref)

            # count how many off target indexes there are and store that data
            off_target = len(set(degen_test)) - 1
            if off_target == 0:
                if repeat == 'nan':
                    repeat = color_dict[off_target]
                    off_target_color = color_dict[off_target]
                else:
                    off_target_color = color_dict[off_target]

                on_target_data.append((name, refidx, igs, row[0], repeat, off_target_color))
            else:
                if repeat == 'nan':
                    repeat = color_dict[off_target]
                    off_target_color = color_dict[off_target]
                else:
                    off_target_color = color_dict[off_target]

                off_target_data.append((name, refidx, igs, row[0], repeat, off_target_color))

            all_data.append((name, refidx, igs, row[0], repeat, off_target_color))

    title_template = 'Score of optimized Ribozyme designs with 100' + '%' + ' coverage'
    key_word = 'optimized'

    make_plot(on_target_data, insert_into_title(title_template, key_word, 'on-target '), file + '/on target data.pdf',
              var_regs,
              color_dict[0])
    make_plot(off_target_data, insert_into_title(title_template, key_word, 'off-target '),
              file + '/off target data.pdf',
              var_regs, 'lightgrey')
    make_plot(all_data, insert_into_title(title_template, key_word, 'all '), file + '/all target data.pdf',
              var_regs, [color_dict[0], 'lightgrey'], splt=True, on_target=on_target_data, off_target=off_target_data)

    # here make updated output csv files that show both the reference species index as well as the actual species index
    # begin with the optimized sequences
    sorted_opti_seqs.insert(5, 'Reference index', sorted_opti_seqs['Indexes'])

    for dfidx, dfrow in sorted_opti_seqs.iterrows():
        referenceindexes = {}
        for species in dfrow['Reference index'].keys():
            referenceindexes[species] = conversion_dicts[species][dfrow['Reference index'][species]]
        sorted_opti_seqs['Reference index'][dfidx] = referenceindexes

    # do the same for raw data
    ranked_and_sorted_igs.insert(5, 'Reference index', ranked_and_sorted_igs['Index of Splice Site'])
    for dfidx, dfrow in ranked_and_sorted_igs.iterrows():
        referenceidx = conversion_dicts[dfrow['Target']][dfrow['Index of Splice Site']]
        ranked_and_sorted_igs['Reference index'][dfidx] = referenceidx

    sorted_opti_seqs.to_csv(file + '/Reference index Ranked Ribozyme Designs with Optimized Guide Sequence Designs.csv',
                            index=False)

    ranked_and_sorted_igs.to_csv(file + '/Reference index Ranked Ribozyme Designs with Raw Guide Sequence Designs.csv',
                                 index=False)

    # if we also want to do the generation data for NGS vs. RiboDesigner scores
    if riboscore_test:
        ngs_data = pd.read_csv(file + '/index counts from NGS data.csv', dtype=int)
        ribo_indexes = pd.Series(data=reference_indexes, name='E. coli indexes', dtype=int)
        ribo_scores = pd.Series(data=reference_scores, name='E. coli scores')

        new_data = pd.DataFrame(data=[ribo_indexes, ribo_scores], index=None).T

        joined_data = new_data.merge(ngs_data, 'right', left_on='E. coli indexes', right_on='Splice Site').sort_values(
            by=['Count', 'E. coli scores'], ascending=[False, False])
        joined_data.to_csv(file + '/NGS vs RiboDesigner scores.csv')

        joined_data['Color for plot'] = np.where((joined_data['Splice Site'] == ribozyme_to_test_index), 'r',
                                                 color_to_use)

        fig, ax = plt.subplots()

        ax.scatter(x=joined_data['Count'], y=joined_data['E. coli scores'], c=joined_data['Color for plot'])
        ax.xscale('log')
        ax.title('RiboDesigner scores vs. NGS output counts')
        ax.xlabel('Log counts')
        ax.ylabel('RiboDesigner Score')

        fig.savefig(file + '/riboDesigner scores.pdf', transparent=True)
    return


def set_params_for_plots(figsize, context):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': figsize}
    cmap = np.append([cmc.batlow.colors[0]], [cmc.batlow.colors[-1]], axis=0)
    sns.set_theme(context=context, style="ticks", rc=custom_params, palette=cmap)
    return


def score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path):
    # ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True]
    m = ribodesigner_settings[0]
    n = ribodesigner_settings[1]
    minlen = ribodesigner_settings[2]
    barcode_seq_file = ribodesigner_settings[3]
    ribobody_file = ribodesigner_settings[4]
    min_true_cov = ribodesigner_settings[5]
    identity_thresh = ribodesigner_settings[6]
    msa_fast = ribodesigner_settings[7]
    data = [None] * len(datasets)
    plt.figure()
    set_params_for_plots((30, 16), 'talk')

    fig, axs = plt.subplots(math.ceil(len(datasets) / 3), 3, sharex=True, sharey=True, layout="constrained")

    for i in range(0, len(datasets)):
        dataset_name = datasets[i].replace('_', ' ').replace('.fasta', '')

        output_path_folder = output_path + datasets[i].replace('.fasta', '') + '_results'
        target_sequences_folder = datasets_path + datasets[i]
        org_nums = len(read_fasta_folder(target_sequences_folder))

        if os.path.isdir(output_path_folder) is False or len(os.listdir(output_path_folder)) == 0:
            if os.path.isdir(output_path_folder) is False:
                os.mkdir(output_path_folder)

            out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, target_sequences_folder,
                                    min_true_cov=0, identity_thresh=0.7, ref_sequence_file=ref_path,
                                    fileout=True, folder_to_save=output_path_folder, msa_fast=msa_fast)
            out_data_df = pd.DataFrame(data=out_data, index=None, columns=['IGS', 'Reference index', 'Score', '% cov',
                                                                           '% on target', 'True % cov',
                                                                           '(Target name, Target idx, Other occurrences of'
                                                                           ' IGS in target sequence)',
                                                                           'Optimized guide',
                                                                           'Optimized guide + G + IGS',
                                                                           'Full Ribozyme design'],
                                       dtype=object).sort_values(by=['True % cov', 'Score'], ascending=[False, False])

        else:

            output_path_file = output_path_folder + '/Ranked Ribozyme Designs with Optimized Guide Sequence Designs quantitative.csv'
            out_data_df = pd.read_csv(output_path_file)

        # append data for heatmap
        try:
            big_df_for_heatmap = pd.concat([big_df_for_heatmap, out_data_df.loc[:, ['Score', 'True % cov']]],
                                           ignore_index=True)

        except:
            big_df_for_heatmap = out_data_df.loc[:, ['Score', 'Reference index', 'True % cov']]

        if '16s' in dataset_name:
            # Also make a separate dataframe for bacteria and archaea for the heatmap
            try:
                big_df_for_heatmap_bac = pd.concat([big_df_for_heatmap_bac,
                                                    out_data_df.loc[:, ['Score', 'Reference index', 'True % cov']]],
                                                   ignore_index=True)
            except:
                big_df_for_heatmap_bac = out_data_df.loc[:, ['Score', 'Reference index', 'True % cov']]

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
    plt.savefig(output_path + 'multi panel figure.png', transparent=False)
    plt.show()

    # # Now make a final plot comparing number of species vs. percentage above threshold
    set_params_for_plots((12, 8), 'talk')
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
    return


def plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path):
    # ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True]
    m = ribodesigner_settings[0]
    n = ribodesigner_settings[1]
    minlen = ribodesigner_settings[2]
    barcode_seq_file = ribodesigner_settings[3]
    ribobody_file = ribodesigner_settings[4]
    min_true_cov = ribodesigner_settings[5]
    identity_thresh = ribodesigner_settings[6]
    msa_fast = ribodesigner_settings[7]

    # variable regions V1-V9 (start, end) 1-based indexing on E. coli:
    var_regs = [(68, 100), (137, 226), (440, 496), (590, 650), (829, 856), (999, 1037), (1119, 1156), (1243, 1295),
                (1435, 1465)]

    # set up plots to look pretty
    set_params_for_plots((12, 8), 'talk')

    for i in range(0, len(datasets)):
        # custom_params = {"axes.spines.right": False, "axes.spines.top": False}
        # sns.set_theme(style="ticks", rc=custom_params, palette='viridis')
        # sns.set_context("paper")

        dataset_name = datasets[i].replace('_', ' ').replace('.fasta', '')
        # if '16s' not in dataset_name:
        #     # only align bacteria and archaea to E. coli
        #     continue

        output_path_folder = output_path + datasets[i].replace('.fasta', '') + '_ecoli_ref_results'
        target_sequences_folder = datasets_path + datasets[i]
        org_nums = len(read_fasta_folder(target_sequences_folder))

        if os.path.isdir(output_path_folder) is False or len(os.listdir(output_path_folder)) == 0:
            if os.path.isdir(output_path_folder) is False:
                os.mkdir(output_path_folder)

            # Same as before, but we have a reference sequence now (E. coli) to plot variable regions
            out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, target_sequences_folder,
                                ref_sequence_file=ref_path, min_true_cov=0, identity_thresh=0.7, fileout=True,
                                folder_to_save=output_path_folder, msa_fast=msa_fast)

            out_data_df = pd.DataFrame(data=out_data, index=None, columns=['IGS', 'Reference index', 'Score', '% cov',
                                                                           '% on target', 'True % cov',
                                                                           '(Target name, Target idx, Other occurrences of'
                                                                           ' IGS in target sequence)',
                                                                           'Optimized guide',
                                                                           'Optimized guide + G + IGS',
                                                                           'Full Ribozyme design'],
                                       dtype=object).sort_values(by=['True % cov', 'Score'], ascending=[False, False])


        else:
            output_path_file = output_path_folder + '/Ranked Ribozyme Designs with Optimized Guide Sequence ' \
                                                    'Designs quantitative.csv'
            out_data_df = pd.read_csv(output_path_file)
            out_data = out_data_df.values.tolist()

        colors = [None] * len(out_data_df.loc[:, ['Score']])

        for j in range(0, len(colors)):
            row = out_data_df.loc[j, ['True % cov', 'Score']]
            if row[0] > 0.7 and row[1] > 0.7:
                colors[j] = True
            else:
                colors[j] = False

        fig, ax = plt.subplots()
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 1600)
        plot_variable_regions(ax, var_regs)
        # plot score vs. site of U (ref--> E.coli)
        sns.scatterplot(data=out_data_df, x='Reference index', y='Score', alpha=0.8, hue=colors, legend=False)

        # Add labels
        ax.set_ylabel('Score')
        ax.set_xlabel('Location of U on E. coli 16s rRNA (base pair location)')
        ax.set_title(dataset_name + '\nScore along 16s for generated designs')

        # make pretty and save
        sns.despine(offset=10, trim=False)
        fig.savefig(output_path_folder + '/designs_along_16s.png', transparent=False)
        plt.show()
        return

