import figure_plotting_functions as fig_plt
import os
import numpy as np
import pandas as pd
from RiboDesigner import RiboDesigner
from oldribodesigner import RiboDesigner as oldRiboDesigner
from playsound import playsound
from ribodesigner_v2 import ribodesigner

if __name__ == '__main__':
    # Figure 2: synthetic community data violin plots
    # Basically, look at the datasets used data and then run RiboDesigner on them.
    # Generate a violin plot using the results.
    # Run RiboDesigner on all datasets we are looking at
    # First, set up base files and parameters
    m = 5
    n = 50
    minlen = 50

    # Barcode sequence is split sfGFP just cuz. This does not affect guide sequence design.
    barcode_seq_file = 'Common_sequences/sfGFP_2_seq_barcode.txt'

    # We'll be using the Ribozyme published in the RAM paper
    ribobody_file = 'Common_sequences/ribozyme_body.txt'

    # Prepare the datasets
    datasets_path = 'Datasets_used/zymo_files/'

    # Output folder
    output_path = 'test_output_files/'

    # Reference sequence - will have to be E. coli to graph the variable regions
    ref_path = 'Common_sequences/e-coli-16s-mg1655.fasta'

    # ########################################################
    # # Test figure making
    # print('Testing figure making...')
    # test_data = 'Datasets_used/Test_data/'
    #
    # # ribodesigner_settings = [barcode_seq_file, ribobody_file, m, n, minlen, min_true_cov, identity_thresh, msa_fast]
    # ribodesigner_settings = [barcode_seq_file, ribobody_file, m, n, minlen, 0, 0.7, True]
    # fig_plt.score_vs_true_coverage(['test_fasta_for_figures.fasta'], test_data, output_path, ribodesigner_settings, ref_path,
    #                                'test data')
    #
    # fig_plt.plot_for_16s_coverage(['test_fasta_for_figures.fasta'], test_data, output_path, ribodesigner_settings, ref_path,
    #                               'test data')
    # print(f'Figure test done!\n########################################################\n')

    # ########################################################
    # # Test datasets
    # print('Running test data...\n')
    #
    # test_data = "/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/My Drive/KRG Thesis/Scripts/" \
    #             "Data files and Outputs/Ribozyme paper dataset/Original files"
    # print('########################################################\nNew less memory - Super5:\n')
    # out_data = RiboDesigner(test_data, barcode_seq_file, ribobody_file, m, n, minlen,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=False, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    #
    # print('########################################################\nOld more momory - Super5:\n')
    # out_data_old = oldRiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=False, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    # print(f'Are old data and new data the same? {out_data==out_data_old}')
    #
    # test_data_archaea = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_1.fasta'
    #
    # print('########################################################\nNew less memory - Super5:\n')
    # out_data = RiboDesigner(test_data, barcode_seq_file, ribobody_file, m, n, minlen,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=False, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    #
    # print('########################################################\nOld more momory - Super5:\n')
    # out_data_old = oldRiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data_archaea,
    #                                min_true_cov=0, identity_thresh=0.7, fileout=False, ref_sequence_file=ref_path,
    #                                folder_to_save=output_path, msa_fast=True)
    #
    # print(f'Are old data and new data the same? {out_data == out_data_old}')
    # playsound('/System/Library/Sounds/Pop.aiff')
    # print(f'Test data done!\n########################################################\n')

    ########################################################
    # test data targeted
    good_targets = "/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/My Drive/KRG Thesis/Scripts/" \
                   "Data files and Outputs/Ribozyme paper dataset/Original files"
    bad_targets = 'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Species/Bacillus_halotolerans.fasta'
    big_data = f'Datasets_used/SILVA_squished_datasets/Enterobacterales_only_squished'

    # Test new RiboDesigner
    out_data_weighted = ribodesigner(target_sequences_folder=bad_targets, barcode_seq_file=barcode_seq_file,
                            ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
                            background_sequences_folder=bad_targets, min_true_cov=0.7, identity_thresh=0.7,
                            fileout=False, msa_fast=True, ref_sequence_file=ref_path, gaps_allowed=False,
                            percent_of_background_seqs_used=0.75, score_type='naive', n_limit=0)

    print('########################################################\n')

    out_data_weighted = RiboDesigner(target_sequences_folder=good_targets, barcode_seq_file=barcode_seq_file,
                            ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
                            background_sequences_folder=bad_targets, min_true_cov=0.7, identity_thresh=0.7,
                            fileout=False, msa_fast=True, ref_sequence_file=ref_path, gaps_allowed=False,
                            percent_of_background_seqs_used=0.75, score_type='naive', n_limit=0)

    out_data_naive = RiboDesigner(target_sequences_folder=good_targets, barcode_seq_file=barcode_seq_file,
                            ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
                            background_sequences_folder=bad_targets, min_true_cov=0.7, identity_thresh=0.7,
                            fileout=False, msa_fast=True, ref_sequence_file=ref_path, gaps_allowed=False,
                            percent_of_background_seqs_used=0.75, score_type='weighted')

    # out_data = RiboDesigner(target_sequences_folder=good_targets, barcode_seq_file=barcode_seq_file,
    #                         ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
    #                         background_sequences_folder=bad_targets, min_true_cov=0.7, identity_thresh=0.7,
    #                         fileout=False, msa_fast=True, ref_sequence_file=ref_path, gaps_allowed=False,
    #                         percent_of_background_seqs_used=0.75, score_type='quantitative')

    playsound('/System/Library/Sounds/Pop.aiff')
    print(f'Test data done!\n########################################################\n')

    # ########################################################
    # # Score vs. True Coverage graphs
    # # Generate data here. If already generated we'll just import it with pandas
    # # For each folder in datasets used run RiboDesigner and keep the data for plotting later
    # # get rid of that pesky .DS_Store file
    # datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
    # datasets.sort()
    # ribodesigner_settings = [barcode_seq_file, ribobody_file, m, n, minlen, min_true_cov, identity_thresh, msa_fast]
    # ribodesigner_settings = [barcode_seq_file, ribobody_file, m, n, minlen, 0, 0.7, True]
    # fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, None)
    # fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)

    # ########################################################
    # # SILVA squished datasets
    # dataset_names = ['Archaea_Only', 'Eukaryota_Only', 'Bacteria_Only', 'All_Kingdoms']
    # # dataset_names = ['All_Kingdoms']
    # output_path = 'SILVA_output_files_Super5/E_coli_ref/'
    # all_above_coverage = []
    #
    # for name in dataset_names:
    #     datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
    #     print(f'Now analyzing data in {datasets_path[:-1]}...')
    #
    #     datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
    #     datasets.sort()
    #     ribodesigner_settings = [barcode_seq_file, ribobody_file, m, n, minlen, 0, 0.7, True]
    #
    #     above_coverage = fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings,
    #                                                     ref_path, name.replace('_', ' '), 'png')
    #     fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path,
    #                                   name.replace('_', ' '), above_coverage, 'png')
    #
    #     all_above_coverage.append(above_coverage)
    #
    #
    # # Now plot the 16s rRNA figures one on top of the other
    # # The largest SSU is the 18s, we need this to keep everything on the same axis!
    # limit_override = 1800
    # fig_plt.plot_for_16s_coverage_multipanel(all_above_coverage, dataset_names, output_path, file_type='png',
    #                                          xlim=limit_override)
    #
    # # Now, align to model organisms and make more data!
    # archea_model_ref_path = 'Common_sequences/Methanobrevibacter smithii 16s.fasta'
    # eukaryota_model_ref_path = 'Common_sequences/Saccharomyces cerevisiae 18s.fasta'
    # archea_output_path = 'SILVA_output_files_Super5/M_smithii_ref/'
    # eukaryota_output_path = 'SILVA_output_files_Super5/S_cerevisiae_ref/'
    #
    # dataset_names = ['Archaea_Only', 'Eukaryota_Only']
    # ref_paths = [archea_model_ref_path, eukaryota_model_ref_path]
    # output_paths = [archea_output_path, eukaryota_output_path]
    #
    # # The variable regions below were obtained by doing a mafft pairwise alignment with MG1655. Let me know if this is
    # # incorrect please and I will fix it! This is base-1 indexing btw.
    # FIX THIS!!!!!! need to align to new variable regions
    # m_smithii_var_regs = [(64, 74), (112, 205), (394, 433), (528, 589), (768, 799), (940, 978), (1056, 1103),
    #                       (1190, 1242), (1384, 1402)]
    # s_cerevisiae_var_regs = [(69, 80), (126, 292), (478, 510), (643, 850), (1048, 1070), (1350, 1400), (1480, 1531),
    #                          (1674, 1730)]  # this is from 1. Reich, M. & Labes, A. How to boost marine fungal research: A first step towards a multidisciplinary approach by combining molecular fungal ecology and natural products chemistry. Marine Genomics 36, 57–75 (2017).
    # var_regs_overrides = [m_smithii_var_regs, s_cerevisiae_var_regs, None]
    #
    # for i, name in enumerate(dataset_names):
    #     datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
    #     print(f'Now analyzing data in {datasets_path[:-1]}...')
    #
    #     datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
    #     datasets.sort()
    #     ribodesigner_settings = [barcode_seq_file, ribobody_file, m, n, minlen, 0, 0.7, True]
    #
    #     above_coverage = fig_plt.score_vs_true_coverage(datasets, datasets_path, output_paths[i], ribodesigner_settings,
    #                                                     ref_paths[i], name.replace('_', ' '), 'png')
    #     fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_paths[i], ribodesigner_settings, ref_paths[i],
    #                                   name.replace('_', ' '), above_coverage, 'png', var_regs_overrides[i])
    #
    #     all_above_coverage[i] = above_coverage
    #
    # dataset_names = ['Archaea_Only', 'Eukaryota_Only', 'Bacteria_Only']
    # big_out_path = 'SILVA_output_files_Super5/'
    # # Now plot the 16s rRNA figures one on top of the other
    # fig_plt.plot_for_16s_coverage_multipanel(all_above_coverage, dataset_names, big_out_path, file_type='png',
    #                                          var_regs_overrides=var_regs_overrides, share_x_lim=False)
    # playsound('/System/Library/Sounds/Pop.aiff')

    # ########################################################
    # # Targeted designs
    # good_targets = 'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_1.fasta'
    # bad_targets = 'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_1.fasta'
    # output_path = 'SILVA_output_files_Super5/Targeted'
    #
    # RiboDesigner(target_sequences_folder=good_targets, barcode_seq_file=barcode_seq_file, ribobody_file=ribobody_file,
    #              igs_length=m, guide_length=n, min_length=minlen, targeted=True,
    #              background_sequences_folder=bad_targets, min_true_cov=0.7, identity_thresh=0.7, fileout=True,
    #              ref_sequence_file=ref_path, folder_to_save=output_path, msa_fast=True)
    #
    # playsound('/System/Library/Sounds/Pop.aiff')

    ########################################################
    # # Targeted designs order level
    # enterobacterales = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Order/Enterobacterales.fasta'
    # pseudomonadales = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Order/Pseudomonadales.fasta'
    # output_path_e = 'SILVA_output_files_Super5/for_paper/targets_enterobacterales_0.2_per_threshold'
    # output_path_p = 'SILVA_output_files_Super5/for_paper/targets_pseudomondales_0.2_per_threshold'
    #
    # RiboDesigner(target_sequences_folder=enterobacterales, barcode_seq_file=barcode_seq_file,
    #              ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
    #              background_sequences_folder=pseudomonadales, min_true_cov=0.7, identity_thresh=0.7, fileout=True,
    #              ref_sequence_file=ref_path, folder_to_save=output_path_e, msa_fast=True, gaps_allowed=False,
    #              percent_of_target_seqs_used=.002, percent_of_background_seqs_used=.05, score_type='naive')
    # print('Dataset 1 complete.\n')
    #
    # RiboDesigner(target_sequences_folder=pseudomonadales, barcode_seq_file=barcode_seq_file,
    #              ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
    #              background_sequences_folder=enterobacterales, min_true_cov=0.7, identity_thresh=0.7, fileout=True,
    #              ref_sequence_file=ref_path, folder_to_save=output_path_p, msa_fast=True, gaps_allowed=False,
    #              percent_of_target_seqs_used=.002, percent_of_background_seqs_used=.05, score_type='naive')
    #
    # playsound('/System/Library/Sounds/Pop.aiff')
    #
    # output_path_e = 'SILVA_output_files_Super5/for_paper/targets_enterobacterales_weighted_0.2_per_threshold'
    # output_path_p = 'SILVA_output_files_Super5/for_paper/targets_pseudomondales_weighted_0.2_per_threshold'
    #
    # RiboDesigner(target_sequences_folder=enterobacterales, barcode_seq_file=barcode_seq_file,
    #              ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
    #              background_sequences_folder=pseudomonadales, min_true_cov=0.7, identity_thresh=0.7, fileout=True,
    #              ref_sequence_file=ref_path, folder_to_save=output_path_e, msa_fast=True, gaps_allowed=False,
    #              percent_of_target_seqs_used=.002, percent_of_background_seqs_used=.05, score_type='weighted')
    # print('Dataset 1 compplete.\n')
    #
    # RiboDesigner(target_sequences_folder=pseudomonadales, barcode_seq_file=barcode_seq_file,
    #              ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen, targeted=True,
    #              background_sequences_folder=enterobacterales, min_true_cov=0.7, identity_thresh=0.7, fileout=True,
    #              ref_sequence_file=ref_path, folder_to_save=output_path_p, msa_fast=True, gaps_allowed=False,
    #              percent_of_target_seqs_used=.002, percent_of_background_seqs_used=.05, score_type='weighted')
    #
    # playsound('/System/Library/Sounds/Pop.aiff')

    ########################################################
    # Targeted designs by order level: run many times and simulate
    output_path = 'SILVA_output_files_Super5/Targeted'
    enterobacterales = f'Datasets_used/SILVA_squished_datasets/Enterobacterales_only_squished'
    pseudomonadales = f'Datasets_used/SILVA_squished_datasets/Pseudomonadales_only_squished'

    background_seqs_no_entero = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished/' \
                                'Background_Bacteria_squished_no_entero.fasta'
    background_seqs_no_pseudo = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished/' \
                                'Background_Bacteria_squished_no_pseudo.fasta'

    # enterobacterales = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Order/Enterobacterales.fasta'
    # pseudomonadales = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Order/Pseudomonadales.fasta'

    e_datasets = np.array([file_name for file_name in os.listdir(enterobacterales) if '.fasta' in file_name])
    p_datasets = np.array([file_name for file_name in os.listdir(pseudomonadales) if '.fasta' in file_name])

    e_designs = []
    e_by_idx = {}
    p_designs = []
    p_by_idx = {}

    i = 1
    for e_set, p_set in zip(e_datasets, p_datasets):
        e_out_folder = f'{output_path}/Enterobacterales_only_squished_{i}'
        e_file = f'{enterobacterales}/{e_set}'
        p_out_folder = f'{output_path}/Pseudomonadales_only_squished_{i}'
        p_file = f'{pseudomonadales}/{p_set}'
        i += 1

        e_out = RiboDesigner(target_sequences_folder=e_file, barcode_seq_file=barcode_seq_file,
                             ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen,
                             targeted=True, background_sequences_folder=background_seqs_no_entero, min_true_cov=0.7,
                             identity_thresh=0.7, fileout=True, folder_to_save=e_out_folder, ref_sequence_file=ref_path,
                             msa_fast=True, gaps_allowed=False, min_delta=0, score_type='weighted', n_limit=2.0/50)
        e_designs.append(e_out)


        for design in e_out:
            # each design is saved as: 'IGS, Reference index, Score, % cov, % on target, True % cov, Composite score,
            # Adjusted score vs. background, Number of species targeted,
            # Optimized guide, Optimized guide + G + IGS, Full Ribozyme design, Delta composite score vs background
            id_code = f'{design[0]}{design[1]}'
            # Save for each design: id, Adjusted score vs. background, Delta composite score vs background, optimized guide
            design_info = (id_code, design[7], design[8], design[10])

            # Now separate by ref id
            try:
                to_extend = e_by_idx[design[1]]
                e_by_idx[design[1]] = to_extend.extend(design_info)
            except:
                e_by_idx[design[1]] = [design_info]


        p_out = RiboDesigner(target_sequences_folder=p_file, barcode_seq_file=barcode_seq_file,
                             ribobody_file=ribobody_file, igs_length=m, guide_length=n, min_length=minlen,
                             targeted=True, background_sequences_folder=background_seqs_no_pseudo, min_true_cov=0.7,
                             identity_thresh=0.7, fileout=True, folder_to_save=p_out_folder, ref_sequence_file=ref_path,
                             msa_fast=True, gaps_allowed=False, min_delta=0, score_type='weighted', n_limit=2.0/50)
        p_designs.append(p_out)

        for design in p_out:
            # each design is saved as: 'IGS, Reference index, Score, % cov, % on target, True % cov, Composite score,
            # Adjusted score vs. background, Number of species targeted,
            # Optimized guide, Optimized guide + G + IGS, Full Ribozyme design, Delta composite score vs background
            id_code = f'{design[0]}{design[1]}'
            # Save for each design: id, Adjusted score vs. background, Delta composite score vs background, optimized guide
            design_info = (id_code, design[7], design[8], design[10])

            # Now separate by ref id
            try:
                to_extend = p_by_idx[design[1]]
                p_by_idx[design[1]] = to_extend.extend(design_info)
            except:
                p_by_idx[design[1]] = [design_info]


    # Then, check results: order by reference index


    # Graph counts per index on 16s rRNA

    # Pick designs that happened the most frequently. Can you make a consensus sequence? (may make this worse but idk)


    # Now, run these designs against all background sequences. Check how they do

    playsound('/System/Library/Sounds/Pop.aiff')


# then test in silico
# 2.4803 sec