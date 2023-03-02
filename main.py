import figure_plotting_functions as fig_plt
import os
import numpy as np
import pandas as pd
from RiboDesigner import RiboDesigner
from oldribodesigner import RiboDesigner as oldRiboDesigner
from playsound import playsound

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
    # ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True, True]
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
    # out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    #
    # print('########################################################\nOld more momory - Super5:\n')
    # out_data_old = oldRiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    # print(f'Are old data and new data the same? {out_data==out_data_old}')
    #
    # test_data_archaea = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_1.fasta'
    #
    # print('########################################################\nNew less memory - Super5:\n')
    # out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data_archaea,
    #                         min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                         folder_to_save=output_path, msa_fast=True)
    #
    # print('########################################################\nOld more momory - Super5:\n')
    # out_data_old = oldRiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data_archaea,
    #                                min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                                folder_to_save=output_path, msa_fast=True)
    #
    # print(f'Are old data and new data the same? {out_data == out_data_old}')
    # playsound('/System/Library/Sounds/Pop.aiff')
    # print(f'Test data done!\n########################################################\n')

    # ########################################################
    # # Score vs. True Coverage graphs
    # # Generate data here. If already generated we'll just import it with pandas
    # # For each folder in datasets used run RiboDesigner and keep the data for plotting later
    # # get rid of that pesky .DS_Store file
    # datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
    # datasets.sort()
    # ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True]
    # fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, None)
    # fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)

    # ########################################################
    # SILVA squished datasets
    dataset_names = ['Archaea_Only', 'Eukaryota_Only', 'Bacteria_Only', 'All']
    dataset_names = ['All']
    output_path = 'SILVA_output_files_Super5/E_coli_ref/'
    all_above_coverage = []

    for name in dataset_names:
        datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
        print(f'Now analyzing data in {datasets_path[:-1]}...')

        datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
        datasets.sort()
        ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0.7, 0.7, True, True]

        above_coverage = fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings,
                                                        ref_path, name.replace('_', ' '), 'png')
        fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path,
                                      name.replace('_', ' '), above_coverage, 'png')

        all_above_coverage.append(above_coverage)


    # Now plot the 16s rRNA figures one on top of the other
    # The largest SSU is the 18s, we need this to keep everything on the same axis!
    limit_override = 1800
    fig_plt.plot_for_16s_coverage_multipanel(all_above_coverage, dataset_names, output_path, file_type='png',
                                             xlim=limit_override)

    # Now, align to model organisms and make more data!
    archea_model_ref_path = 'Common_sequences/Methanobrevibacter smithii 16s.fasta'
    eukaryota_model_ref_path = 'Common_sequences/Saccharomyces cerevisiae 18s.fasta'
    archea_output_path = 'SILVA_output_files_Super5/M_smithii_ref/'
    eukaryota_output_path = 'SILVA_output_files_Super5/S_cerevisiae_ref/'

    dataset_names = ['Archaea_Only', 'Eukaryota_Only']
    ref_paths = [archea_model_ref_path, eukaryota_model_ref_path]
    output_paths = [archea_output_path, eukaryota_output_path]

    # The variable regions below were obtained by doing a mafft pairwise alignment with MG1655. Let me know if this is
    # incorrect please and I will fix it! This is base-1 indexing btw.
    m_smithii_var_regs = [(64, 74), (112, 205), (394, 433), (528, 589), (768, 799), (940, 978), (1056, 1103),
                          (1190, 1242), (1384, 1402)]
    s_cerevisiae_var_regs = [(69, 80), (126, 292), (478, 510), (643, 850), (1048, 1070), (1350, 1400), (1480, 1531),
                             (1674, 1730)]
    var_regs_overrides = [m_smithii_var_regs, s_cerevisiae_var_regs, None]

    for i, name in enumerate(dataset_names):
        datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
        print(f'Now analyzing data in {datasets_path[:-1]}...')

        datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
        datasets.sort()
        ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True, True]

        above_coverage = fig_plt.score_vs_true_coverage(datasets, datasets_path, output_paths[i], ribodesigner_settings,
                                                        ref_paths[i], name.replace('_', ' '), 'png')
        fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_paths[i], ribodesigner_settings, ref_paths[i],
                                      name.replace('_', ' '), above_coverage, 'png', var_regs_overrides[i])

        all_above_coverage[i] = above_coverage

    dataset_names = ['Archaea_Only', 'Eukaryota_Only', 'Bacteria_Only']
    big_out_path = 'SILVA_output_files_Super5/'
    # Now plot the 16s rRNA figures one on top of the other
    fig_plt.plot_for_16s_coverage_multipanel(all_above_coverage, dataset_names, big_out_path, file_type='png',
                                             var_regs_overrides=var_regs_overrides, share_x_lim=False)
    playsound('/System/Library/Sounds/Pop.aiff')

