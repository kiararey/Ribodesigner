import figure_plotting_functions as fig_plt
import os
import numpy as np
import pandas as pd
from RiboDesigner import RiboDesigner
from v2RiboDesigner import RiboDesigner as oldRiboDesigner
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

    ########################################################
    # # Test datasets
    # print('Running test data...\n')

    # test_data = "/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/My Drive/KRG Thesis/Scripts/" \
    #             "Data files and Outputs/Ribozyme paper dataset/Original files"
    # print('########################################################\nNew multiprocessing - Super5:\n')
    # out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    #
    # print('########################################################\nOld multiprocessing - Super5:\n')
    # out_data_old = oldRiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data,
    #                                 min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                                     folder_to_save=output_path, msa_fast=True)
    #
    # test_data_archaea = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_1.fasta'
    #
    # print('########################################################\nNew multiprocessing - Super5:\n')
    # out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data_archaea,
    #                         min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                         folder_to_save=output_path, msa_fast=True)
    #
    # print('########################################################\nOld multiprocessing - Super5:\n')
    # out_data_old = oldRiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data_archaea,
    #                         min_true_cov=0, identity_thresh=0.7, fileout=True, ref_sequence_file=ref_path,
    #                         folder_to_save=output_path, msa_fast=True)
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
    dataset_names = ['Archaea_Only', 'Eukaryota_Only', 'Bacteria_Only']
    output_path = 'SILVA_output_files_Super5/E_coli_ref'
    all_above_coverage = []

    for name in dataset_names:
        datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
        print(f'Now analyzing data in {datasets_path[:-1]}...')

        datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
        datasets.sort()
        ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True, True]

        above_coverage = fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path,
                                       name.replace('_', ' '), 'png')
        # fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path,
        #                               name.replace('_', ' '), above_coverage, 'png')
        playsound('/System/Library/Sounds/Pop.aiff')

        all_above_coverage.append(above_coverage)

        # fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path,
        #                               name.replace('_', ' '))


    # Now plot the 16s rRNA figures one on top of the other
    fig_plt.plot_for_16s_coverage_multipanel(all_above_coverage, dataset_names, output_path, file_type='png')

    # Now, align to model organisms and make more data!
    archea_model_ref_path = 'Common_sequences/Methanobrevibacter smithii 16s.fasta'
    eukaryota_model_ref_path = 'Common_sequences/Saccharomyces cerevisiae 18s.fasta'
    archea_output_path = 'SILVA_output_files_Super5/M_smithii_ref'
    eukaryota_output_path = 'SILVA_output_files_Super5/S_cerevisiae_ref'

    dataset_names = ['Archaea_Only', 'Eukaryota_Only']
    ref_paths = [archea_model_ref_path, eukaryota_model_ref_path]
    output_paths = [archea_output_path, eukaryota_output_path]

    for i, name in enumerate(dataset_names):
        datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
        print(f'Now analyzing data in {datasets_path[:-1]}...')

        datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
        datasets.sort()
        ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True, True]

        above_coverage = fig_plt.score_vs_true_coverage(datasets, datasets_path, output_paths[i], ribodesigner_settings,
                                                        ref_paths[i],
                                                        name.replace('_', ' '), 'png')
        fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_paths[i], ribodesigner_settings, ref_paths[i],
                                      name.replace('_', ' '), above_coverage, 'png')
        playsound('/System/Library/Sounds/Pop.aiff')

