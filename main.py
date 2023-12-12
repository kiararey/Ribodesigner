import sys

from ribodesigner import (ribodesigner, ribo_checker, couple_designs_to_test_seqs, prepare_test_seqs)
from graph_making import make_graphs, make_sequence_logo_graph

if __name__ == '__main__':
    # Run RiboDesigner on all datasets we are looking at
    # First, set up base files and parameters
    m = 5
    n = 50
    minlen = 35

    try:
        worker_number = sys.argv[1]
        number_of_workers = sys.argv[2]
    except:
        worker_number = 0
        number_of_workers = 50
        pass

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

    # Our random sequence!
    random_seq_path = 'Common_sequences/Lactobacillus_casei_example.fasta'

    u1376 = 'CAACCCACTCCCATGGTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGTgTTCAC'

    # this is from Reich, M. & Labes, A. How to boost marine fungal research: A first step towards a multidisciplinary
    # approach by combining molecular fungal ecology and natural products chemistry. Marine Genomics 36, 57–75 (2017).
    s_cerevisiae_var_regs = [(69, 80), (126, 292), (478, 510), (643, 850), (1048, 1070), (1350, 1400), (1480, 1531),
                             (1674, 1730)]

    # Chakravorty, S., Helb, D., Burday, M., Connell, N. & Alland, D. A detailed analysis of 16S ribosomal RNA gene
    # segments for the diagnosis of pathogenic bacteria. J Microbiol Methods 69, 330–339 (2007).
    e_coli_var_regs = [(69, 99), (137, 242), (433, 497), (576, 682), (822, 879), (986, 1043), (1117, 1173),
                       (1243, 1294), (1435, 1465)]

    # ########################################################
    # test data targeted
    path = 'Datasets_used/SILVA_squished_datasets_1_per_genus/'
    bad_targets = 'Datasets_used/Bacillus_halotolerans.fasta'
    universal_data_1 = path + 'SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_2.fasta'
    universal_data_2 = path + 'SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_2.fasta'
    universal_data_3 = path + 'SILVA_squished_datasets_Eukaryota_Only/Eukaryota_Only_by_Genus_2.fasta'
    universal_data_4 = path + 'SILVA_squished_datasets_All_Kingdoms/All_by_Genus_2.fasta'
    big_data_entero_only = path + 'Enterobacterales_only_squished/Enterobacterales_only_by_Genus_1.fasta'
    big_data_pseudo_only = path + 'Pseudomonadales_only_squished/Pseudomonadales_only_by_Genus_1.fasta'
    big_data_no_entero_or_pseudo = path + 'Background_Bacteria_squished/' \
                                          'Background_Bacteria_squished_no_pseudo_or_entero.fasta'
    big_data_no_pseudo = path + 'Background_Bacteria_squished/Background_Bacteria_squished_no_pseudo.fasta'
    big_data_no_entero = path + 'Background_Bacteria_squished/Background_Bacteria_squished_no_entero.fasta'
    big_data_only_entero_and_pseudo = path + 'Pseudo_and_entero_only_squished/Pseudo_and_entero_only_by_Genus_1.fasta'
    background_data_bac = path + 'SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_1.fasta'
    background_data_arc = path + 'SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_1.fasta'
    background_data_euk = path + 'SILVA_squished_datasets_Eukaryota_Only/Eukaryota_Only_by_Genus_1.fasta'
    background_data_all = path + 'SILVA_squished_datasets_All_Kingdoms/All_by_Genus_1.fasta'
    test_output_folder = 'test_output_files/test_outputs_parallelizing'
    test_file = 'test_dataset_for_graphs.csv'
    big_data_file_for_output = 'large_dataset.csv'
    ref_analysis_folder = 'test_output_files/test_outputs_parallelizing/native ecoli mg1655 designs'

    test_data_folders = [background_data_euk, background_data_arc, background_data_bac, background_data_all]
    universal_data_files = [universal_data_1, universal_data_2, universal_data_3, universal_data_4]
    test_data_folders_test = [bad_targets, big_data_entero_only]
    # Test new RiboDesigner for images
    universal_datasets = []
    selective_datasets = []

    test_data_pickles = ['test_sequences_All_by_Genus_1.pickle', 'test_sequences_Archaea_Only_by_Genus_1.pickle',
                         'test_sequences_Bacillus_halotolerans.pickle',
                         'test_sequences_Bacteria_Only_by_Genus_1.pickle',
                         'test_sequences_Eukaryota_Only_by_Genus_1.pickle']
    test_data_pickles = [test_output_folder + '/' + item for item in test_data_pickles]

    universal_data_pickles = ['designs_Eukaryota_Only_by_Genus_2_universal.pickle',
                              'designs_All_by_Genus_2_universal.pickle',
                              'designs_Archaea_Only_by_Genus_2_universal.pickle',
                              'designs_Bacteria_Only_by_Genus_2_universal.pickle']
    universal_data_pickles = [test_output_folder + '/' + item for item in universal_data_pickles]

    ref_seq_pickle = test_output_folder + '/designs_e-coli-16s-mg1655_universal.pickle'
    random_seq_pickle = test_output_folder + '/designs_Lactobacillus_casei_example_universal.pickle'

    # # Here we make the designs
    # for test_data in test_data_folders:
    #     test_seqs_pickle_file_name = prepare_test_seqs(test_folder=test_data, ref_sequence_file=ref_path,
    #                                                    guide_length=n, igs_length=m, min_length=minlen,
    #                                                    folder_to_save=test_output_folder, graph_results=True,
    #                                                    var_regs=e_coli_var_regs, graph_file_type='png')
    #     test_data_pickles.append(test_seqs_pickle_file_name)
    #
    # # Here, we're using ribodesigner functions to see what would happen if we used the native sequences after each
    # # U site as guides in E. coli MG1655
    # ref_seq_pickle_file_name = ribodesigner(target_sequences_folder=ref_path,
    #                                         ref_sequence_file=ref_path, igs_length=m,
    #                                         guide_length=n, min_length=minlen, selective=False, min_true_cov=0,
    #                                         msa_fast=True,
    #                                         score_type='weighted', n_limit=1, percent_of_target_seqs_used=1,
    #                                         gaps_allowed=False, fileout=False, random_guide_sample_size=10,
    #                                         test_folders=test_data_folders, folder_to_save=test_output_folder)
    # Similarly, what if we chose another sequence
    # random_seq_pickle_file_name = ribodesigner(target_sequences_folder=random_seq_path,
    #                                         ref_sequence_file=ref_path, igs_length=m,
    #                                         guide_length=n, min_length=minlen, selective=False, min_true_cov=0,
    #                                         msa_fast=True,
    #                                         score_type='weighted', n_limit=1, percent_of_target_seqs_used=1,
    #                                         gaps_allowed=False, fileout=False, random_guide_sample_size=10,
    #                                         test_folders=test_data_folders, folder_to_save=test_output_folder)
    # # finally, here is our universal designs
    # for universal_data in universal_data_files:
    #     universal_pickle_file_name = ribodesigner(target_sequences_folder=universal_data,
    #                                             ref_sequence_file=ref_path, igs_length=m,
    #                                             guide_length=n, min_length=minlen, selective=False, min_true_cov=0,
    #                                             msa_fast=True,
    #                                             score_type='weighted', n_limit=1, percent_of_target_seqs_used=1,
    #                                             gaps_allowed=False, fileout=False, random_guide_sample_size=10,
    #                                             test_folders=test_data_folders, folder_to_save=test_output_folder)
    #     universal_data_pickles.append(universal_pickle_file_name)

    # # Now we couple the designs with their test sequences to later test them
    # for test_pickle in test_data_pickles:
    #     random_seq_coupled_file_name = couple_designs_to_test_seqs(designs_input=random_seq_pickle,
    #                                                             test_seqs_input=test_pickle,
    #                                                             flexible_igs=True, igs_len=m, ref_idx_u_loc=1376,
    #                                                             score_type='weighted',
    #                                                             file_to_save=test_output_folder)
    #
    #     ref_seq_coupled_file_name = couple_designs_to_test_seqs(designs_input=ref_seq_pickle,
    #                                                            test_seqs_input=test_pickle,
    #                                                            flexible_igs=True, igs_len=m, ref_idx_u_loc=1376,
    #                                                            score_type='weighted',
    #                                                            file_to_save=test_output_folder)
    #
    #     control_design_coupled_file_name = couple_designs_to_test_seqs(designs_input=u1376,
    #                                                                   test_seqs_input=test_pickle,
    #                                                                   flexible_igs=True, igs_len=m, ref_idx_u_loc=1376,
    #                                                                   score_type='weighted',
    #                                                                   file_to_save=test_output_folder)
    #
    #     for universal_design_pickle in universal_data_pickles:
    #         universal_designs_coupled_file_name = couple_designs_to_test_seqs(designs_input=universal_design_pickle,
    #                                                                        test_seqs_input=test_pickle,
    #                                                                        flexible_igs=True,
    #                                                                        file_to_save=test_output_folder)

    # # Below we're making coupled files for testing purposes =
    # test_pickle = test_output_folder + '/test_sequences_Bacteria_Only_by_Genus_1.pickle'
    # ref_seq_coupled_file_name = couple_designs_to_test_seqs(designs_input=ref_seq_pickle, test_seqs_input=test_pickle,
    #                                                         flexible_igs=True, igs_len=m, ref_idx_u_loc=1376,
    #                                                         score_type='weighted', file_to_save=test_output_folder)
    #
    # control_design_coupled_file_name = couple_designs_to_test_seqs(designs_input=u1376, test_seqs_input=test_pickle,
    #                                                                flexible_igs=True, igs_len=m, ref_idx_u_loc=1376,
    #                                                                score_type='weighted',
    #                                                                file_to_save=test_output_folder)
    #
    # universal_designs_pickle = test_output_folder + '/designs_Bacteria_Only_by_Genus_2_universal.pickle'
    # universal_designs_coupled_file_name = couple_designs_to_test_seqs(designs_input=universal_designs_pickle,
    #                                                                   test_seqs_input=test_pickle, flexible_igs=True,
    #                                                                   file_to_save=test_output_folder)

    # # finally, we test! Below is for local
    files_to_test = test_output_folder + '/coupled'
    for i in range(number_of_workers):
        ribo_checker(coupled_folder=files_to_test, number_of_workers=number_of_workers, worker_number=0,
                     n_limit=1)
    # # This is for NOTS
    # files_to_test = test_output_folder + '/coupled'
    # ribo_checker(coupled_folder='/scratch/kpr1/RiboDesigner/' + files_to_test, number_of_workers=number_of_workers,
    #              worker_number=worker_number, n_limit=1)

    # print(f'Test data done!\n########################################################\n')

    control_design_results_file_name = test_output_folder + '/Coupled for testing/for_testing/results/1_designs_TTCAC1376_vs_test_sequences_Bacteria_Only_by_Genus_1_worker_0_results.txt'
    universal_designs_results_file_name = test_output_folder + '/Coupled for testing/for_testing/results/77193_designs_designs_Bacteria_Only_by_Genus_2_universal_vs_test_sequences_Bacteria_Only_by_Genus_1_worker_0_results.txt'
    ref_seq_results_file_name = test_output_folder + '/Coupled for testing/for_testing/results/304_designs_designs_e-coli-16s-mg1655_universal_vs_test_sequences_Bacteria_Only_by_Genus_1_worker_0_results.txt'
    random_seq_results_file_name = test_output_folder + '/Coupled for testing/for_testing/results/326_designs_designs_Lactobacillus_casei_example_universal_vs_test_sequences_Bacteria_Only_by_Genus_1_worker_0_results.txt'

    folder_for_ref_seq_results = 'test_output_files/test_outputs_parallelizing/Coupled for testing/for_testing/figures/ref_seq is MG1655'
    folder_for_random_seq_results = 'test_output_files/test_outputs_parallelizing/Coupled for testing/for_testing/figures/ref_seq is L shirota'
    # make_graphs(control_designs_path=control_design_results_file_name,
    #             universal_designs_path=universal_designs_results_file_name,
    #             ref_seq_designs_path=ref_seq_results_file_name, var_regs=e_coli_var_regs, save_fig=False,
    #             save_file_loc='test_output_files/test_outputs_parallelizing/coupled/figures')
    #
    # ########################################################
    make_graphs(control_designs_path=control_design_results_file_name,
                universal_designs_path=universal_designs_results_file_name,
                ref_seq_designs_path=ref_seq_results_file_name, var_regs=e_coli_var_regs,
                save_fig=True, save_file_loc=folder_for_ref_seq_results)

    make_graphs(control_designs_path=control_design_results_file_name,
                universal_designs_path=universal_designs_results_file_name,
                ref_seq_designs_path=random_seq_results_file_name, var_regs=e_coli_var_regs, save_fig=True,
                save_file_loc=folder_for_random_seq_results)

    make_sequence_logo_graph(test_data_path=test_data_pickles[3], design_data_path=universal_designs_results_file_name,
                             ref_data_path=ref_seq_results_file_name)

    make_sequence_logo_graph(test_data_path=test_data_pickles[3], design_data_path=universal_designs_results_file_name,
                             ref_data_path=random_seq_results_file_name)

    print(f'Graphs done!\n########################################################\n')

