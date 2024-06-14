import sys
import os
from alive_progress import alive_bar
import multiprocessing as mp
from ribodesigner import (ribodesigner, ribo_checker, couple_designs_to_test_seqs, prepare_test_seqs, combine_data,
                          select_designs, adjust_var_regs)
from graph_making import (make_graphs, make_sequence_logo_graph, make_violin_plots, graphs_multiple_guide_lengths,
                          graphs_multiple_conditions, get_fungi_designs)



def check_checkpoint_file(coupled_folder):
    print('Checking that there are files to analyze or if checkpoint file is corrupted...')
    # Check if we have anything to test
    analysis_files = [file for file in os.listdir(coupled_folder) if file.endswith('.coupled')]

    if len(analysis_files) == 0:
        print('Please make sure to couple designs with the appropriate test sequences')
        return -1

    # check the amount of work by summing the number of designs on each file in the coupled folder
    lengths = [int(name.split('/')[-1].split('\\')[-1].split('_')[0]) for name in analysis_files]
    total_work = sum(lengths)

    # read each line as a tuple of three values - big index, file name, small index
    with open(os.path.normpath(coupled_folder + '/big_checkpoint.txt'), 'r') as handle:
        work_to_do_list = handle.read().splitlines()

    # Check the big work done checkpoint file and remove any indexes that are already there
    work_done_files = [os.path.normpath(f'{coupled_folder}/{f}') for f in os.listdir(coupled_folder)
                       if f.startswith('work_done_')]
    work_done = []
    for work_done_file in work_done_files:
        with open(work_done_file) as handle:
            work_done.extend(handle.read().splitlines())

    work_done = set(work_done)
    #  Basically remove make an array that contains all lines NOT in big work done
    work_to_do_list = [tuple(item.strip('\n').split('\t')) for item in work_to_do_list if
                       item.split('\t')[0] not in work_done]
    work_to_do = len(work_to_do_list)
    work_completed = len(work_done)
    if total_work != work_to_do + work_completed:
        if total_work == work_completed:
            decision = input('All work to do has been done, but checkpoint file is corrupted. Recommending deleting '
                            'checkpoint file and associated work / coupled files to re-do analysis. Would you like'
                            'to delete these? [Y/N]: ')
        else:
            decision = input(f'big_checkpoint.txt does not match amount of work to do. '
                             f'Would you like to delete checkpoint file and associated work / coupled files '
                             f'to start analysis over? [Y/N]: ')
        while decision != 'Y' and decision != 'N' and decision != 'y' and decision != 'n':
            decision = input(
                f'Please enter either Y or N: ')
        if decision == 'Y' or decision == 'y':
            print(f'Deleting files...')
            # delete files
            for file in os.listdir(coupled_folder):
                if file.startswith('work_done') and file.endswith('.txt'):
                    os.remove(os.path.normpath(coupled_folder + '/' + file))
                elif file == 'big_checkpoint.txt':
                    os.remove(os.path.normpath(coupled_folder + '/' + file))
                elif file.endswith('.coupled'):
                    os.remove(os.path.normpath(coupled_folder + '/' + file))
            print('Done! Please make sure to check any results files in this folder too to make sure they are correct.'
                  'Then rerun ribodesigner_routine to start that analysis again.')
        else:
            print('Ok! I won\'t delete any files. Now exiting...')
        return -1

    else:
        print('All looks good!')
        return 0

def run_local(output_folder, guide_len, num_of_workers:int = mp.cpu_count(), get_tm_nn:bool = True):
    # finally, we test! Below is for local
    # First check if checkpoint file is corrupted:
    result = check_checkpoint_file(coupled_folder=os.path.normpath(output_folder+ '/coupled'))
    if result == -1:
        return -1
    print(f'\n now testing {output_folder}... \n')
    coupled_file = os.path.normpath(output_folder + '/coupled')
    in_data = [(coupled_file, num_of_workers, j, 1, False, guide_len, get_tm_nn)
               for j in range(num_of_workers)]
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        with mp.Pool(processes=len(in_data)) as pool:
            out_data = pool.starmap(ribo_checker, in_data)
        bar()
    combine_data(os.path.normpath(output_folder + '/coupled/results'))
    return out_data


def run_remote(output_folder, guide_len, n_limit, scratch_path: str = None, number_of_workers: int = mp.cpu_count(),
               worker_number: int = 0):
    # This is for NOTS (make sure to upload coupled data with globus before you do this!)
    if scratch_path is None:
        scratch_path = input('What is the scratch path? ')
    print(f'\n now testing {output_folder}... \n')
    get_tm_nn = True
    # coupled_file = '/scratch/kpr1/RiboDesigner/' + output_folder + '/coupled/'
    coupled_file = os.path.normpath(scratch_path + output_folder + '/coupled')
    ribo_checker(coupled_folder=coupled_file, number_of_workers=number_of_workers, worker_number=worker_number,
                 n_limit=n_limit, opti_len=guide_len, get_tm_nn=get_tm_nn)
    combine_data(coupled_file + 'results')
    return

def ribodesigner_routine(target_seqs_to_process: list, test_seqs_to_process: list, out_path: str, ref_seq_file: str,
                         guide_len: int = 50, igs_len: int = 5, min_len: int = 35, graph_results: bool = True,
                         var_regs=None, graph_type: str = 'png', get_consensus_batches: bool = True,
                         batch_num: int = 10, score_type: str = 'weighted', msa_fast: bool = True,
                         remove_x_dupes_in_graph: bool = True, var_regs_lim: int = 1580, min_true_cov: float = 0,
                         percent_of_target_seqs_used: float = 1, gaps_allowed: bool = False,
                         random_guide_sample_size: int = 10, flexible_igs: bool = True):
    if var_regs is None:
        var_regs = [(69, 99), (137, 242), (433, 497), (576, 682), (822, 879), (986, 1043), (1117, 1173), (1243, 1294),
                    (1435, 1465)]
    all_test_file_names = []
    for test_file in test_seqs_to_process:
        title = test_file.split('.')[0].split('/')[-1].split('\\')[-1]
        test_save_file_name = os.path.normpath(f'{out_path}/test_sequences_{title}.pickle')
        all_test_file_names.append(test_save_file_name)
        if not os.path.exists(test_save_file_name):
            _ = prepare_test_seqs(test_folder=test_file, ref_sequence_file=ref_seq_file, guide_length=guide_len,
                                  igs_length=igs_len, min_length=min_len, folder_to_save=out_path,
                                  graph_results=graph_results, var_regs=var_regs, graph_file_type=graph_type,
                                  get_consensus_batches=get_consensus_batches, batch_num=batch_num,
                                  score_type=score_type, msa_fast=msa_fast,
                                  remove_x_dupes_in_graph=remove_x_dupes_in_graph, lim=var_regs_lim)
        else:
            print(f'{test_save_file_name} exists already! Moving on...')
    all_target_file_names = []
    for target_file in target_seqs_to_process:
        target_title = target_file.split('.')[0].split('/')[-1].split('\\')[-1]
        target_save_file_name = os.path.normpath(f'{out_path}/designs_{target_title}_universal.pickle')
        all_target_file_names.append(target_save_file_name)

        if not os.path.exists(target_save_file_name):
            _ = ribodesigner(target_sequences_folder=target_file, ref_sequence_file=ref_seq_file,
                             guide_length=guide_len, igs_length=igs_len, min_length=min_len, fileout=True,
                             folder_to_save=out_path, min_true_cov=min_true_cov, msa_fast=msa_fast,
                             score_type=score_type, percent_of_target_seqs_used=percent_of_target_seqs_used,
                             gaps_allowed=gaps_allowed, random_guide_sample_size=random_guide_sample_size)
        else:
            print(f'{target_save_file_name} exists already! Moving on...')

        for test_outfile in all_test_file_names:
            _ = couple_designs_to_test_seqs(designs_input=target_save_file_name, test_seqs_input=test_outfile,
                                            flexible_igs=flexible_igs, file_to_save=out_path)
        return all_target_file_names, all_test_file_names

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
        number_of_workers = mp.cpu_count()
        pass

    # Barcode sequence is split sfGFP just cuz. This does not affect guide sequence design.
    barcode_seq_file = os.path.normpath('Common_sequences/sfGFP_2_seq_barcode.txt')

    # We'll be using the Ribozyme published in the RAM paper
    ribobody_file = os.path.normpath('Common_sequences/ribozyme_body.txt')

    # Prepare the datasets
    datasets_path = os.path.normpath('Datasets_used/zymo_files/')

    # Output folder
    output_path = os.path.normpath('test_output_files/')

    # Reference sequence - will have to be E. coli to graph the variable regions
    ref_path = os.path.normpath('Common_sequences/e-coli-16s-mg1655.fasta')
    ref_path_arc = os.path.normpath('Common_sequences/Methanobrevibacter smithii 16s.fasta')
    ref_path_euk = os.path.normpath('Common_sequences/Saccharomyces cerevisiae 18s.fasta')

    # Our random sequence!
    random_seq_path = os.path.normpath('Common_sequences/Lactobacillus_casei_example.fasta')

    u1376 = 'CAACCCACTCCCATGGTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGTgTTCAC'

    # this is from Reich, M. & Labes, A. How to boost marine fungal research: A first step towards a multidisciplinary
    # approach by combining molecular fungal ecology and natural products chemistry. Marine Genomics 36, 57-75 (2017).
    s_cerevisiae_var_regs = [(69, 80), (126, 292), (478, 510), (643, 850), (1048, 1070), (1350, 1400), (1480, 1531),
                             (1674, 1730)]

    # Chakravorty, S., Helb, D., Burday, M., Connell, N. & Alland, D. A detailed analysis of 16S ribosomal RNA gene
    # segments for the diagnosis of pathogenic bacteria. J Microbiol Methods 69, 330-339 (2007).
    e_coli_var_regs = [(69, 99), (137, 242), (433, 497), (576, 682), (822, 879), (986, 1043), (1117, 1173),
                       (1243, 1294), (1435, 1465)]

    # get archaeal 16s variable regions by aligning to e_coli variable regions. This way we can graph anything aligned
    # to m_smithii with its own index
    m_smithii_var_regs = adjust_var_regs(known_seq_file=ref_path, known_var_regs=e_coli_var_regs,
                                         unknown_seq_file=ref_path_arc)

    # find equivalent location of our original design in archaea and eukaryota
    arch_1376_idx = adjust_var_regs(known_seq_file=ref_path, known_var_regs=1376, unknown_seq_file=ref_path_arc)
    euk_1376_idx = adjust_var_regs(known_seq_file=ref_path, known_var_regs=1376, unknown_seq_file=ref_path_euk)

    # ########################################################
    # # test data targeted
    # path = 'Datasets_used/SILVA_squished_datasets_1_per_genus/'
    # bad_targets = 'Datasets_used/Bacillus_halotolerans.fasta'
    # universal_data_1 = path + 'SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_2.fasta'
    # universal_data_2 = path + 'SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_2.fasta'
    # universal_data_3 = path + 'SILVA_squished_datasets_Eukaryota_Only/Eukaryota_Only_by_Genus_2.fasta'
    # universal_data_4 = path + 'SILVA_squished_datasets_All_Kingdoms/All_by_Genus_2.fasta'
    # big_data_entero_only = path + 'Enterobacterales_only_squished/Enterobacterales_only_by_Genus_1.fasta'
    # big_data_pseudo_only = path + 'Pseudomonadales_only_squished/Pseudomonadales_only_by_Genus_1.fasta'
    # big_data_no_entero_or_pseudo = path + 'Background_Bacteria_squished/' \
    #                                       'Background_Bacteria_squished_no_pseudo_or_entero.fasta'
    # big_data_no_pseudo = path + 'Background_Bacteria_squished/Background_Bacteria_squished_no_pseudo.fasta'
    # big_data_no_entero = path + 'Background_Bacteria_squished/Background_Bacteria_squished_no_entero.fasta'
    # big_data_only_entero_and_pseudo = path + 'Pseudo_and_entero_only_squished/Pseudo_and_entero_only_by_Genus_1.fasta'
    # big_data_gram_pos_only = path + 'Gram_positives_only/Gram_positives_only.fasta'
    # big_data_no_gram_pos = path + 'No_Gram_positives/No_Gram_positives.fasta'
    # background_data_bac = path + 'SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_1.fasta'
    # background_data_arc = path + 'SILVA_squished_datasets_Archaea_Only/Archaea_Only_by_Genus_1.fasta'
    # background_data_euk = path + 'SILVA_squished_datasets_Eukaryota_Only/Eukaryota_Only_by_Genus_1.fasta'
    # background_data_all = path + 'SILVA_squished_datasets_All_Kingdoms/All_by_Genus_1.fasta'
    # test_output_folder = 'test_output_files/test_outputs_parallelizing'
    # test_file = 'test_dataset_for_graphs.csv'
    # big_data_file_for_output = 'large_dataset.csv'
    # ref_analysis_folder = 'test_output_files/test_outputs_parallelizing/native ecoli mg1655 designs'
    #
    # test_data_folders = [background_data_euk, background_data_arc, background_data_bac, background_data_all]
    # selective_data_folders_targets = [big_data_no_entero_or_pseudo, big_data_entero_only, big_data_pseudo_only,
    #                                   big_data_gram_pos_only]
    # selective_data_folders_tests = [big_data_only_entero_and_pseudo, big_data_no_entero, big_data_no_pseudo,
    #                                 big_data_no_gram_pos]
    # universal_data_files = [universal_data_1, universal_data_2, universal_data_3, universal_data_4]
    # test_data_folders_test = [bad_targets, big_data_entero_only]
    # # Test new RiboDesigner for images
    # universal_datasets = []
    # selective_datasets = []
    #
    # test_data_pickles = ['test_sequences_All_by_Genus_1.pickle', 'test_sequences_Archaea_Only_by_Genus_1.pickle',
    #                      'test_sequences_Bacteria_Only_by_Genus_1.pickle',
    #                      'test_sequences_Eukaryota_Only_by_Genus_1.pickle',
    #                      'test_sequences_All_by_Genus_1_batched.pickle',
    #                      'test_sequences_Archaea_Only_by_Genus_1_batched.pickle',
    #                      'test_sequences_Bacteria_Only_by_Genus_1_batched.pickle',
    #                      'test_sequences_Eukaryota_Only_by_Genus_1_batched.pickle']
    # test_data_pickles = [test_output_folder + '/' + item for item in test_data_pickles]
    #
    # universal_data_pickles = ['designs_Eukaryota_Only_by_Genus_2_universal.pickle',
    #                           'designs_All_by_Genus_2_universal.pickle',
    #                           'designs_Archaea_Only_by_Genus_2_universal.pickle',
    #                           'designs_Bacteria_Only_by_Genus_2_universal.pickle']
    # universal_data_pickles = [test_output_folder + '/' + item for item in universal_data_pickles]
    #
    # ref_seq_pickle = test_output_folder + '/designs_e-coli-16s-mg1655_universal.pickle'
    # random_seq_pickle = test_output_folder + '/designs_Lactobacillus_casei_example_universal.pickle'

    # # Here we make the designs batched
    # for test_data in test_data_folders:
    #     test_data = background_data_bac
    #     test_seqs_pickle_file_name = prepare_test_seqs(test_folder=test_data, ref_sequence_file=ref_path,
    #                                                    guide_length=n, igs_length=m, min_length=minlen,
    #                                                    folder_to_save=test_output_folder, graph_results=True,
    #                                                    var_regs=e_coli_var_regs, graph_file_type='png',
    #                                                    get_consensus_batches=True, batch_num=10, score_type='weighted',
    #                                                    msa_fast=True, remove_x_dupes_in_graph=True)
    #     test_data_pickles.append(test_seqs_pickle_file_name)
    # # #
    # # unbatched
    # for test_data in test_data_folders:
    #     test_seqs_pickle_file_name = prepare_test_seqs(test_folder=test_data, ref_sequence_file=ref_path,
    #                                                    guide_length=n, igs_length=m, min_length=minlen,
    #                                                    folder_to_save=test_output_folder, graph_results=True,
    #                                                    var_regs=e_coli_var_regs, graph_file_type='png',
    #                                                    get_consensus_batches=False,  batch_num=10,
    #                                                    score_type='weighted', msa_fast=True)
    #     test_data_pickles.append(test_seqs_pickle_file_name)

    # # Here, we're using ribodesigner functions to see what would happen if we used the native sequences after each
    # # U site as guides in E. coli MG1655
    # ref_seq_pickle_file_name = ribodesigner(target_sequences_folder=ref_path,
    #                                         ref_sequence_file=ref_path, igs_length=m,
    #                                         guide_length=n, min_length=minlen, selective=False, min_true_cov=0,
    #                                         msa_fast=True,
    #                                         score_type='weighted', n_limit=1, percent_of_target_seqs_used=1,
    #                                         gaps_allowed=False, fileout=False, random_guide_sample_size=10,
    #                                         folder_to_save=test_output_folder)
    # # Similarly, what if we chose another sequence
    # random_seq_pickle_file_name = ribodesigner(target_sequences_folder=random_seq_path,
    #                                            ref_sequence_file=ref_path, igs_length=m,
    #                                            guide_length=n, min_length=minlen, selective=False, min_true_cov=0,
    #                                            msa_fast=True,
    #                                            score_type='weighted', n_limit=1, percent_of_target_seqs_used=1,
    #                                            gaps_allowed=False, fileout=False, random_guide_sample_size=10,
    #                                            folder_to_save=test_output_folder)
    # # finally, here are our universal designs
    # for universal_data in universal_data_files:
    #     universal_pickle_file_name = ribodesigner(target_sequences_folder=universal_data,
    #                                               ref_sequence_file=ref_path, igs_length=m,
    #                                               guide_length=n, min_length=minlen, selective=False, min_true_cov=0,
    #                                               msa_fast=True,
    #                                               score_type='weighted', n_limit=1, percent_of_target_seqs_used=1,
    #                                               gaps_allowed=False, fileout=False, random_guide_sample_size=10,
    #                                               folder_to_save=test_output_folder)
    #     universal_data_pickles.append(universal_pickle_file_name)

    # # Now we couple the designs with their test sequences to later test them
    # if os.path.exists(test_output_folder + '/coupled/big_checkpoint.txt'):
    #     os.remove(test_output_folder + '/coupled/big_checkpoint.txt')
    # for test_pickle in test_data_pickles:
    #     random_seq_coupled_file_name = couple_designs_to_test_seqs(designs_input=random_seq_pickle,
    #                                                                test_seqs_input=test_pickle,
    #                                                                flexible_igs=True, igs_len=m,
    #                                                                score_type='weighted',
    #                                                                file_to_save=test_output_folder)
    #
    #     ref_seq_coupled_file_name = couple_designs_to_test_seqs(designs_input=ref_seq_pickle,
    #                                                             test_seqs_input=test_pickle,
    #                                                             flexible_igs=True, igs_len=m,
    #                                                             score_type='weighted',
    #                                                             file_to_save=test_output_folder)
    #
    #     control_design_coupled_file_name = couple_designs_to_test_seqs(designs_input=u1376,
    #                                                                    test_seqs_input=test_pickle,
    #                                                                    flexible_igs=True, igs_len=m, ref_idx_u_loc=1376,
    #                                                                    score_type='weighted',
    #                                                                    file_to_save=test_output_folder)
    #
    #     for universal_design_pickle in universal_data_pickles:
    #         universal_designs_coupled_file_name = couple_designs_to_test_seqs(designs_input=universal_design_pickle,
    #                                                                           test_seqs_input=test_pickle,
    #                                                                           flexible_igs=True,
    #                                                                           file_to_save=test_output_folder)

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
    # get_tm_nn = True
    # files_to_test = test_output_folder + '/coupled'
    # in_data = [(files_to_test, number_of_workers, i, 1, False, n, get_tm_nn) for i in range(number_of_workers)]
    # with alive_bar(unknown='fish', spinner='fishes') as bar:
    #     with mp.Pool(processes=len(in_data)) as pool:
    #         out_data = pool.starmap(ribo_checker, in_data)
    #     bar()

    # this is not parallelized but still for local
    # for i in range(number_of_workers):
    #     print(f'Worker number {i}')
    #     ribo_checker(coupled_folder=files_to_test, number_of_workers=number_of_workers, worker_number=i,
    #                  n_limit=1, opti_len=n, get_tm_nn=get_tm_nn)
    # # This is for NOTS
    # files_to_test = test_output_folder + '/coupled'
    # ribo_checker(coupled_folder='/scratch/kpr1/RiboDesigner/' + files_to_test, number_of_workers=number_of_workers,
    #              worker_number=worker_number, n_limit=1, get_tm_nn=get_tm_nn)

    # print(f'Test data done!\n########################################################\n')

    # results_folder = test_output_folder + '/coupled/results'
    # # combine_data(results_folder)
    #
    # control_design_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                      if f.startswith('1_designs')]
    # control_design_results_file_names.sort()
    # batched_control_design_results_file_names = [name for name in control_design_results_file_names if
    #                                              'batched' in name]
    # unbatched_control_design_results_file_names = [name for name in control_design_results_file_names if
    #                                                not 'batched' in name]
    #
    # ref_design_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                  if f.startswith('304_designs')]
    # ref_design_results_file_names.sort()
    # batched_ref_design_results_file_names = [name for name in ref_design_results_file_names if 'batched' in name]
    # unbatched_ref_design_results_file_names = [name for name in ref_design_results_file_names if not 'batched' in name]
    #
    # random_design_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                     if f.startswith('326_designs')]
    # random_design_results_file_names.sort()
    # batched_random_design_results_file_names = [name for name in random_design_results_file_names if 'batched' in name]
    # unbatched_random_design_results_file_names = [name for name in random_design_results_file_names if
    #                                               not 'batched' in name]
    #
    # universal_design_archaea_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                                if f.startswith('14652_designs')]
    # universal_design_archaea_results_file_names.sort()
    # batched_universal_design_archaea_results_file_names = [name for name in universal_design_archaea_results_file_names
    #                                                        if 'batched' in name]
    # unbatched_universal_design_archaea_results_file_names = [name for name in
    #                                                          universal_design_archaea_results_file_names if
    #                                                          not 'batched' in name]
    #
    # universal_design_eukarya_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                                if f.startswith('72490_designs')]
    # universal_design_eukarya_results_file_names.sort()
    # batched_universal_design_eukarya_results_file_names = [name for name in universal_design_eukarya_results_file_names
    #                                                        if 'batched' in name]
    # unbatched_universal_design_eukarya_results_file_names = [name for name in
    #                                                          universal_design_eukarya_results_file_names if
    #                                                          not 'batched' in name]
    #
    # universal_design_bacteria_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                                 if f.startswith('77193_designs')]
    # universal_design_bacteria_results_file_names.sort()
    # batched_universal_design_bacteria_results_file_names = [name for name in
    #                                                         universal_design_bacteria_results_file_names if
    #                                                         'batched' in name]
    # unbatched_universal_design_bacteria_results_file_names = [name for name in
    #                                                           universal_design_bacteria_results_file_names if
    #                                                           not 'batched' in name]
    #
    # universal_design_all_results_file_names = [f'{results_folder}/{f}' for f in os.listdir(results_folder)
    #                                            if f.startswith('137620_designs')]
    # universal_design_all_results_file_names.sort()
    # batched_universal_design_all_results_file_names = [name for name in universal_design_all_results_file_names if
    #                                                    'batched' in name]
    # unbatched_universal_design_all_results_file_names = [name for name in universal_design_all_results_file_names if
    #                                                      not 'batched' in name]
    #
    # batched_all_targets_universal_design_file_names = [*batched_universal_design_archaea_results_file_names,
    #                                                    *batched_universal_design_eukarya_results_file_names,
    #                                                    *batched_universal_design_bacteria_results_file_names,
    #                                                    *batched_universal_design_all_results_file_names]
    # unbatched_all_targets_universal_design_file_names = [*unbatched_universal_design_archaea_results_file_names,
    #                                                      *unbatched_universal_design_eukarya_results_file_names,
    #                                                      *unbatched_universal_design_bacteria_results_file_names,
    #                                                      *unbatched_universal_design_all_results_file_names]

    # # All data
    # make_graphs(control_designs_path=control_design_results_file_names,
    #             universal_designs_path=universal_design_bacteria_results_file_names,
    #             ref_seq_designs_path=ref_design_results_file_names,
    #             random_seq_designs_path=random_design_results_file_names, var_regs=e_coli_var_regs,
    #             save_fig=True, save_file_loc=results_folder)

    # combined_results_folder = test_output_folder + '/coupled/results/combined'
    #
    # Batched
    # make_graphs(control_designs_path=batched_control_design_results_file_names,
    #             universal_designs_path=batched_all_targets_universal_design_file_names,
    #             ref_seq_designs_path=batched_ref_design_results_file_names,
    #             random_seq_designs_path=batched_random_design_results_file_names, var_regs=e_coli_var_regs,
    #             save_fig=True, save_file_loc=results_folder + '/batched')
    # folder_for_ref_seq_results = results_folder + '/batched' + '/ref_seq is MG1655'
    # folder_for_random_seq_results = results_folder + '/batched' + '/ref_seq is L casei'
    # #
    # batched_bacteria = combined_results_folder + '/' + batched_universal_design_bacteria_results_file_names[0].split('/')[-1][:-21] + '.txt'
    # batched_ref_seq = combined_results_folder + '/' + batched_ref_design_results_file_names[0].split('/')[-1][:-21] + '.txt'
    # batched_random_seq = combined_results_folder + '/' + batched_random_design_results_file_names[0].split('/')[-1][:-21] + '.txt'
    # #
    # make_sequence_logo_graph(test_data_path=test_data_pickles[2], design_data_path=[batched_bacteria],
    #                          ref_data_path=[batched_ref_seq], save_fig=True, save_file_loc=folder_for_ref_seq_results)
    #
    # make_sequence_logo_graph(test_data_path=test_data_pickles[2], design_data_path=[batched_bacteria],
    #                          ref_data_path=[batched_random_seq], save_fig=True,
    #                          save_file_loc=folder_for_random_seq_results)
    # # Unbatched
    # make_graphs(control_designs_path=unbatched_control_design_results_file_names,
    #             universal_designs_path=unbatched_all_targets_universal_design_file_names,
    #             ref_seq_designs_path=unbatched_ref_design_results_file_names,
    #             random_seq_designs_path=unbatched_random_design_results_file_names, var_regs=e_coli_var_regs,
    #             save_fig=True, save_file_loc=results_folder + '/unbatched')
    #
    # folder_for_ref_seq_results = results_folder + '/unbatched' + '/ref_seq is MG1655'
    # folder_for_random_seq_results = results_folder + '/unbatched' + '/ref_seq is L casei'
    #
    # unbatched_bacteria = combined_results_folder + '/' + unbatched_universal_design_bacteria_results_file_names[0].split('/')[-1][:-21] + '.txt'
    # unbatched_ref_seq = combined_results_folder + '/' + unbatched_ref_design_results_file_names[0].split('/')[-1][:-21] + '.txt'
    # unbatched_random_seq = combined_results_folder + '/' + unbatched_random_design_results_file_names[0].split('/')[-1][:-21] + '.txt'
    #
    # make_sequence_logo_graph(test_data_path=test_data_pickles[2], design_data_path=[unbatched_bacteria],
    #                          ref_data_path=[unbatched_ref_seq], save_fig=True, save_file_loc=folder_for_ref_seq_results)
    #
    # make_sequence_logo_graph(test_data_path=test_data_pickles[2], design_data_path=[unbatched_bacteria],
    #                          ref_data_path=[unbatched_random_seq], save_fig=True,
    #                          save_file_loc=folder_for_random_seq_results)

    # control_design_results_file_name = test_output_folder + '/coupled/results/1_designs_TTCAC1376_vs_test_sequences_Eukaryota_Only_by_Genus_1_worker_0_results.txt'
    # universal_designs_results_file_name = test_output_folder + '/coupled/results/72490_designs_designs_Eukaryota_Only_by_Genus_2_universal_vs_test_sequences_All_by_Genus_1_worker_0_results.txt'
    # ref_seq_results_file_name = test_output_folder + '/coupled/results/304_designs_designs_e-coli-16s-mg1655_universal_vs_test_sequences_Eukaryota_Only_by_Genus_1_worker_0_results.txt'
    # random_seq_results_file_name_1 = test_output_folder + '/coupled/results/326_designs_designs_Lactobacillus_casei_example_universal_vs_test_sequences_Bacteria_Only_by_Genus_1_worker_0_results.txt'
    # random_seq_results_file_name_2 = test_output_folder + '/coupled/results/326_designs_designs_Lactobacillus_casei_example_universal_vs_test_sequences_Eukaryota_Only_by_Genus_1_worker_0_results.txt'
    # random_seq_results_file_name_3 = test_output_folder + '/coupled/results/326_designs_designs_Lactobacillus_casei_example_universal_vs_test_sequences_All_Only_by_Genus_1_worker_0_results.txt'

    # combined_results_folder = test_output_folder + '/coupled/results/combined'
    # control_design_results_file_names = [f'{combined_results_folder}/{f}' for f in os.listdir(combined_results_folder)
    #                                      if f.startswith('1_designs')]
    # ref_design_results_file_names = [f'{combined_results_folder}/{f}' for f in os.listdir(combined_results_folder)
    #                                  if f.startswith('304_designs')]
    # random_design_results_file_names = [f'{combined_results_folder}/{f}' for f in os.listdir(combined_results_folder)
    #                                     if f.startswith('326_designs')]
    # universal_design_results_file_names = [f'{combined_results_folder}/{f}' for f in os.listdir(combined_results_folder)
    #                                        if f.startswith('77193_designs')]

    # folder_for_ref_seq_results = 'test_output_files/test_outputs_parallelizing/Coupled for testing/for_testing/figures/ref_seq is MG1655'
    # folder_for_random_seq_results = 'test_output_files/test_outputs_parallelizing/Coupled for testing/for_testing/figures/ref_seq is L casei'
    # make_graphs(control_designs_path=control_design_results_file_name,
    #             universal_designs_path=universal_designs_results_file_name,
    #             ref_seq_designs_path=ref_seq_results_file_name, var_regs=e_coli_var_regs, save_fig=False,
    #             save_file_loc='test_output_files/test_outputs_parallelizing/coupled/figures')
    #
    # ########################################################
    # # All data
    # make_graphs(control_designs_path=control_design_results_file_names,
    #             universal_designs_path=universal_design_bacteria_results_file_names,
    #             ref_seq_designs_path=ref_design_results_file_names, var_regs=e_coli_var_regs,
    #             save_fig=True, save_file_loc=folder_for_ref_seq_results)
    #
    # make_graphs(control_designs_path=control_design_results_file_names,
    #             universal_designs_path=universal_design_bacteria_results_file_names,
    #             ref_seq_designs_path=random_design_results_file_names, var_regs=e_coli_var_regs, save_fig=True,
    #             save_file_loc=folder_for_random_seq_results)
    # # Bacteria only
    # make_graphs(control_designs_path=[control_design_results_file_names[2]],
    #             universal_designs_path=[universal_design_results_file_names[5]],
    #             ref_seq_designs_path=[ref_design_results_file_names[2]], var_regs=e_coli_var_regs,
    #             save_fig=True, save_file_loc=folder_for_ref_seq_results)
    #
    # make_graphs(control_designs_path=[control_design_results_file_names[2]],
    #             universal_designs_path=[universal_design_results_file_names[5]],
    #             ref_seq_designs_path=[random_design_results_file_names[4]], var_regs=e_coli_var_regs, save_fig=True,
    #             save_file_loc=folder_for_random_seq_results)
    #
    # make_sequence_logo_graph(test_data_path=test_data_pickles[2],
    #                          design_data_path=[universal_design_results_file_names[5]],
    #                          ref_data_path=[ref_design_results_file_names[0]], save_fig=True,
    #                          save_file_loc=folder_for_ref_seq_results)
    #
    # make_sequence_logo_graph(test_data_path=test_data_pickles[2],
    #                          design_data_path=[universal_design_results_file_names[5]],
    #                          ref_data_path=[random_design_results_file_names[4]], save_fig=True,
    #                          save_file_loc=folder_for_random_seq_results)
    #
    # print(f'Graphs done!\n########################################################\n')

    ################################################################################################################

    # # The following code is to analyze the effects of changing the guide length on bacterial designs
    # out_folder_guide_lengths = 'test_output_files/varying_guide_lengths'
    # test_data_guide_len_pickles = []
    # start = 10
    #
    # for i in range(start, 210, 10):
    #     test_out_file = out_folder_guide_lengths + f'/{i}_bp'
    #     if not os.path.exists(test_out_file):
    #         os.mkdir(test_out_file)
    #     # test_seqs_pickle_file_name = prepare_test_seqs(test_folder=background_data_bac, ref_sequence_file=ref_path,
    #     #                                                guide_length=i, igs_length=m, min_length=i,
    #     #                                                folder_to_save=test_out_file, graph_results=True,
    #     #                                                var_regs=e_coli_var_regs, graph_file_type='png',
    #     #                                                get_consensus_batches=True, batch_num=10, score_type='weighted',
    #     #                                                msa_fast=True)
    #     # test_data_guide_len_pickles.append(test_seqs_pickle_file_name)
    #     test_seqs_pickle_file_name = test_out_file + '/test_sequences_Bacteria_Only_by_Genus_1.pickle'
    #
    #     # Now, we design bacterial sequences to go with these designs
    #     universal_pickle_file_name = ribodesigner(target_sequences_folder=universal_data_1, ref_sequence_file=ref_path,
    #                                               igs_length=m, guide_length=i, min_length=i, selective=False,
    #                                               min_true_cov=0, msa_fast=True, score_type='weighted', n_limit=1,
    #                                               percent_of_target_seqs_used=1, gaps_allowed=False, fileout=False,
    #                                               random_guide_sample_size=10, folder_to_save=test_out_file)
    #     universal_data_pickles.append(universal_pickle_file_name)
    #     # and we pair these designs
    #     universal_designs_coupled_file_name = couple_designs_to_test_seqs(designs_input=universal_pickle_file_name,
    #                                                                       test_seqs_input=test_seqs_pickle_file_name,
    #                                                                       flexible_igs=True, file_to_save=test_out_file)
    #     print(f'{i} bp datasets generated.\n########################################################\n')

    # # This is for NOTS (make sure to upload coupled data with globus before you do this!)
    # for i in range(10, 210, 10):
    #     test_out_file = f'test_output_files/varying_guide_lengths/{i}_bp/coupled'
    #     if not os.path.exists(test_out_file):
    #         # If we have not uploaded the file yet, skip this length
    #         print(f'No data found for guide length of {i} bp!')
    #         continue
    #     ribo_checker(coupled_folder='/scratch/kpr1/RiboDesigner/' + test_out_file, number_of_workers=number_of_workers,
    #                  worker_number=worker_number, n_limit=1, opti_len=i, get_tm_nn=True)
    #     print(f'{i} bp datasets analyzed.\n########################################################\n')

    # # finally, we test! Below is for local
    # get_tm_nn = True
    # files_to_test = test_output_folder + '/coupled'
    # for i in range(10, 210, 10):
    #     test_out_file = f'test_output_files/varying_guide_lengths/{i}_bp/coupled'
    #     if not os.path.exists(test_out_file):
    #         # If we have not uploaded the file yet, skip this length
    #         print(f'No data found for guide length of {i} bp!')
    #         continue
    #     in_data = [(test_out_file, number_of_workers, j, 1, False, i, get_tm_nn) for j in range(number_of_workers)]
    #     with alive_bar(unknown='fish', spinner='fishes') as bar:
    #         with mp.Pool(processes=len(in_data)) as pool:
    #             out_data = pool.starmap(ribo_checker, in_data)
    #         bar()
    #     print(f'{i} bp datasets analyzed.\n########################################################\n')

    # # post_nots_output_folder = 'test_output_files/varying_guide_lengths_NOTS_output'
    # post_local_output_folder = 'test_output_files/varying_guide_lengths'
    #
    # make_violin_plots(post_local_output_folder, vars_to_plot=['test_score', 'tm_nn_vs_test'], file_type='png')

    ################################################################################################################

    # # The following code is to analyze the effects of changing the guide length on our reference sequence
    # out_folder_guide_lengths_ref_seq = 'test_output_files/varying_guide_lengths_e_coli'
    # out_folder_guide_lengths_test_seqs = 'test_output_files/varying_guide_lengths'
    # test_data_guide_len_pickles = []
    # start = 10
    #
    # for i in range(start, 210, 10):
    #     ref_target_out_file = out_folder_guide_lengths_ref_seq + f'/{i}_bp'
    #     test_seqs_out_file = out_folder_guide_lengths_test_seqs + f'/{i}_bp'
    #     if not os.path.exists(ref_target_out_file):
    #         os.mkdir(ref_target_out_file)
    #
    #     # Now, we design bacterial sequences to go with these designs
    #     ref_pickle_file_name = ribodesigner(target_sequences_folder=ref_path, ref_sequence_file=ref_path, igs_length=m,
    #                                         guide_length=i, min_length=i, selective=False, min_true_cov=0,
    #                                         msa_fast=True, score_type='weighted', n_limit=1,
    #                                         percent_of_target_seqs_used=1, gaps_allowed=False, fileout=False,
    #                                         random_guide_sample_size=10, folder_to_save=ref_target_out_file)
    #     universal_data_pickles.append(ref_pickle_file_name)
    #     test_seqs_out_file = test_seqs_out_file + '/test_sequences_Bacteria_Only_by_Genus_1.pickle'
    #     # and we pair these designs
    #     universal_designs_coupled_file_name = couple_designs_to_test_seqs(designs_input=ref_pickle_file_name,
    #                                                                       test_seqs_input=test_seqs_out_file,
    #                                                                       flexible_igs=True,
    #                                                                       file_to_save=ref_target_out_file)
    #     print(f'{i} bp datasets generated.\n########################################################\n')
    #
    # # finally, we test! Below is for local
    # get_tm_nn = True
    # for i in range(10, 210, 10):
    #     test_out_file = f'test_output_files/varying_guide_lengths_e_coli/{i}_bp/coupled'
    #     if not os.path.exists(test_out_file):
    #         # If we have not uploaded the file yet, skip this length
    #         print(f'No data found for guide length of {i} bp!')
    #         continue
    #     in_data = [(test_out_file, number_of_workers, j, 1, False, i, get_tm_nn) for j in range(number_of_workers)]
    #     with alive_bar(unknown='fish', spinner='fishes') as bar:
    #         with mp.Pool(processes=len(in_data)) as pool:
    #             out_data = pool.starmap(ribo_checker, in_data)
    #         bar()
    #     print(f'{i} bp datasets analyzed.\n########################################################\n')
    #
    # post_local_output_folder_e_coli = 'test_output_files/varying_guide_lengths_e_coli'
    # post_local_output_folder = 'test_output_files/varying_guide_lengths'
    #
    # make_violin_plots([('E. coli', post_local_output_folder_e_coli),
    #                    ('Bacterial designs', post_local_output_folder)],
    #                   vars_to_plot=['test_score', 'tm_nn_vs_test'], folder_to_save=post_local_output_folder_e_coli,
    #                   file_type='png')

    # ################################################################################################################
    # # Selective designs! Testing guide lengths from 10 to 40 to complete my 10 to 50 dataset
    # # selective_background_data_pickles = []
    # '/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/My Drive/KRG Thesis/test_outputs_selective'
    # selective_output_folder_root = 'test_output_files/test_outputs_selective'
    # selective_background_data_pickles = []
    # for i in range(10, 50, 10):
    #     selective_output_folder = selective_output_folder_root + f'/{i}_bp'
    #     for test_data in selective_data_folders_tests:
    #         print(f'Now testing {test_data}')
    #         test_seqs_pickle_file_name = prepare_test_seqs(test_folder=test_data, ref_sequence_file=ref_path,
    #                                                        guide_length=i, igs_length=m, min_length=minlen,
    #                                                        folder_to_save=selective_output_folder, graph_results=True,
    #                                                        var_regs=e_coli_var_regs, graph_file_type='png',
    #                                                        get_consensus_batches=True, batch_num=10,
    #                                                        score_type='weighted',  msa_fast=True,
    #                                                        remove_x_dupes_in_graph=True)
    #         selective_background_data_pickles.append(test_seqs_pickle_file_name)
    #     selective_target_data_pickles = []
    #     for target_path in selective_data_folders_targets:
    #         print(f'Now generating designs for {target_path}')
    #         target_seq_file_name = ribodesigner(target_sequences_folder=target_path, ref_sequence_file=ref_path,
    #                                             igs_length=m, guide_length=i, min_length=minlen, selective=False,
    #                                             min_true_cov=0, msa_fast=True, score_type='weighted', n_limit=1,
    #                                             percent_of_target_seqs_used=1, gaps_allowed=False, fileout=False,
    #                                             random_guide_sample_size=10, folder_to_save=selective_output_folder)
    #         selective_target_data_pickles.append(target_seq_file_name)
    #     # selective_target_data_pickles = ['designs_Background_Bacteria_squished_no_pseudo_or_entero_universal.pickle',
    #     #                                  'designs_Enterobacterales_only_by_Genus_1_universal.pickle',
    #     #                                  'designs_Pseudomonadales_only_by_Genus_1_universal.pickle',
    #     #                                  'designs_Gram_positives_only_universal.pickle']
    #     # selective_target_data_pickles = [f'{selective_output_folder}/{file}' for file in selective_target_data_pickles]
    #     # selective_background_data_pickles = ['test_sequences_Pseudo_and_entero_only_by_Genus_1.pickle',
    #     #                                      'test_sequences_Background_Bacteria_squished_no_entero.pickle',
    #     #                                      'test_sequences_Background_Bacteria_squished_no_pseudo.pickle',
    #     #                                      'test_sequences_No_Gram_positives.pickle']
    #     # selective_background_data_pickles = [f'{selective_output_folder}/{file}' for file in
    #     #                                      selective_background_data_pickles]
    #
    #     for target, background in zip(selective_target_data_pickles, selective_background_data_pickles):
    #         print(f'Now coupling designs for {target} vs. {background} test sequences')
    #         _ = couple_designs_to_test_seqs(designs_input=target, test_seqs_input=background, flexible_igs=True,
    #                                         file_to_save=selective_output_folder)
    #         # Also go ahead and evaluate against all bacteria dataset
    #         _ = couple_designs_to_test_seqs(designs_input=target, test_seqs_input=test_data_pickles[2],
    #                                         flexible_igs=True, file_to_save=selective_output_folder)
    #
    #     # finally, we test! Below is for local
    #     print(f'\n now testing {i} bp... \n')
    #     get_tm_nn = True
    #     coupled_file = selective_output_folder + '/coupled/'
    #     in_data = [(coupled_file, number_of_workers, j, 1, False, i, get_tm_nn) for j in range(number_of_workers)]
    #     with alive_bar(unknown='fish', spinner='fishes') as bar:
    #         with mp.Pool(processes=len(in_data)) as pool:
    #             out_data = pool.starmap(ribo_checker, in_data)
    #         bar()
    # test_data_selective = ['Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished/Background_Bacteria_squished_no_pseudo_or_entero_2.fasta',
    #  'Datasets_used/SILVA_squished_datasets_1_per_genus/Enterobacterales_only_squished/Enterobacterales_only_by_Genus_2.fasta',
    #  'Datasets_used/SILVA_squished_datasets_1_per_genus/Pseudomonadales_only_squished/Pseudomonadales_only_by_Genus_2.fasta',
    #  'Datasets_used/SILVA_squished_datasets_1_per_genus/Gram_positives_only/Gram_positives_only_2.fasta']
    # for i in range(20, 60, 10):
    #     selective_output_folder = f'test_output_files/test_outputs_selective/{i}_bp'
    #     selective_background_data_pickles = []
    #     for test_data in test_data_selective:
    #         test_seqs_pickle_file_name = prepare_test_seqs(test_folder=test_data, ref_sequence_file=ref_path,
    #                                                        guide_length=i, igs_length=m, min_length=minlen,
    #                                                        folder_to_save=selective_output_folder, graph_results=True,
    #                                                        var_regs=e_coli_var_regs, graph_file_type='png',
    #                                                        get_consensus_batches=True, batch_num=10,
    #                                                        score_type='weighted',  msa_fast=True,
    #                                                        remove_x_dupes_in_graph=True)
    #         selective_background_data_pickles.append(test_seqs_pickle_file_name)
    #     selective_target_data_pickles = ['designs_Background_Bacteria_squished_no_pseudo_or_entero_universal.pickle',
    #                                          'designs_Enterobacterales_only_by_Genus_1_universal.pickle',
    #                                          'designs_Pseudomonadales_only_by_Genus_1_universal.pickle',
    #                                          'designs_Gram_positives_only_universal.pickle']
    #     selective_target_data_pickles = [f'{selective_output_folder}/{file}' for file in selective_target_data_pickles]
    #     selective_background_data_pickles = ['test_sequences_Background_Bacteria_squished_no_pseudo_or_entero_2.pickle',
    #                                              'test_sequences_Enterobacterales_only_by_Genus_2.pickle',
    #                                              'test_sequences_Pseudomonadales_only_by_Genus_2.pickle',
    #                                              'test_sequences_Gram_positives_only_2.pickle']
    #     selective_background_data_pickles = [f'{selective_output_folder}/{file}' for file in
    #                                              selective_background_data_pickles]
    #     out_file = f'test_output_files/test_outputs_selective/out_file/{i}_bp'
    #     for target, background in zip(selective_target_data_pickles, selective_background_data_pickles):
    #         print(f'Now coupling designs for {target} vs. {background} test sequences')
    #         _ = couple_designs_to_test_seqs(designs_input=target, test_seqs_input=background, flexible_igs=True,
    #                                         file_to_save=out_file)
    #     # combine_data(f'test_output_files/test_outputs_selective/{i}_bp/coupled/results')
    #     # finally, we test! Below is for local
    #     print(f'\n now testing {i} bp... \n')
    #     get_tm_nn = True
    #     coupled_file = out_file + '/coupled/'
    #     in_data = [(coupled_file, number_of_workers, j, 1, False, i, get_tm_nn) for j in range(number_of_workers)]
    #     with alive_bar(unknown='fish', spinner='fishes') as bar:
    #         with mp.Pool(processes=len(in_data)) as pool:
    #             out_data = pool.starmap(ribo_checker, in_data)
    #         bar()
    # # selective_designs_results_file_names = [f'{coupled_file}results/{f}' for f in os.listdir(coupled_file + 'results')
    # #                                         if f.endswith('results.txt')]
    # for i in range(10, 60, 10):
    #     out_file = f'test_output_files/test_outputs_selective/out_file/{i}_bp'
    #     combine_data(out_file + '/coupled/results/')
    # make_graphs(control_designs_path=batched_control_design_results_file_names,
    #             universal_designs_path=batched_all_targets_universal_design_file_names,
    #             ref_seq_designs_path=batched_ref_design_results_file_names,
    #             random_seq_designs_path=batched_random_design_results_file_names,
    #             selective_designs_path=selective_designs_results_file_names,
    #             var_regs=e_coli_var_regs,
    #             save_fig=True, save_file_loc=results_folder + '/batched')

    # for i in range(50, 60, 10):
    #     output_folder = f'test_output_files/varying_guide_lengths_universal_all/{i}_bp'
    #     # test_seqs_pickle_file_name = prepare_test_seqs(test_folder=background_data_all, ref_sequence_file=ref_path,
    #     #                                                guide_length=i, igs_length=m, min_length=minlen,
    #     #                                                folder_to_save=output_folder, graph_results=True,
    #     #                                                var_regs=e_coli_var_regs, graph_file_type='png',
    #     #                                                get_consensus_batches=True, batch_num=10, score_type='weighted',
    #     #                                                msa_fast=True, remove_x_dupes_in_graph=True)
    #     # target_seq_file_name = ribodesigner(target_sequences_folder=universal_data_4, ref_sequence_file=ref_path,
    #     #                                     igs_length=m, guide_length=i, min_length=minlen, selective=False,
    #     #                                     min_true_cov=0, msa_fast=True, score_type='weighted', n_limit=1,
    #     #                                     percent_of_target_seqs_used=1, gaps_allowed=False, fileout=False,
    #     #                                     random_guide_sample_size=10, folder_to_save=output_folder)
    #     #
    #     # _ = couple_designs_to_test_seqs(designs_input=target_seq_file_name, test_seqs_input=test_seqs_pickle_file_name, flexible_igs=True,
    #     #                                 file_to_save=output_folder)
    #     # finally, we test! Below is for local
    #     print(f'\n now testing {i} bp... \n')
    #     get_tm_nn = True
    #     coupled_file = output_folder + '/coupled/'
    #     in_data = [(coupled_file, number_of_workers, j, 1, False, i, get_tm_nn) for j in range(number_of_workers)]
    #     with alive_bar(unknown='fish', spinner='fishes') as bar:
    #         with mp.Pool(processes=len(in_data)) as pool:
    #             out_data = pool.starmap(ribo_checker, in_data)
    #         bar()
    #     combine_data(output_folder + '/coupled/results/')

    # graphs_multiple_guide_lengths(universal_path='test_output_files/varying_guide_lengths',
    #                               selective_path='test_output_files/test_outputs_selective',
    #                               output_folder='test_output_files/best_designs', add_overhangs=True)

    # print('Selecting best designs...')
    # order_these = []
    # out_file = 'test_output_files/best_designs'
    # # Get two best universal designs from bacteria at v8-9
    # test_files = batched_universal_design_bacteria_results_file_names
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=2,
    #                                                   results_folder=out_file, design_type='universal', igs_min=0.7,
    #                                                   guide_min=0.7, choose_var_reg_site=True,
    #                                                   file_extra_text='_best_v8_to_9_bacteria', add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # # Get best bacteria design at v3-v4
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=1,
    #                                                   results_folder=out_file,
    #                                   design_type='universal', igs_min=0.7, guide_min=0.7, choose_var_reg_site=True,
    #                                   start_idx=497, end_idx=576, file_extra_text='_best_v3_to_4_bacteria',
    #                                   add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # # Get three best universal designs at v8-9
    # test_files = batched_universal_design_all_results_file_names
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=3, results_folder=out_file,
    #                                   design_type='universal', igs_min=0.7, guide_min=0.7, choose_var_reg_site=True,
    #                                   file_extra_text='_best_true_universal', add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # # Get one best selective design for each testing set
    # test_files = [f'{coupled_file}results/{f}' for f in os.listdir(coupled_file + 'results') if
    #               'Enterobacterales_only_by_Genus_1_universal_vs_test_sequences_Background_Bacteria_squished_no_entero'
    #               in f]
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=1, results_folder=out_file,
    #                                   design_type='selective', igs_min=0.7, guide_min=0, choose_var_reg_site=False,
    #                                   file_extra_text='_best_for_entero_good_igs', add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # test_files = [f'{coupled_file}results/{f}' for f in os.listdir(coupled_file + 'results') if
    #               'no_pseudo_or_entero_universal_vs_test_sequences_Pseudo_and_entero_only' in f]
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=1, results_folder=out_file,
    #                                   design_type='selective', igs_min=0.7, guide_min=0, choose_var_reg_site=False,
    #                                   file_extra_text='_best_for_all_but_entero_or_pseudo_good_igs', add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # test_files = [f'{coupled_file}results/{f}' for f in os.listdir(coupled_file + 'results') if
    #               'Pseudomonadales_only_by_Genus_1_universal_vs_test_sequences_Background_Bacteria_squished_no_pseudo'
    #               in f]
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=1, results_folder=out_file,
    #                                   design_type='selective', igs_min=0.7, guide_min=0, choose_var_reg_site=False,
    #                                   file_extra_text='_best_for_pseudo_good_igs', add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # test_files = [f'{coupled_file}results/{f}' for f in os.listdir(coupled_file + 'results') if 'Gram_positives' in f]
    # design_with_overhangs, design_df = select_designs(tested_to_targets_path=test_files, designs_required=1, results_folder=out_file,
    #                                   design_type='selective', igs_min=0.7, guide_min=0, choose_var_reg_site=False,
    #                                   file_extra_text='_best_for_Gram_positives', add_overhangs=True)
    # order_these.append(design_with_overhangs)
    #
    # print('Best designs selected')
    # print(order_these)

    # # Ok here we're going to generate archaeal designs w archaeal ref seq and eukaryotic designs w eukaryotic ref seq
    # files = [(background_data_arc, ref_path_arc, m_smithii_var_regs, 'Archaea', universal_data_2, 1450),
    #          (background_data_euk, ref_path_euk, s_cerevisiae_var_regs, 'Eukaryota', universal_data_3, 1850),
    #          (background_data_bac, ref_path, e_coli_var_regs, 'Bacteria', universal_data_1, 1580),
    #          (background_data_all, ref_path, e_coli_var_regs, 'All', universal_data_4, 1580)]
    # for file, ref_seq, var_regs, a, target, lim in files:
    #     out = f'{output_path}universal_diff_var_regs/{a}'
    #     if a == 'All':
    #         test_save_file_name = f'{out}/test_sequences_{a}_by_Genus_1.pickle'
    #     else:
    #         test_save_file_name = f'{out}/test_sequences_{a}_Only_by_Genus_1.pickle'
    #     if not os.path.exists(test_save_file_name):
    #         test_seqs_pickle_file_name = prepare_test_seqs(test_folder=file, ref_sequence_file=ref_seq, guide_length=n,
    #                                                        igs_length=m, min_length=minlen, folder_to_save=out,
    #                                                        graph_results=True, var_regs=var_regs, graph_file_type='png',
    #                                                        get_consensus_batches=True, batch_num=10,
    #                                                        score_type='weighted', msa_fast=True,
    #                                                        remove_x_dupes_in_graph=True, lim=lim)
    #     else:
    #         print(f'{test_save_file_name} exists already! Moving on...')
    #     # make appropriate designs
    #     if a == 'All':
    #         target_save_file_name = f'{out}/designs_{a}_by_Genus_2_universal.pickle'
    #     else:
    #         target_save_file_name = f'{out}/designs_{a}_Only_by_Genus_2_universal.pickle'
    #     if not os.path.exists(target_save_file_name):
    #         design_pickle_name = ribodesigner(target_sequences_folder=target, ref_sequence_file=ref_seq, guide_length=n,
    #                                           igs_length=m, min_length=minlen, folder_to_save=out,
    #                                           selective=False, min_true_cov=0, msa_fast=True, score_type='weighted',
    #                                           n_limit=1, percent_of_target_seqs_used=1, gaps_allowed=False, fileout=False,
    #                                           random_guide_sample_size=10)
    #     else:
    #         print(f'{target_save_file_name} exists already! Moving on...')
    #     _ = couple_designs_to_test_seqs(designs_input=target_save_file_name, test_seqs_input=test_save_file_name,
    #                                     flexible_igs=True, file_to_save=out)
    #     output = run_local(output_folder=out, guide_len=n)
    #
    # Generate test sequences with taxonomy
    # taxonomy_levels_all = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
    # to_generate = {'Phylum': ['Proteobacteria', 'Firmicutes', 'Actinobacteriota_and_Firmicutes'],
    #                'Class': ['Gammaproteobacteria', 'Bacilli'],
    #                'Order': ['Enterobacterales', 'Pseudomonadales', 'Bacillales'],
    #                'Family': ['Enterobacteriaceae', 'Pseudomonadaceae', 'Bacillaceae'],
    #                'Genus': ['Escherichia-Shigella', 'Pseudomonas', 'Bacillus']}
    # test_seq_path = f'Datasets_used/SILVA_squished_datasets_1_per_genus/Selective datasets per taxonomy/'
    # # test_seq_path = f'Datasets_used/SILVA_squished_datasets_5_per_genus/Selective datasets per taxonomy/'
    # save_file_path = output_path + 'selective_by_taxonomy/'
    #
    # # Graph!
    # test_seqs_to_process = [([test_seq_path + f'{taxonomy}_{include}_included/{taxonomy}_{include}_included_1.fasta',
    #                           test_seq_path + f'{taxonomy}_{include}_excluded/{taxonomy}_{include}_excluded_1.fasta'],
    #                          save_file_path + f'{taxonomy}_{include}')
    #                         for taxonomy in taxonomy_levels_all for include in to_generate[taxonomy]]
    # # test_seqs_to_process.reverse()
    # target_seqs_to_process = [[test_seq_path + f'{taxonomy}_{include}_included/{taxonomy}_{include}_included_2.fasta',
    #                           test_seq_path + f'{taxonomy}_{include}_excluded/{taxonomy}_{include}_excluded_2.fasta']
    #                         for taxonomy in taxonomy_levels_all for include in to_generate[taxonomy]]
    # for (test_files, out_path), target_files in zip(test_seqs_to_process, target_seqs_to_process):
    #     all_test_file_names = []
    #     all_target_file_names = []
    #     for file in test_files:
    #         title = file.split('.')[0].split('/')[-1]
    #         test_save_file_name = f'{out_path}/test_sequences_{title}.pickle'
    #         all_test_file_names.append(test_save_file_name)
    #         if not os.path.exists(test_save_file_name):
    #             test_seqs_pickle_file_name = prepare_test_seqs(test_folder=file, ref_sequence_file=ref_path,
    #                                                            guide_length=n, igs_length=m, min_length=minlen,
    #                                                            folder_to_save=out_path,
    #                                                            graph_results=True, var_regs=e_coli_var_regs,
    #                                                            graph_file_type='png', get_consensus_batches=True,
    #                                                            batch_num=10, score_type='weighted', msa_fast=True,
    #                                                            remove_x_dupes_in_graph=True, lim=1580)
    #         else:
    #             print(f'{test_save_file_name} exists already! Moving on...')
    #     for target_file in target_files:
    #         target_title = target_file.split('.')[0].split('/')[-1]
    #         target_save_file_name = f'{out_path}/designs_{target_title}_universal.pickle'
    #         all_target_file_names.append(target_save_file_name)
    #
    #         if not os.path.exists(target_save_file_name):
    #             design_pickle_name = ribodesigner(target_sequences_folder=target_file, ref_sequence_file=ref_path,
    #                                                   guide_length=n,
    #                                                   igs_length=m, min_length=minlen, fileout=False,
    #                                                   folder_to_save=out_path,
    #                                                   selective=False, min_true_cov=0, msa_fast=True, score_type='weighted',
    #                                                   n_limit=1, percent_of_target_seqs_used=1, gaps_allowed=False,
    #                                                   random_guide_sample_size=10)
    #         else:
    #             print(f'{target_save_file_name} exists already! Moving on...')
    #             continue
    #         for test_outfile in all_test_file_names:
    #             _ = couple_designs_to_test_seqs(designs_input=target_save_file_name, test_seqs_input=test_outfile,
    #                                             flexible_igs=True, file_to_save=out_path)
    #     output = run_local(output_folder=out_path, guide_len=n)
    # all_test_file_names = []
    # all_target_file_names = []
    # test_seq_path = f'Datasets_used/SILVA_squished_datasets_3000_per_order/'
    # test_seqs_to_process = [([test_seq_path + f'Order_{include}_included/Order_{include}_included_1.fasta'],
    #                          save_file_path + f'Order_{include}_all_seqs')
    #                         for include in ['Enterobacterales', 'Pseudomonadales']]
    # target_seqs_to_process = [([test_seq_path + f'Order_{include}_included/Order_{include}_included_2.fasta'],
    #                            save_file_path + f'Order_{include}_all_seqs')
    #                           for include in ['Enterobacterales', 'Pseudomonadales']]
    # for (test_files, out_path) in test_seqs_to_process:
    #     for file in test_files:
    #         title = file.split('.')[0].split('/')[-1]
    #         test_save_file_name = f'{out_path}/test_sequences_{title}.pickle'
    #         all_test_file_names.append(test_save_file_name)
    #         if not os.path.exists(test_save_file_name):
    #             test_seqs_pickle_file_name = prepare_test_seqs(test_folder=file, ref_sequence_file=ref_path,
    #                                                            guide_length=n, igs_length=m, min_length=minlen,
    #                                                            folder_to_save=out_path,
    #                                                            graph_results=True, var_regs=e_coli_var_regs,
    #                                                            graph_file_type='png', get_consensus_batches=True,
    #                                                            batch_num=10, score_type='weighted', msa_fast=True,
    #                                                            remove_x_dupes_in_graph=True, lim=1580)
    #         else:
    #             print(f'{test_save_file_name} exists already! Moving on...')
    #
    # for (target_files, out_path) in target_seqs_to_process:
    #     for target_file in target_files:
    #         target_title = target_file.split('.')[0].split('/')[-1]
    #         target_save_file_name = f'{out_path}/designs_{target_title}_universal.pickle'
    #         all_target_file_names.append(target_save_file_name)
    #
    #         if not os.path.exists(target_save_file_name):
    #             design_pickle_name = ribodesigner(target_sequences_folder=target_file, ref_sequence_file=ref_path,
    #                                                   guide_length=n,
    #                                                   igs_length=m, min_length=minlen, fileout=False,
    #                                                   folder_to_save=out_path,
    #                                                   selective=False, min_true_cov=0, msa_fast=True, score_type='weighted',
    #                                                   n_limit=1, percent_of_target_seqs_used=1, gaps_allowed=False,
    #                                                   random_guide_sample_size=10)
    #         else:
    #             print(f'{target_save_file_name} exists already! Moving on...')
    #         for test_outfile in all_test_file_names:
    #             _ = couple_designs_to_test_seqs(designs_input=target_save_file_name, test_seqs_input=test_outfile,
    #                                             flexible_igs=True, file_to_save=out_path)
    #     output = run_local(output_folder=out_path, guide_len=n)
    #
    # graphs_multiple_conditions(universal_path='test_output_files/universal_diff_var_regs',
    #                            selective_path='test_output_files/selective_by_taxonomy',
    #                            output_folder='test_output_files/best_designs/For experimental selection',
    #                            add_overhangs=True, m_smithii_var_regs=m_smithii_var_regs)

    # # fungiii
    # test_seqs_to_process = ['Datasets_used/SILVA_squished_datasets_fungi/Saccharomyces cerevisiae_Only_by_Species_test.fasta',
    #                         'Datasets_used/SILVA_squished_datasets_fungi/Ascomycota_Basidiomycota_Only_by_Family_unique_Species_test.fasta']
    # target_seqs_to_process = ['Datasets_used/SILVA_squished_datasets_fungi/Ascomycota_Basidiomycota_Only_by_Family_unique_Species_target.fasta']
    # out_path = 'test_output_files/fungi'
    # all_test_file_names = []
    # all_target_file_names = []
    # for test_file in test_seqs_to_process:
    #     title = test_file.split('.')[0].split('/')[-1]
    #     test_save_file_name = f'{out_path}/test_sequences_{title}.pickle'
    #     all_test_file_names.append(test_save_file_name)
    #     if not os.path.exists(test_save_file_name):
    #         test_seqs_pickle_file_name = prepare_test_seqs(test_folder=test_file, ref_sequence_file=ref_path_euk,
    #                                                        guide_length=n, igs_length=m, min_length=n,
    #                                                        folder_to_save=out_path, graph_results=True,
    #                                                        var_regs=s_cerevisiae_var_regs, graph_file_type='png',
    #                                                        get_consensus_batches=True, batch_num=10,
    #                                                        score_type='weighted', msa_fast=True,
    #                                                        remove_x_dupes_in_graph=True, lim=1800)
    #     else:
    #         print(f'{test_save_file_name} exists already! Moving on...')
    # for target_file in target_seqs_to_process:
    #     target_title = target_file.split('.')[0].split('/')[-1]
    #     target_save_file_name = f'{out_path}/designs_{target_title}_universal.pickle'
    #     all_target_file_names.append(target_save_file_name)
    #
    #     if not os.path.exists(target_save_file_name):
    #         design_pickle_name = ribodesigner(target_sequences_folder=target_file, ref_sequence_file=ref_path_euk,
    #                                           guide_length=n, igs_length=m, min_length=n,
    #                                           folder_to_save=out_path, min_true_cov=0, msa_fast=True,
    #                                           score_type='weighted', percent_of_target_seqs_used=1,
    #                                           gaps_allowed=False, random_guide_sample_size=10)
    #     else:
    #         print(f'{target_save_file_name} exists already! Moving on...')
    #
    #     for test_outfile in all_test_file_names:
    #         _ = couple_designs_to_test_seqs(designs_input=target_save_file_name, test_seqs_input=test_outfile,
    #                                         flexible_igs=True, file_to_save=out_path)


    # Here is an example of how to run ribodesigner with all options to make designs locally and test them locally.
    test_seqs_to_process = ['Datasets_used/SILVA_squished_datasets_fungi/Saccharomyces cerevisiae_Only_by_Species_test.fasta',
                            'Datasets_used/SILVA_squished_datasets_fungi/Ascomycota_Basidiomycota_Only_by_Family_unique_Species_test.fasta']
    target_seqs_to_process = ['Datasets_used/SILVA_squished_datasets_fungi/Ascomycota_Basidiomycota_Only_by_Family_unique_Species_target.fasta']
    out_path = 'test_output_files/fungi'

    # If when generating graphs it says it is not done cooking, comment out the ribodesigner routine and just run
    # run_local or else you risk having to couple everything again
    target_file_names, test_file_names = (
        ribodesigner_routine(target_seqs_to_process=target_seqs_to_process, test_seqs_to_process=test_seqs_to_process,
                             out_path=out_path, ref_seq_file=ref_path_euk, guide_len=n, igs_len=m, min_len=n,
                             graph_results=True,var_regs=s_cerevisiae_var_regs, graph_type='png',
                             get_consensus_batches= True, batch_num=10, score_type='weighted', msa_fast=True,
                             remove_x_dupes_in_graph=True, var_regs_lim=1800, min_true_cov=0,
                             percent_of_target_seqs_used=1, gaps_allowed=False, random_guide_sample_size=10,
                             flexible_igs=True))

    # If it says that big_checkpoint is corrupted, follow the prompts on the command line and run everything again!
    # It will skip analyzing files that have already been made and remake the checkpoint file if needed.
    output = run_local(output_folder=out_path, guide_len=n, num_of_workers=number_of_workers)

    # # if run_remotely:
    # #     run_remote(out_path, guide_len, n_limit=1, scratch_path='/scratch/kpr1/RiboDesigner/',
    # #                number_of_workers=mp.cpu_count(), worker_number=0)
    # # else:
    # #     output = run_local(output_folder=out_path, guide_len=guide_len)
    #
    stringent_path = 'test_output_files/fungi/coupled/results/combined/75016_designs_designs_Ascomycota_Basidiomycota_Only_by_Family_unique_Species_target_universal_vs_test_sequences_Saccharomyces cerevisiae_Only_by_Species_test.txt'
    to_compare_path = 'test_output_files/fungi/coupled/results/combined/75016_designs_designs_Ascomycota_Basidiomycota_Only_by_Family_unique_Species_target_universal_vs_test_sequences_Ascomycota_Basidiomycota_Only_by_Family_unique_Species_test.txt'
    get_fungi_designs(results_stringent_path=stringent_path, results_to_compare_path=to_compare_path,
                 output_folder=out_path, file_type='png', add_overhangs=True)

    # # For SEED poster
    # graphs_multiple_conditions(universal_path='test_output_files_v1/universal_diff_var_regs',
    #                            selective_path='test_output_files_v1/selective_by_taxonomy',
    #                            output_folder='/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/'
    #                                          'My Drive/KRG Thesis/Posters/SEED 2024/Data files from Ribodesigner',
    #                            add_overhangs=True, m_smithii_var_regs=m_smithii_var_regs)
