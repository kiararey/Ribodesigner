import multiprocessing as mp
import sys
from ribodesigner import (import_data_to_df, run_local, run_remote, ribodesigner_routine)

if __name__ == '__main__':
    # Here is an example of how to run ribodesigner with all options to make designs locally and test them locally.

    # --------- Edit everything between these lines so that your correct sequences are being analyzed!
    # a list of .fasta files including organisms that you want to make designs for
    target_seqs_to_process = ['<Path to your target sequences file>', '<Another path to your target sequences file>']
    # a list of .fasta files including organsims that you want to test your designs against
    test_seqs_to_process = ['<Path to your test sequences file>', '<Another path to your test sequences file>']
    # the path to the folder where you want your files to be saved
    out_path = '<path to save your files>'
    # the path to the refrence sequence you will use to align everything to
    ref_path = '<Path to your reference fasta file>'
    # ---------

    # m is your IGS/ pentanucleotide length
    m = 5
    # n is your ideal guide length
    n = 50
    # minlen is your minimum tolerated guide length. We kept it in our study as 50
    minlen = n

    # If you are using a cluster, set this to true
    run_remotely = False

    try:
        worker_number = sys.argv[1]
        number_of_workers = sys.argv[2]
    except:
        worker_number = 0
        number_of_workers = mp.cpu_count()
        pass

    # this is from Reich, M. & Labes, A. How to boost marine fungal research: A first step towards a multidisciplinary
    # approach by combining molecular fungal ecology and natural products chemistry. Marine Genomics 36, 57-75 (2017).
    s_cerevisiae_var_regs = [(69, 80), (126, 292), (478, 510), (643, 850), (1048, 1070), (1350, 1400), (1480, 1531),
                             (1674, 1730)]

    # Chakravorty, S., Helb, D., Burday, M., Connell, N. & Alland, D. A detailed analysis of 16S ribosomal RNA gene
    # segments for the diagnosis of pathogenic bacteria. J Microbiol Methods 69, 330-339 (2007).
    e_coli_var_regs = [(69, 99), (137, 242), (433, 497), (576, 682), (822, 879), (986, 1043), (1117, 1173),
                       (1243, 1294), (1435, 1465)]

    # # I got archaeal 16s variable regions by aligning to e_coli variable regions. This way we can graph anything
    # # aligned to m_smithii with its own index
    # m_smithii_var_regs = adjust_var_regs(known_seq_file=ref_path, known_var_regs=e_coli_var_regs,
    #                                      unknown_seq_file=ref_path_arc)
    # This should result in the following variable regions:
    m_smithii_var_regs = [(64, 75), (113, 223), (389, 437), (514, 622), (760, 822), (923, 982), (1052, 1119),
                          (1190, 1241), (1388, 1413)]

    var_regs_dict = {'All': e_coli_var_regs, 'Bacteria': e_coli_var_regs, 'Eukaryota': s_cerevisiae_var_regs,
                     'Archaea': m_smithii_var_regs}

    # If when generating graphs it says it is not done cooking, comment out the ribodesigner routine and just run
    # run_local or else you risk having to couple everything again
    
    # For an explanation of what each of these options does, please check out the ribodesigner_routine function in
    # ribodesigner.py as it has explanations on each one.
    target_file_names, test_file_names = (
        ribodesigner_routine(target_seqs_to_process=target_seqs_to_process, test_seqs_to_process=test_seqs_to_process,
                             out_path=out_path, ref_seq_file=ref_path, guide_len=n, igs_len=m, min_len=minlen,
                             graph_results=True, var_regs=var_regs_dict['Bacteria'], graph_type='png',
                             get_consensus_batches_test=True, get_consensus_batches_designs=False, batch_num=10,
                             score_type='weighted', msa_fast=True, var_regs_lim=1800,
                             min_true_cov=0, percent_of_target_seqs_used=1,
                             random_guide_sample_size=10, flexible_igs=True))

    # If it says that big_checkpoint is corrupted, follow the prompts on the command line and run everything again!
    # It will skip analyzing files that have already been made and remake the checkpoint file if needed.
    if run_remotely:
        run_remote(out_path, n, n_limit=1, scratch_path=out_path,
                   number_of_workers=mp.cpu_count(), worker_number=number_of_workers)
    else:
        output, output_files = run_local(output_folder=out_path, guide_len=n, num_of_workers=number_of_workers)
        df_data = import_data_to_df(output_files)
        # This will save all of your data as a single .csv file. I recommend using pandas to analyze the df though.
        df_data.to_csv(out_path + 'results.csv')