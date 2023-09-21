from playsound import playsound
from ribodesigner import ribodesigner, test_ribo_design, make_graphs

if __name__ == '__main__':
    # Figure 2: synthetic community data violin plots
    # Basically, look at the datasets used data and then run RiboDesigner on them.
    # Generate a violin plot using the results.
    # Run RiboDesigner on all datasets we are looking at
    # First, set up base files and parameters
    m = 5
    n = 50
    minlen = 35

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
    good_targets = "/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/My Drive/KRG Thesis/Scripts/" \
                   "Data files and Outputs/Ribozyme paper dataset/Original files"
    bad_targets = 'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Species/Bacillus_halotolerans.fasta'
    universal_data_1 = 'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_2.fasta'
    universal_data_2 = 'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_3.fasta'
    universal_data_3 = 'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_4.fasta'
    big_data_entero_only = 'Datasets_used/SILVA_squished_datasets/Enterobacterales_only_squished/Enterobacterales_1.fasta'
    big_data_pseudo_only = 'Datasets_used/SILVA_squished_datasets/Pseudomonadales_only_squished/Pseudomonadales_1.fasta'
    big_data_no_pseudo = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished/Background_Bacteria_squished_no_pseudo.fasta'
    big_data_no_entero = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished/Background_Bacteria_squished_no_entero.fasta'
    big_data_background = 'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_Bacteria_Only/Bacteria_Only_by_Genus_1.fasta'
    test_output_folder = 'test_output_files/test_outputs_ribodesigner_v2'
    test_file = 'test_dataset_for_graphs.csv'
    big_data_file_for_output = 'large_dataset.csv'

    # Test new RiboDesigner for images
    universal_datasets = []
    selective_datasets = []

    # control_design = test_ribo_design(design=u1376, target_folder=big_data_background, ref_seq_folder=ref_path, igs_len=m,
    #                                   score_type='weighted', thresh=0.5, msa_fast=True, gaps_allowed=False,
    #                                   file_out=True, folder_to_save=test_output_folder + f'/control dataset')

    for i, dataset in enumerate([universal_data_1, universal_data_2, universal_data_3]):
        out_data_temp = ribodesigner(target_sequences_folder=dataset, ref_sequence_file=ref_path, igs_length=m,
                                     guide_length=n, min_length=n, selective=False, min_true_cov=0.3,
                                     background_sequences_folder=big_data_background, identity_thresh=0.5,
                                     msa_fast=True, percent_of_background_seqs_used=1, score_type='weighted', n_limit=0,
                                     percent_of_target_seqs_used=0.01, gaps_allowed=False, fileout=True,
                                     random_guide_sample_size=3,
                                     folder_to_save=test_output_folder + f'/universal dataset {i + 1}')
        universal_datasets.append(out_data_temp)

    for i, datasets in enumerate([(big_data_entero_only, big_data_no_entero), (big_data_pseudo_only, big_data_no_pseudo), big_data_3]):
        out_data_temp = ribodesigner(target_sequences_folder=datasets[0], ref_sequence_file=ref_path, igs_length=m,
                                     guide_length=n, min_length=n, selective=True, min_true_cov=0.3,
                                     background_sequences_folder=datasets[1], identity_thresh=0.5,
                                     msa_fast=True, percent_of_background_seqs_used=1, score_type='weighted', n_limit=0,
                                     percent_of_target_seqs_used=1, gaps_allowed=False, fileout=True,
                                     random_guide_sample_size=10,
                                     folder_to_save=test_output_folder + f'/selective dataset {i + 1}')
        selective_datasets.append(out_data_temp)

    make_graphs(control_designs=control_design, selective_designs=selective_datasets,
                universal_designs=universal_datasets, var_regs=e_coli_var_regs, file_loc=test_output_folder + '/' + big_data_file_for_output)


    # # This is using the csv made with the code on top of this one
    # make_graphs(control_designs=[], selective_designs=[],
    #             universal_designs=[], var_regs=e_coli_var_regs, data_file=test_output_folder + '/' + test_file)

    # This is using the csv made with the code on top of this one
    make_graphs(control_designs=[], selective_designs=[],
                universal_designs=[], var_regs=e_coli_var_regs, data_file=test_output_folder + '/' + big_data_file_for_output)

    playsound('/System/Library/Sounds/Pop.aiff')
    print(f'Test data done!\n########################################################\n')