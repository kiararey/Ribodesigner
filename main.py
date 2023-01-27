# This is a sample Python script.
from RiboDesigner import RiboDesigner
import figure_plotting_functions as fig_plt
import os
import numpy as np
from RiboDesigner import RiboDesigner
import time

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
class OptimizedRibozyme:
    def __init__(self, full_design: str, igs: str, optimized_guide_g_igs: str, ref_idx: int, perc_cov: float,
                 perc_on_target:float, true_perc_cov: float, barcode: str, targets: list):
        self.seq = full_design
        self.igs = igs
        self.guide_g_igs = optimized_guide_g_igs
        self.ref_idx = ref_idx
        self.perc_cov = perc_cov
        self.perc_on_target = perc_on_target
        self.true_perc_cov = true_perc_cov
        self.barcode = barcode
        self.targets = targets

    def __str__(self):
        return f'{self.seq})'

    def __repr__(self):
        return f'{self.seq})'


class BasicRibozyme:
    def __init__(self, full_design: str, igs: str, guide_g_igs: str, ref_idx: int, perc_cov: float,
                 perc_on_target:float, true_perc_cov: float, barcode: str, target: str):
        self.seq = full_design
        self.igs = igs
        self.guide_g_igs = guide_g_igs
        self.ref_idx = ref_idx
        self.perc_cov = perc_cov
        self.perc_on_target = perc_on_target
        self.true_perc_cov = true_perc_cov
        self.barcode = barcode
        self.target = target

    def __str__(self):
        return f'{self.seq})'

    def __repr__(self):
        return f'{self.seq})'


class TargetSeq:
    def __init__(self, name: str, sequ: str, igs_and_ribo_seqs: list, igs_and_guide_seqs: list, all_igs: list, guides: list, indexes: list):
        self.name = name
        self.sequ = sequ
        self.igs_and_ribo_seqs = igs_and_ribo_seqs
        self.igs_and_guide_seqs = igs_and_guide_seqs
        self.all_igs = all_igs
        self.guides = guides
        self.indexes = indexes

    def __str__(self):
        return f'{self.name}({self.seq})'

    def __repr__(self):
        return f'{self.name}({self.seq})'

#
# class AlignedTarget:
#     [name, sequ, (temp_ribozyme_sequences, temp_IGS_and_guide_sequences, temp_IGSes,
#                   temp_guide_sequences, temp_og_and_ref_idexes)]


# Press the green button in the gutter to run the script.
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
    output_path = 'Figure_2_output_files/'

    # Reference sequence - will have to be E. coli to graph the variable regions
    ref_path = 'Common_sequences/e-coli-16s-mg1655.fasta'

    ########################################################
    start = time.time()
    # Test datasets
    print('Running test data...\n')
    test_data = "/Users/kiarareyes/Library/CloudStorage/GoogleDrive-kpr1@rice.edu/My Drive/KRG Thesis/Scripts/" \
                "Data files and Outputs/Ribozyme paper dataset/Original files"

    out_data = RiboDesigner(m, n, minlen, barcode_seq_file, ribobody_file, test_data,
                                    min_true_cov=0, identity_thresh=0.7, fileout=False, ref_sequence_file=ref_path,
                                        folder_to_save=output_path)
    end = time.time()
    print(f'Test data done! Time taken: {end - start}s\n')


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

    ########################################################
    # SILVA squished datasets
    dataset_names = ['Archaea_Only', 'Eukaryota_Only', 'Bacteria_Only', 'All']
    # dataset_names = ['Bacteria_Only', 'All']
    output_path = 'SILVA_figure_output_files/'

    for name in dataset_names:
        datasets_path = f'Datasets_used/SILVA_squished_datasets/SILVA_squished_datasets_{name}/'
        print(f'Now analyzing data in {datasets_path[:-1]}...')

        datasets = np.array([file_name for file_name in os.listdir(datasets_path) if file_name != '.DS_Store'])
        datasets.sort()
        ribodesigner_settings = [m, n, minlen, barcode_seq_file, ribobody_file, 0, 0.7, True]
        start = time.time()
        fig_plt.score_vs_true_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)
        end = time.time()
        print(f'{name} figures done! Time taken: {end - start}s\n')
        fig_plt.plot_for_16s_coverage(datasets, datasets_path, output_path, ribodesigner_settings, ref_path)