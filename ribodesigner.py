import glob
import pickle
import multiprocessing as mp
import os
import subprocess
import time
from collections import defaultdict
from math import exp
import json

import Bio.motifs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import AlignIO
from Bio.Align import AlignInfo, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from alive_progress import alive_bar
from icecream import ic

ic.configureOutput(prefix='debug | --> ')


# Expanded version of SilvaSequence
class TargetSeq:
    # Important: TargetSeq id must be formatted like a Silva sequence if we want the give_taxonomy function to work!
    def __init__(self, id_attr: str = '', full_name_attr: str = '', seq_attr: Seq = '',
                 cat_sites_attr: np.ndarray = None):
        self.id = id_attr
        self.full_name = full_name_attr
        self.seq = seq_attr
        self.cat_sites = cat_sites_attr

    def __str__(self):
        return f'{self.id}({self.seq})'

    def __repr__(self):
        return f'{self.id}, cat_sites={self.cat_sites}'

    # noinspection DuplicatedCode
    def give_taxonomy(self, level: str):
        level_names = {'Domain': 0, 'Phylum': 1, 'Class': 2, 'Order': 3, 'Family': 4,
                       'Genus': 5, 'Species': 6, 'Taxon': 7}
        if level not in level_names:
            print('Taxonomic level not found')
        # gives us the name at a certain taxonomic level
        all_levels = self.full_name.split(';')
        try:
            here_is_taxonomy = all_levels[level_names[level]]
        except IndexError:
            here_is_taxonomy = ''
        return here_is_taxonomy

    def reset_cat_sites_and_seq(self):
        # To minimize data storage
        self.seq = ''
        self.cat_sites = []
        return self


class RibozymeDesign:
    """
    Accepted kwargs:
    'g': list[Seq]  # guides to use to later
    """

    def __init__(self, id_attr: str, guides_to_use_attr: list[Seq] = None, targets_attr: set = None,
                 guide_attr: Seq = '', score_attr: float = None, score_type_attr: str = '', perc_cov_attr: float = None,
                 perc_on_target_attr: float = None, true_perc_cov_attr: float = None,
                 background_score_attr: float = None, perc_cov_background_attr: float = None,
                 perc_on_target_background_attr: float = None, true_perc_cov_background_attr: float = None,
                 background_guides_attr: list[Seq] = None, anti_guide_attr: Seq = '', anti_guide_score_attr: int = None,
                 background_targets_attr: set = None, igs_attr: str = '', ref_idx_attr: int = None,
                 u_consv_background_attr: float = None, tested_design_attr: bool = False,
                 perc_cov_test_attr: float = None,
                 perc_on_target_test_attr: float = None, true_perc_cov_test_attr: float = None):
        self.id = id_attr
        if igs_attr:
            self.igs = igs_attr
        else:
            self.igs = id_attr[:5]
        if ref_idx_attr:
            self.ref_idx = ref_idx_attr
        else:
            self.ref_idx = int(id_attr[5:])
        self.guide = guide_attr
        self.guides_to_use = guides_to_use_attr
        self.targets = targets_attr
        if targets_attr:
            self.number_of_targets = len(targets_attr)
        else:
            self.number_of_targets = None
        self.score = score_attr
        self.score_type = score_type_attr

        # calculate percent coverage, percent on target, or true percent coverage if needed and possible.
        self.calc_percent_coverages(perc_cov_attr=perc_cov_attr, perc_on_target_attr=perc_on_target_attr,
                                    true_perc_cov_attr=true_perc_cov_attr)

        # If we initialized using an optimized design, set these attributes now.
        if self.score and self.true_perc_cov:
            self.composite_score = self.true_perc_cov * self.score
        else:
            self.composite_score = 0
        if self.guide and self.score and self.score_type:
            self.optimized_to_targets = True
        else:
            self.optimized_to_targets = False
        self.tested_design = tested_design_attr

        # For initialized targeted ribozyme designs only:
        self.background_score = background_score_attr
        self.background_targets = background_targets_attr
        if background_targets_attr:
            self.number_of_targets_background = len(background_targets_attr)
        else:
            self.number_of_targets_background = None
        self.background_guides = background_guides_attr
        self.u_conservation_background = u_consv_background_attr
        self.anti_guide = anti_guide_attr
        self.anti_guide_score = anti_guide_score_attr
        self.calc_background_percent_coverages(perc_cov_background_attr=perc_cov_background_attr,
                                               perc_on_target_background_attr=perc_on_target_background_attr,
                                               true_perc_cov_background_attr=true_perc_cov_background_attr)

        if self.true_perc_cov_background and background_score_attr:
            self.composite_background_score = self.true_perc_cov_background * background_score_attr
        else:
            self.composite_background_score = None

        if self.composite_score and self.composite_background_score:
            self.delta_igs_vs_background = self.true_perc_cov - self.true_perc_cov_background
            self.delta_guide_vs_background = self.score - self.background_score
            self.delta_vs_background = self.composite_score - self.composite_background_score
        elif self.composite_score:
            self.delta_igs_vs_background = self.true_perc_cov
            self.delta_guide_vs_background = self.score
            self.delta_vs_background = self.composite_score
        else:
            self.delta_igs_vs_background = 0
            self.delta_guide_vs_background = 0
            self.delta_vs_background = 0

        if self.composite_background_score:
            self.optimized_to_background = True
        else:
            self.optimized_to_background = False
        self.guide_batches = None
        self.test_targets = None
        self.u_conservation_test = None
        self.perc_cov_test = perc_cov_test_attr
        self.perc_on_target_test = perc_on_target_test_attr
        self.true_perc_cov_test = true_perc_cov_test_attr
        self.test_score = None
        self.composite_test_score = None

        if self.composite_score and self.composite_test_score:
            self.delta_igs_vs_test = self.true_perc_cov - self.true_perc_cov_test
            self.delta_guide_vs_test = self.score - self.test_score
            self.delta_vs_test = self.composite_score - self.composite_test_score
        elif self.composite_score:
            self.delta_igs_vs_test = self.true_perc_cov
            self.delta_guide_vs_test = self.score
            self.delta_vs_test = self.composite_score
        else:
            self.delta_igs_vs_test = None
            self.delta_guide_vs_test = None
            self.delta_vs_test = None
        self.tested = False
        self.number_of_targets_test = None
        self.name_of_test_dataset = None

    def to_dict(self, taxonomy: str | list = '', all_data: bool = False):
        inputs = {}
        if not taxonomy:
            taxonomy = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Taxon']
        if type(taxonomy) == list:
            for target in taxonomy:
                if self.targets:
                    target_input = (str(give_taxonomy(str(self.targets), level=target))
                                    .replace('[[', '{').replace('\']]', '}')
                                    .replace(', \'', ':').replace('\'], [', ';')
                                    .replace('\n', ''))
                    inputs[f'target_{target}'] = target_input
                else:
                    inputs[f'target_{target}'] = ''
                if self.background_targets:
                    back_target_input = (str(give_taxonomy(str(self.background_targets), level=target))
                                         .replace('[[', '{').replace('\']]', '}')
                                         .replace(', \'', ':').replace('\'], [', ';')
                                         .replace('\n', ''))
                    inputs[f'target_{target}_background'] = back_target_input
                else:
                    inputs[f'target_{target}_background'] = ''
                if self.test_targets:
                    test_target_input = (str(give_taxonomy(str(self.test_targets), level=target))
                                         .replace('[[', '{').replace('\']]', '}')
                                         .replace(', \'', ':').replace('\'], [', ';')
                                         .replace('\n', ''))
                    inputs[f'target_{target}_test'] = test_target_input
                else:
                    inputs[f'target_{target}_test'] = ''
        else:
            if self.targets:
                target_input = (str(give_taxonomy(str(self.targets), level=taxonomy))
                                .replace('[[', '{').replace('\']]', '}')
                                .replace(', \'', ':').replace('\'], [', ';')
                                .replace('\n', ''))
            else:
                target_input = ''
            if self.background_targets:
                background_target_input = (str(give_taxonomy(str(self.background_targets), level=taxonomy))
                                           .replace('[[', '{').replace('\']]', '}')
                                           .replace(', \'', ':').replace('\'], [', ';')
                                           .replace('\n', ''))
            else:
                background_target_input = ''
            if self.test_targets:
                test_target_input = (str(give_taxonomy(str(self.test_targets), level=taxonomy))
                                     .replace('[[', '{').replace('\']]', '}')
                                     .replace(', \'', ':').replace('\'], [', ';')
                                     .replace('\n', ''))
            else:
                test_target_input = ''
            inputs[taxonomy] = target_input
            inputs[f'{taxonomy}_background'] = background_target_input
            inputs[f'{taxonomy}_test'] = test_target_input
        return_dict = {
            'id': self.id,
            'igs': self.igs,
            'reference_idx': self.ref_idx,
            'optimized_to_targets': self.optimized_to_targets,
            'optimized_to_background': self.optimized_to_background,
            'tested': self.tested,
            'tested_design': self.tested_design,
            'guide': str(self.guide),
            'num_of_targets': self.number_of_targets,
            'score_type': self.score_type,
            'score': self.score,
            '%_coverage': self.perc_cov,
            '%_on target': self.perc_on_target,
            'true_%_cov': self.true_perc_cov,
            'composite_score': self.composite_score,
            'num_of_targets_background': self.number_of_targets_background,
            'u_conservation_background': self.u_conservation_background,
            'background_score': self.background_score,
            '%_coverage_background': self.perc_cov_background,
            '%_on target_background': self.perc_on_target_background,
            'true_%_cov_background': self.true_perc_cov_background,
            'composite_background_score': self.composite_background_score,
            'delta_igs_vs_background': self.delta_igs_vs_background,
            'delta_guide_vs_background': self.delta_guide_vs_background,
            'delta_vs_background': self.delta_vs_background,
            'name_of_test_dataset': self.name_of_test_dataset,
            'num_of_targets_test': self.number_of_targets_test,
            'u_conservation_test': self.u_conservation_test,
            'test_score': self.test_score,
            '%_coverage_test': self.perc_cov_test,
            '%_on target_test': self.perc_on_target_test,
            'true_%_cov_test': self.true_perc_cov_test,
            'composite_test_score': self.composite_test_score,
            'delta_igs_vs_test': self.delta_igs_vs_test,
            'delta_guide_vs_test': self.delta_guide_vs_test,
            'delta_vs_test': self.delta_vs_test

        }
        if all_data:
            guide_tuples = [f'({seq};{val})' for seq, val in self.guide_batches]
            guide_batches_input = ';'.join(guide_tuples)
            more_data = {
                'guides_to_use': self.guides_to_use,
                'background_guides': self.background_guides,
                'anti_guide': str(self.anti_guide),
                'anti_guide_score': self.anti_guide_score,
                'guide_batches': guide_batches_input
            }
            return_dict.update(more_data)
        return_dict.update(inputs)
        return return_dict

    def __lt__(self, other):
        if self.optimized_to_background:
            return self.delta_vs_background < other.delta_vs_background
        else:
            return self.composite_score < other.composite_score

    def __str__(self):
        # ID: igs + ref_idx, then also display the guide sequence
        return f'ID:{self.igs}{self.ref_idx}, Guide:{self.guide}'

    def __repr__(self):
        text_basic = f'{type(self).__name__} {self.id}: (IGS={self.igs}, ref_idx={self.ref_idx}, ' \
                     f'targets_to_use={self.number_of_targets}, perc_cov={self.perc_cov}, ' \
                     f'perc_on_target={self.perc_on_target}, true_perc_cov={self.true_perc_cov}, ' \
                     f'optimized_to_targets={self.optimized_to_targets}, ' \
                     f'optimized_to_background={self.optimized_to_background}'
        target_info = f'\nGuide_sequence={self.guide}, score={self.score}, score_type={self.score_type}, ' \
                      f'composite_score={self.composite_score}'
        background_info = f'\nbackground_score={self.background_score}, ' \
                          f'perc_cov_background={self.perc_cov_background}, ' \
                          f'perc_on_target_background={self.perc_on_target_background}, ' \
                          f'true_perc_cov_background={self.true_perc_cov_background}, ' \
                          f'composite_background_score={self.composite_background_score}, ' \
                          f'delta_composite_scores={self.delta_vs_background}'
        if self.optimized_to_targets and self.optimized_to_background:
            return text_basic + target_info + background_info + ')'
        elif self.optimized_to_targets:
            return text_basic + target_info + ')'
        elif self.optimized_to_background:
            return text_basic + background_info + ')'
        else:
            return text_basic + ')'

    def calc_percent_coverages(self, perc_cov_attr: float = None, true_perc_cov_attr: float = None,
                               perc_on_target_attr: float = None):
        self.perc_cov = perc_cov_attr
        self.true_perc_cov = true_perc_cov_attr
        self.perc_on_target = perc_on_target_attr
        if perc_cov_attr is not None and true_perc_cov_attr is not None and perc_on_target_attr is not None:
            # calculate percent coverage, percent on target, or true percent coverage if needed and possible.
            if not self.perc_cov and true_perc_cov_attr and perc_on_target_attr:
                self.perc_cov = true_perc_cov_attr / perc_on_target_attr

            if not self.perc_on_target and true_perc_cov_attr and perc_cov_attr:
                self.perc_on_target = true_perc_cov_attr / perc_cov_attr

            if not self.true_perc_cov and perc_on_target_attr and perc_cov_attr:
                self.true_perc_cov = perc_on_target_attr * perc_cov_attr

    def calc_background_percent_coverages(self, perc_cov_background_attr: float = None,
                                          true_perc_cov_background_attr: float = None,
                                          perc_on_target_background_attr: float = None):
        self.perc_cov_background = perc_cov_background_attr
        self.true_perc_cov_background = true_perc_cov_background_attr
        self.perc_on_target_background = perc_on_target_background_attr

        if perc_on_target_background_attr is not None and true_perc_cov_background_attr is not None and \
                perc_on_target_background_attr is not None:
            # calculate percent coverage, percent on target, or true percent coverage if needed and possible.
            if not self.perc_cov_background and true_perc_cov_background_attr and perc_on_target_background_attr:
                self.perc_cov_background = true_perc_cov_background_attr / perc_on_target_background_attr

            if not self.perc_on_target_background and true_perc_cov_background_attr and perc_cov_background_attr:
                self.perc_on_target_background = true_perc_cov_background_attr / perc_cov_background_attr

            if not self.true_perc_cov_background and perc_on_target_background_attr and perc_cov_background_attr:
                self.true_perc_cov_background = perc_on_target_background_attr * perc_cov_background_attr

    def update_after_optimizing(self, score_attr, guide_attr, score_type_attr, reset_guides: bool = False):
        self.score = score_attr
        self.guide = guide_attr
        if reset_guides:
            self.guides_to_use = None
            self.anti_guide = None
        self.score_type = score_type_attr
        self.optimized_to_targets = True
        self.composite_score = self.true_perc_cov * score_attr

    def update_to_background(self, background_score_attr: float, new_guide: Seq, new_score: float,
                             reset_guides: bool = False):
        self.optimized_to_background = True
        self.background_score = background_score_attr
        self.score = new_score
        if self.background_targets:
            self.number_of_targets_background = len(self.background_targets)
        elif self.background_guides:
            self.number_of_targets_background = len(self.background_guides)
        else:
            self.number_of_targets_background = 0
        if reset_guides:
            self.background_guides = None
        self.guide = new_guide
        self.composite_background_score = self.true_perc_cov_background * background_score_attr
        self.delta_igs_vs_background = self.true_perc_cov - self.true_perc_cov_background
        self.delta_guide_vs_background = self.score - self.background_score
        self.delta_vs_background = self.composite_score - self.composite_background_score

    def update_to_test(self, test_score_attr: float, name_of_test_dataset_attr: str):
        self.tested = True
        self.name_of_test_dataset = name_of_test_dataset_attr
        self.test_score = test_score_attr
        if self.test_targets:
            self.number_of_targets_test = len(self.test_targets)
        else:
            self.number_of_targets_test = 0
        self.composite_test_score = self.true_perc_cov_test * test_score_attr
        if self.true_perc_cov and self.true_perc_cov_test:
            self.delta_igs_vs_test = self.true_perc_cov - self.true_perc_cov_test
        if self.score and self.test_score:
            self.delta_guide_vs_test = self.score - self.test_score
        if self.composite_score and self.composite_test_score:
            self.delta_vs_test = self.composite_score - self.composite_test_score

    def print_attributes(self, taxonomy='Order'):
        if self.optimized_to_targets and self.optimized_to_background:
            text = f'{self.igs},{self.ref_idx},{self.number_of_targets},{self.score},{self.score_type},' \
                   f'{self.perc_cov},{self.perc_on_target},{self.true_perc_cov},{self.composite_score},' \
                   f'{self.background_score},{self.perc_cov_background},{self.perc_on_target_background},' \
                   f'{self.true_perc_cov_background},{self.composite_background_score},{self.delta_vs_background},' \
                   f'{self.guide},{self.guide}G{self.igs}\n'
        elif self.optimized_to_targets:
            text = f'{self.igs},{self.ref_idx},{self.number_of_targets},{self.score},{self.score_type},' \
                   f'{self.perc_cov},{self.perc_on_target},{self.true_perc_cov},{self.composite_score},' \
                   f'{self.guide},{self.guide}G{self.igs}\n'
        elif self.tested_design:
            background_target_input = str(give_taxonomy(str(self.background_targets),
                                                        level=taxonomy)).replace(',', ';').replace('\'', '')
            # first_row = 'IGS,Reference index,Number of tested species targeted,Score type,Design score,
            # Score on targets,% cov background,% on target background,true % cov background,Composite score background,
            # Optimized guide,Optimized guide + G + IGS,Ideal guide,Ideal guide + G + IGS\n'
            text = f'{self.igs},{self.ref_idx},{self.number_of_targets_background},{background_target_input},' \
                   f'{self.score_type},{self.score},{self.background_score},{self.perc_cov_background},' \
                   f'{self.perc_on_target_background},{self.true_perc_cov_background},' \
                   f'{self.composite_background_score},{self.guide},{self.guide}G{self.igs},{self.anti_guide},' \
                   f'{self.anti_guide}G{self.igs}\n'

        else:
            print('Please optimize to targets or optimize to targets then to background first.')
            return
        return text


def ribodesigner(target_sequences_folder: str, igs_length: int = 5,
                 guide_length: int = 50, min_length: int = 35, ref_sequence_file=None, selective: bool = False,
                 background_sequences_folder: str = '', min_delta: float = 0, min_true_cov: float = 0.7,
                 fileout: bool = False, folder_to_save: str = '', n_limit: float = 0.0,
                 score_type: str = 'weighted', msa_fast: bool = False, gaps_allowed: bool = True,
                 percent_of_target_seqs_used: float = 1.0, percent_of_background_seqs_used: float = 1,
                 random_guide_sample_size: int = 10, test_folders: list[str] = None, store_batch_results: bool = False):
    """Generates ribozyme designs to target a set of sequences.
    :param percent_of_background_seqs_used: In case background data is very large we can get a random sample of the
    sequences used without replacement
    :param percent_of_target_seqs_used: In case target data is very large we can get a random sample of the sequences
    used without replacement
    :param gaps_allowed:
    :param min_delta: when making targeted designs, disregard any that have a composite score less than this
    :param target_sequences_folder: folder containing all the sequences you want to design ribozymes for. Can be either
    a folder of fasta files or a single fasta file.
    :param igs_length: how long you want your IGS to be in base pairs. Default is 5 bp.
    :param guide_length: how long you want your guide sequence to be in base pairs. Default is 5 bp.
    :param min_length: minimum guide sequence length from 3' end. Must be smaller or equal to guide_length.
    Default is 35 bp. Ex: if you want your guide sequence to bind to at least 35 nt at the 3' end of the target
    sequence, set min_length = 35.
    :param ref_sequence_file:
    :param background_sequences_folder: folder containing all the sequences you do NOT want to design ribozymes for.
    Can be either a folder of fasta files or a single fasta file.
    :param min_true_cov: minimum percentage of targets you want to hit at a conserved location with a single optimized
    design. Default is 0.7 (70% of targets).
    :param identity_thresh: How much sequence identity do you want to use for the MSA. Only applicable for designs made
    without considering IUPAC ambiguity codes.
    :param fileout: whether we want a csv file output or not. Default is False.
    :param folder_to_save: the path where the folder we will save our outputs in if fileout = True
    :param score_type:
    :param msa_fast: whether to use super5 MUSCLE MSA or just regular MUSCLE MSA. Recommended for large datasets (over
    300 sequences) for faster data processing.
    :param n_limit: What percentage of Ns are acceptable in the final design? This is applied when doing targeted
    designs only for the time being.
    """

    start = time.perf_counter()
    target_names_and_seqs = read_silva_fasta(in_file=target_sequences_folder)

    try:
        target_names_and_seqs[0]
    # if we hit an index error we've run out of sequence and
    # should not add new residues
    except IndexError:
        print(f'No sequences found in {target_sequences_folder}. Please make sure your files are not empty!\n')
        return None

    if not msa_fast and target_names_and_seqs.size > 300:
        print(f'Consider setting msa_fast=True for large datasets over 300 sequences for faster data processing!!\n')

    if not ref_sequence_file:
        # If we do not have a reference sequence, just choose one randomly
        print('No reference sequence provided. Picking a random sequence as reference...')
        ref_name_and_seq = np.random.choice(target_names_and_seqs)
    else:
        ref_name_and_seq = read_fasta(ref_sequence_file)[0]

    print(f'Found {target_names_and_seqs.size} total target sequences to analyze.')

    if percent_of_target_seqs_used < 1:
        target_names_and_seqs = np.random.choice(target_names_and_seqs,
                                                 size=round(target_names_and_seqs.size * percent_of_target_seqs_used),
                                                 replace=False)
        print(f'Randomly sampling {target_names_and_seqs.size} sequences to analyze.\n')

    time1 = time.perf_counter()
    print(f'Now finding catalytic sites and aligning to reference {ref_name_and_seq[0].replace("_", " ")}')
    process_num = mp.cpu_count()
    # with alive_bar(unknown='fish', spinner='fishes') as bar:
    #     fn_align_to_ref = np.vectorize(align_to_ref, otypes=[TargetSeq],
    #                                    excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])
    #     # prepare data for starmap
    #     in_data = [(sub_array, ref_name_and_seq, igs_length, guide_length, min_length) for sub_array in
    #                np.array_split(target_names_and_seqs, (min(process_num, len(target_names_and_seqs))))]
    #     # multithread
    #     with mp.Pool(processes=len(in_data)) as pool:
    #         out_data = pool.starmap(fn_align_to_ref, in_data)
    #     target_names_and_seqs = np.concatenate([sub_array for sub_array in out_data])
    #     bar()

    with alive_bar(unknown='fish', spinner='fishes') as bar:
        fn_align_to_ref = np.vectorize(align_to_ref, otypes=[TargetSeq],
                                       excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])
        fn_align_to_ref(target_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length,
                        guide_length=guide_length, min_length=min_length)
        bar()
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'finding catalytic sites and indexing target sequences to reference')
    time1 = time.perf_counter()

    # first, filter through possible designs and get the ones that meet our threshold
    to_optimize = filter_igs_candidates(aligned_targets=target_names_and_seqs, min_true_cov=min_true_cov,
                                        igs_len=igs_length)
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='finding putative ribozyme sites')

    # Now, optimize all of the possible guide sequences for each IGS:
    print('Now optimizing guide sequences...')
    # # Test optimizing guide sequences
    # for i in range(0, 10):
    #     optimize_designs(to_optimize[i], score_type=score_type, msa_fast=msa_fast,
    #                      gaps_allowed=gaps_allowed, compare_to_background=False,
    #                      random_sample_size=random_guide_sample_size, random_sample=True)
    #     ic()
    # time1 = time.perf_counter()
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        # vectorize function. For the multiprocessing module I made it return the updated RibozymeDesign

        fn_msa = np.vectorize(optimize_designs, otypes=[RibozymeDesign],
                              excluded=['score_type', 'guide_len', 'msa_fast', 'gaps_allowed', 'compare_to_background',
                                        'random_sample_size', 'random_sample', 'store_batch_results'])
        # prepare data for starmap
        in_data = [(sub_array, score_type, guide_length, msa_fast, gaps_allowed, False,
                    random_guide_sample_size, True, store_batch_results)
                   for sub_array in np.array_split(to_optimize, process_num)]

        # multithread
        with mp.Pool(processes=process_num) as pool:
            out_data = pool.starmap(fn_msa, in_data)
        optimized_seqs = np.concatenate([sub_array for sub_array in out_data])
        bar()

    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='generating optimized designs')

    if selective:
        print('please wait till selective ribozymes are re-implemented')
        # background_names_and_seqs = read_silva_fasta(in_file=background_sequences_folder)
        #
        # print(f'Found {background_names_and_seqs.size} total background sequences to analyze.')
        #
        # if percent_of_background_seqs_used < 1:
        #     background_names_and_seqs = np.random.choice(
        #         background_names_and_seqs, size=round(background_names_and_seqs.size * percent_of_target_seqs_used),
        #         replace=False)
        #     print(f'Randomly sampling {target_names_and_seqs.size} sequences to analyze.\n')
        #
        # # first find the IGSes and locations of the background sequences. We do not want to hit these.
        # # Also align these to the reference sequence
        # with alive_bar(unknown='fish', spinner='fishes') as bar:
        #     fn_align_to_ref(background_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length,
        #                     guide_length=guide_length, min_length=min_length)
        #     time2 = time.perf_counter()
        #     bar()
        # round_convert_time(start=time1, end=time2, round_to=4,
        #                    task_timed=f'finding catalytic sites and indexing background sequences to reference')
        #
        # print('Now applying designed ribozymes with background sequences and getting statistics...')
        # optimized_seqs = ribo_checker(optimized_seqs, background_names_and_seqs, len(ref_name_and_seq[1]),
        #                               guide_length=guide_length, flexible_igs=True, n_limit=n_limit,
        #                               target_background=not selective, refine_selective=True, for_test=False)
        #
        # time2 = time.perf_counter()
        # round_convert_time(start=time1, end=time2, round_to=4, task_timed='comparing designs against background '
        #                                                                   'sequences')
        # # Remove any .fasta
        # pickle_file_name = target_sequences_folder.split('.')[0].split('/')[-1] + '_selective_vs_' + \
        #                    background_sequences_folder.split('.')[0].split('/')[-1]
    else:
        # Remove any .fasta
        pickle_file_name = target_sequences_folder.split('.')[0].split('/')[-1] + '_universal'
    with open(f'{folder_to_save}/designs_{pickle_file_name}.pickle', 'wb') as handle:
        pickle.dump(optimized_seqs, handle)

    if fileout:
        out_file = folder_to_save + '/designs_' + pickle_file_name
        write_output_file(designs=optimized_seqs, folder_path=out_file, all_data=store_batch_results)

    end = time.perf_counter()
    round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    print('########################################################\n')
    return f'{folder_to_save}/designs_{pickle_file_name}.pickle'


def prepare_test_seqs(test_folder, ref_sequence_file, guide_length, igs_length, min_length, folder_to_save):
    start = time.perf_counter()
    test_names_and_seqs = read_silva_fasta(in_file=test_folder)
    print(f'Found {test_names_and_seqs.size} total test sequences to analyze.')

    if not ref_sequence_file:
        # If we do not have a reference sequence, just choose one randomly
        print('No reference sequence provided. Picking a random sequence as reference...')
        ref_name_and_seq = np.random.choice(test_names_and_seqs)
    else:
        ref_name_and_seq = read_fasta(ref_sequence_file)[0]

    fn_align_to_ref = np.vectorize(align_to_ref, otypes=[TargetSeq],
                                   excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])

    # first find the IGSes and locations and also align these to the reference sequence
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        fn_align_to_ref(test_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length,
                        guide_length=guide_length, min_length=min_length)
        time2 = time.perf_counter()
        bar()
    round_convert_time(start=start, end=time2, round_to=4,
                       task_timed=f'finding catalytic sites and indexing background sequences to reference')

    time1 = time.perf_counter()
    test_seqs_dict = filter_igs_candidates(aligned_targets=test_names_and_seqs, min_true_cov=0, igs_len=igs_length,
                                           return_test_seqs=True)
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='finding putative ribozyme sites')

    # remove .fasta from any files for pickle naming
    pickle_file_name = test_folder.split('.')[0].split('/')[-1]
    with open(f'{folder_to_save}/test_sequences_{pickle_file_name}.pickle', 'wb') as handle:
        pickle.dump(test_seqs_dict, handle)

    end = time.perf_counter()
    round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    print('########################################################\n')
    return f'{folder_to_save}/test_sequences_{pickle_file_name}.pickle'


def back_transcribe_seq_file(seq_file: str) -> Seq:
    """
    Converts DNA sequences into RNA transcripts

    :param seq_file: .txt file containing one or more DNA sequences
    :return: Seq object
    """
    # seq_file must be a .txt file
    with open(seq_file) as f:
        for i in f:
            out_seq = Seq(i).upper().back_transcribe()
    return out_seq


def read_silva_fasta(in_file: str, file_type: str = 'fasta', exclude: list = None, include: list = None,
                     exclude_include_taxonomy_level='') -> np.ndarray[TargetSeq]:
    """
    Reads in a single .fasta file or several fasta files from a directory. These files are assumed to be formatted as
    in the SILVA database as it makes finding taxonomy much easier.
    Include and exclude are mutually exclusive, and the function will take into account exclude first if both are given.
    If the taxonomy level is not given include and exclude will be ignored.
    If your files are NOT formatted as SILVA database files, please use read_fasta() instead.

    :param in_file: file path to single .fasta file or directory of .fasta files
    :param file_type: file extension
    :param exclude: a list of strings of taxonomy names to avoid. All else will be included.
    :param include: a list of strings of taxonomy names to include. All else will be excluded.
    :param exclude_include_taxonomy_level: the level of taxonomy needed to include or exclude sequences.
    :return: A list of tuples formatted as (ID, sequence)
    """

    def parse_file(file_to_parse, exclude_list, include_list, taxonomy_level=''):
        out_seqs = []
        with open(file_to_parse) as f:
            for title, seq in SimpleFastaParser(f):
                seq_id = title.split(' ')[0]
                putative_sequence = TargetSeq(id_attr=seq_id, full_name_attr=title,
                                              seq_attr=Seq(seq).upper().back_transcribe())
                if taxonomy_level and exclude_list or include_list:
                    if exclude_list and putative_sequence.give_taxonomy(level=taxonomy_level) not in exclude_list:
                        out_seqs.append(putative_sequence)
                    elif include_list and putative_sequence.give_taxonomy(level=taxonomy_level) in include_list:
                        out_seqs.append(putative_sequence)
                else:
                    out_seqs.append(putative_sequence)
        return out_seqs

    if in_file[-6:] == f'.{file_type}':
        # filepath here is a fasta file
        # will return a list of TargetSeqs
        target_seqs_and_names = parse_file(in_file, exclude, include, exclude_include_taxonomy_level)
    else:
        target_seqs_and_names = []
        for filename in glob.glob(os.path.join(in_file, '*.' + file_type)):
            target_seqs_and_names.extend(parse_file(filename, exclude, include, exclude_include_taxonomy_level))
    return np.array(target_seqs_and_names, dtype=TargetSeq)


def read_fasta(in_file: str, file_type: str = 'fasta') -> np.ndarray[TargetSeq]:
    """
    Reads in a single .fasta file or several fasta files from a directory

    :param in_file: file path to single .fasta file or directory of .fasta files
    :param file_type: file extension
    :return: A list of tuples formatted as (ID, sequence)
    """
    target_seqs_and_names = []
    if in_file[-6:] == f'.{file_type}':
        # filepath here is a fasta file
        # will return a list of tuples as (ID, sequence)
        with open(in_file) as f:
            for title, seq in SimpleFastaParser(f):
                target_seqs_and_names.append((title, Seq(seq).upper().back_transcribe()))
    else:
        for filename in glob.glob(os.path.join(in_file, '*.' + file_type)):
            with open(filename) as f:
                for title, seq in SimpleFastaParser(f):
                    target_seqs_and_names.append((title, Seq(seq).upper().back_transcribe()))
    return np.array(target_seqs_and_names, dtype=TargetSeq)


def round_convert_time(start: float, end: float, round_to: int = 4, task_timed: str = ''):
    time_to_round = abs(end - start)

    # Hours
    hours = time_to_round / 3600
    # Min
    min_time = hours % 1 * 60
    # Sec
    sec = min_time % 1 * 60

    if hours > 1:
        text = f'Time taken {task_timed}: {int(hours)} hrs, {int(min_time)} min, {round(sec, round_to)} sec\n'
    elif min_time > 1:
        text = f'Time taken {task_timed}: {int(min_time)} min, {round(sec, round_to)} sec\n'
    else:
        text = f'Time taken {task_timed}: {round(time_to_round, round_to)} sec\n'

    print(text)
    return


def find_cat_sites(target_sequence: TargetSeq, igs_length: int = 5, guide_length: int = 50,
                   min_length: int = 35):
    """
    Finds all instances of a U or T in a set of sequences and creates ribozyme designs for these sites.

    :param igs_length: desired IGS sequence length.
    :param guide_length: desired guide binding sequence length.
    :param min_length: minimum guide sequence length from 3' end. Must be smaller than guide_length.
    :return:
    """

    def get_data(i, target_sequence_data: TargetSeq):
        # Implement like:
        # fn_get_data = np.vectorize(get_data, otypes=[int], excluded=target_sequence)
        # data = fn_get_data(idx, target_sequence=target_sequence)
        # target_sequence.cat_sites = [(i, igs, guide) for i, igs, guide in zip(data[0], data[1], data[2])]
        guide = target_sequence_data.seq[i + 1:i + guide_length + 1].reverse_complement()
        # generate complementary IGS sequence igs_length bases long *upstream* of U site
        igs = target_sequence_data.seq[i - igs_length:i].reverse_complement()
        # Store as (idx, igs, guide)
        return i + 1, igs, guide

    if target_sequence.cat_sites:
        if len(target_sequence.cat_sites[0]) < 4:
            print('TargetSeq object already has catalytic sites.')
        else:
            print('TargetSeq object was already aligned and has catalytic sites.')
        return

    #   (a.k.a. minimum guide sequence length from 3' end) ex: if you want your guide sequence to
    #   bind to at least 35 nt at the 3' end of the target sequence, set min_length = 35.

    # run function to find the index of all U residues of each sequence in target_seqs
    # find all possible splice sites
    idx = np.array([i for i, ltr in enumerate(target_sequence.seq) if ltr == 'T'])
    # remove indexes that are < guide_length or > len(sequ) - igs_length (must have enough residues to attach to)
    min_idx = np.searchsorted(idx, igs_length) + 1
    max_idx = np.searchsorted(idx, len(target_sequence.seq) - min_length, side='right')
    idx = idx[min_idx:max_idx]

    if not idx.any():
        print(f'No viable catalytic sites in {target_sequence.id}')
        target_sequence.cat_sites = None
    else:
        fn_get_data = np.vectorize(get_data, otypes=[tuple], excluded=[target_sequence])
        data = fn_get_data(idx, target_sequence_data=target_sequence)
        target_sequence.cat_sites = data
    return


def align_to_ref(target_sequence: TargetSeq, ref_name_and_seq, igs_length: int = 5, guide_length: int = 50,
                 min_length: int = 35):
    """
    :param ref_name_and_seq: TargetSequence object
    :return:
    """
    if not target_sequence.cat_sites:
        find_cat_sites(target_sequence, igs_length, guide_length, min_length)
        if not target_sequence.cat_sites.any():
            return
    elif len(target_sequence.cat_sites[0]) == 4:
        print(f'TargetSeq object {target_sequence.id} has already been aligned to a reference sequence.')
        return

    # will have to keep in mind the potential lengths of the sequences and add length igs_length to our final
    # E.coli index
    aligner = PairwiseAligner(mode='global')
    alignments = aligner.align(target_sequence.seq, ref_name_and_seq[1])[0]

    # seq_a is the test sequence, seq_b is the reference sequence
    seq_a, seq_b = alignments
    # seq_a is the test sequence, seq_b is the reference sequence

    # obtain index of new U. Adding 1 because python range end is exclusive
    idx_seq_a = [i + 1 for i, ltr in enumerate(seq_a) if ltr == 'T']

    data = []
    current_spot = 0
    for og_idx, igs, guide in target_sequence.cat_sites:
        # Use our initial guess to see if we hit the current catalytic site
        seq_a_idx = len(seq_a[:idx_seq_a[current_spot]].replace('-', ''))
        while seq_a_idx != og_idx:
            # continue guessing until we finally hit the current catalytic site in our target sequence
            current_spot += 1
            seq_a_idx = len(seq_a[:idx_seq_a[current_spot]].replace('-', ''))
        # find what index that is based on the reference sequence by cuttingg the Seq at that spot...
        ref_string = seq_b[:idx_seq_a[current_spot]]
        # ... then removing all gaps to see what the actual location is on the reference sequence
        ref_idx = len(ref_string.replace('-', ''))
        data.append((og_idx, ref_idx, igs, guide))

    target_sequence.cat_sites = data
    return target_sequence


def filter_igs_candidates(aligned_targets: np.ndarray[TargetSeq], min_true_cov: float = 0, igs_len: float = 5,
                          return_test_seqs: bool = False) -> np.ndarray[RibozymeDesign]:
    """
    returns a dictionary of possible designs that meet our needs: checks each target sequence and filters the IGSes that
    work for the most targets. Returns a dictionary for each IGSid with the list of guides it needs to optimize, as
    well as the true percent coverage of these and the organism number that IGS exists in (this is the index in
    aligned_targets and we can use this later to calculate percent coverage and get taxonomy information of which
    sequences a particular target can hit.
    :param igs_len:
    :param aligned_targets:
    :param min_true_cov:
    :return:
    """

    print('Finding IGS sequences...')

    all_igs_ids = []
    counts_of_igs = defaultdict(int)
    counts_of_ref_idx = defaultdict(int)

    # Extract all the IGS id numbers - that's the IGS sequence and the reference index number
    with alive_bar(aligned_targets.size, spinner='fishes') as bar:
        for seq in aligned_targets:
            igses = set()
            ref_ids = set()
            for item in seq.cat_sites:
                all_igs_ids.append(f'{item[2]}{item[1]}')
                igses.add(f'{item[2]}')
                ref_ids.add(item[1])

            for igs in igses:
                if counts_of_igs[igs]:
                    counts_of_igs[igs] += 1
                else:
                    counts_of_igs[igs] = 1
            for ref_id in ref_ids:
                if counts_of_ref_idx[ref_id]:
                    counts_of_ref_idx[ref_id] += 1
                else:
                    counts_of_ref_idx[ref_id] = 1
            bar()

    # Measure true percent coverage of these and keep IGSes that are at least the minimum true percent coverage needed
    igs_ids_counted, igs_ids_counts = np.unique(all_igs_ids, return_counts=True)
    # get all the guides that have good IGSes
    if return_test_seqs:
        # IGS_id: [guides, targets, perc cov, perc on target, true perc cov]
        igs_over_min_true_cov = {igs: [[], set(), counts_of_igs[igs[:igs_len]] / aligned_targets.size,
                                       counts / counts_of_igs[igs[:igs_len]], counts / aligned_targets.size]
                                 for igs, counts in zip(igs_ids_counted, igs_ids_counts)}
        # Here each item in the list is an IGS id (IGS + reference index), a guide, and the location of the target sequence
        # in our initial aligned_targets array
        igs_subsets = [(f'{item[2]}{item[1]}', item[3], i) for i, seq in enumerate(aligned_targets) for item in
                       seq.cat_sites]
        print('Getting stats for putative target locations...')
        output_dict = defaultdict(dict)
        with alive_bar(len(igs_subsets), spinner='fishes') as bar:
            for igs_id, guide, target_num in igs_subsets:
                if igs_id in igs_over_min_true_cov:
                    # Add guide to list
                    igs_over_min_true_cov[igs_id][0].append(guide)
                    # Add organism number to set. We can calculate the percent coverage with this number later on.
                    igs_over_min_true_cov[igs_id][1].add(aligned_targets[target_num].full_name)
                    output_dict[int(igs_id[igs_len:])][igs_id] = igs_over_min_true_cov[igs_id]
                bar()

        # Now organize by ref_idx:
        for ref_id, count_of_ref in counts_of_ref_idx.items():
            output_dict[ref_id] = (output_dict[ref_id], count_of_ref)
        return output_dict
    else:
        # this gives us a dictionary where the ID is matched to the true percent coverage
        igs_over_min_true_cov = {igs: [[], set(), counts / aligned_targets.size] for igs, counts
                                 in zip(igs_ids_counted, igs_ids_counts)
                                 if counts / aligned_targets.size >= min_true_cov}

        # Here each item in the list is an IGS id (IGS + reference index), a guide, and the location of the target sequence
        # in our initial aligned_targets array
        igs_subsets = [(f'{item[2]}{item[1]}', item[3], i) for i, seq in enumerate(aligned_targets) for item in
                       seq.cat_sites]

        print(f'{len(igs_over_min_true_cov)} putative designs found.')

        print('Filtering IGSes that meet our min true % coverage...')
        with alive_bar(len(igs_subsets), spinner='fishes') as bar:
            for igs_id, guide, target_num in igs_subsets:
                if igs_id in igs_over_min_true_cov:
                    # Add guide to list
                    igs_over_min_true_cov[igs_id][0].append(guide)
                    # Add organism number to set. We can calculate the percent coverage with this number later on.
                    igs_over_min_true_cov[igs_id][1].add(aligned_targets[target_num].full_name)
                bar()
        # Now make an array of all of the putative designs for later use.
        to_optimize = np.array([RibozymeDesign(id_attr=igs_id, guides_to_use_attr=item[0], targets_attr=item[1],
                                               true_perc_cov_attr=item[2],
                                               perc_cov_attr=counts_of_igs[igs_id[:igs_len]] / aligned_targets.size)
                                for igs_id, item in igs_over_min_true_cov.items()])
        return to_optimize


def optimize_designs(to_optimize: RibozymeDesign, score_type: str, guide_len: int = 50,
                     msa_fast: bool = True, gaps_allowed: bool = True, compare_to_background: bool = False,
                     random_sample_size: int = 10, random_sample: bool = False, store_batch_results: bool = False):
    """
    Takes in a RibozymeDesign with guide list assigned and uses MUSCLE msa to optimize the guides.
    when compare_to_background is set to true, will do an MSA on background seqs and will update the anti-guide sequence

    """
    if not compare_to_background:
        random_indices = np.random.permutation(np.arange(len(to_optimize.guides_to_use)))
        sub_sample = to_optimize.guides_to_use
    else:
        random_indices = np.random.permutation(np.arange(len(to_optimize.background_guides)))
        sub_sample = to_optimize.background_guides

    if random_sample and (len(sub_sample) / random_sample_size) >= 2:
        # Extract guide batches: technically it won't split into all even chunk sizes, but it will split arrays into
        # len(sub_sample) % sample_size arrays of sample_size + 1 size, then the rest will be split into arrays of
        # the expected size
        batch_indices = np.array_split(random_indices, len(sub_sample) / random_sample_size)
        batches = np.array([(sub_sample[index] for index in indices) for indices in batch_indices])

        # Prepare to run msa
        names = np.array([f'{to_optimize.id}_{i}' for i in range(len(batches))])
        fn_muscle = np.vectorize(muscle_msa_routine, otypes=[np.ndarray, str], excluded=['muscle_exe_name', 'msa_fast',
                                                                                         'score_type', 'guide_len'])
        # Should return an array of tuples - (truncated_seq, score)
        seqs_and_scores = fn_muscle(batches, names, muscle_exe_name='muscle5', msa_fast=msa_fast,
                                    score_type=score_type, guide_len=guide_len)

        # Pick the sequence with the best score: this will be our guide sequence to optimize
        index = np.argmax(seqs_and_scores[1])
        best_guide, best_score = (seqs_and_scores[0][index], float(seqs_and_scores[1][index]))

        if store_batch_results:
            to_optimize.guide_batches = [(guide, score) for guide, score in zip(seqs_and_scores[0], seqs_and_scores[1])]

    else:
        best_guide, best_score = muscle_msa_routine(sequences_to_align=sub_sample, name_of_file=to_optimize.id,
                                                    muscle_exe_name='muscle5', msa_fast=msa_fast, score_type=score_type,
                                                    guide_len=guide_len)
        if not compare_to_background:
            to_optimize.update_after_optimizing(score_attr=best_score, guide_attr=best_guide,
                                                score_type_attr=score_type,
                                                reset_guides=True)
        else:
            to_optimize.anti_guide = best_guide[-len(to_optimize.guide):]
            to_optimize.anti_guide_score = best_score

        if store_batch_results:
            to_optimize.guide_batches = [(best_guide, best_score)]
        return to_optimize

    # Now update our initial design
    if not compare_to_background:
        # Set initial guide
        to_optimize.guide = best_guide
        to_optimize.score = best_score
        to_optimize.score_type = score_type

        # Prepare the rest of the consensus sequences to do a refinement step!
        seqs_for_refine_step = np.delete(seqs_and_scores[0], index)

        # Do a pairwise analysis: keep decreasing the ambiguity of our sequence
        for other_guide in seqs_for_refine_step:
            if best_score == 1:
                break

            to_optimize.anti_guide = other_guide
            new_seq_score, new_seq = replace_ambiguity(sequence_to_fix=to_optimize, target_background=True, n_limit=1,
                                                       update_score=False)

            if new_seq_score > best_score:
                best_guide = new_seq
                best_score = new_seq_score
                to_optimize.guide = best_guide
                to_optimize.score = best_score

        to_optimize.update_after_optimizing(score_attr=best_score, guide_attr=best_guide,
                                            score_type_attr=score_type, reset_guides=True)
    else:
        # In case the sequence is using min_len that is different than the guide length
        to_optimize.anti_guide = best_guide[-len(to_optimize.guide):]
        to_optimize.anti_guide_score = best_score
    return to_optimize


def get_naive_score(opti_seq):
    score = 1 - str(opti_seq).count('N') / len(opti_seq)
    return score


def get_quantitative_score(left_seq, alignments, chars_to_ignore: list[str] = None, count_gaps: bool = True,
                           penalize_trailing_gaps: bool = False):
    """
    :param alignments:
    :param left_seq:
    :param chars_to_ignore:
    :param count_gaps:
    :param penalize_trailing_gaps:
    :return:
    """
    # a modified version of BioPython'string_to_analyze PSSM function that counts gaps by default and also gives us
    # the quantitative score. This quantitative score is as follows:

    if chars_to_ignore is None:
        chars_to_ignore = []
    if not isinstance(chars_to_ignore, list):
        raise TypeError('chars_to_ignore should be a list.')

    gap_char = '-'
    if not count_gaps:
        chars_to_ignore.append(gap_char)

    sum_probabilities = 0
    # now start looping through all of the sequences and getting info
    seq_num = len(alignments)
    for residue_num in range(len(left_seq)):
        pos_prob = 0
        for record in alignments:
            try:
                this_residue = record[residue_num]
            # if we hit an index error we've run out of sequence and
            # should not add new residues
            except IndexError:
                this_residue = None
            if this_residue and this_residue not in chars_to_ignore:
                try:
                    # If the current residue matches the consensus residue
                    if this_residue == left_seq[residue_num]:
                        pos_prob += 1
                except KeyError:
                    raise ValueError(
                        f'Residue %s not found {this_residue}'
                    ) from None

        sum_probabilities += pos_prob / seq_num

    # Replace gaps with N since we can't order gaps when ordering oligos
    opti_seq = Bio.motifs.create(alignments, alphabet='GATCRYWSMKHBVDN-').degenerate_consensus.strip('-').replace('-',
                                                                                                                  'N')
    score = sum_probabilities / len(left_seq)

    return score, opti_seq


# noinspection DuplicatedCode
def get_weighted_score(opti_seq):
    # prepare scoring matrix a(x)
    a = {'A': 1, 'T': 1, 'G': 1, 'C': 1, 'U': 1, 'R': 0.5, 'Y': 0.5, 'M': 0.5, 'K': 0.5, 'S': 0.5, 'W': 0.5,
         'H': 1 / 3, 'B': 1 / 3, 'D': 1 / 3, 'V': 1 / 3, 'N': 0, '-': 0}
    n = len(opti_seq)

    prelim_score = 0
    for i, x in enumerate(str(opti_seq)):
        prelim_score += a[x]
    score = prelim_score / n

    return score


# noinspection DuplicatedCode
def get_directional_score(opti_seq):
    # prepare scoring matrix a(x)
    a = {'A': 1, 'T': 1, 'G': 1, 'C': 1, 'U': 1, 'R': 0.5, 'Y': 0.5, 'M': 0.5, 'K': 0.5, 'S': 0.5, 'W': 0.5,
         'H': 1 / 3, 'B': 1 / 3, 'D': 1 / 3, 'V': 1 / 3, 'N': 0, '-': 0}
    n = len(opti_seq)
    # This is broken, have not fixed yet
    # score based on extended identity and position of the bases
    weight_sum = 0
    weight_score_sum = 0

    for i, x in enumerate(str(opti_seq)):
        w = exp(n - i)  # Calculate weight based on proximity to IGS
        weight_sum += w
        weight_score_sum += w * a[x]
        i += 1  # move downstream
    score = weight_score_sum / weight_sum
    return score


def couple_designs_to_test_seqs(designs_input: str, test_seqs_input: str, flexible_igs: bool = False):
    # Check if we've pre_processed designs correctly
    if (designs_input[-7:] != '.pickle') | (test_seqs_input[-7:] != '.pickle'):
        print('Make sure to prepare designs and test sequences properly first.')
        return

    # Extract design sequence data
    with open(designs_input, 'rb') as handle:
        designs = pickle.load(handle)
    # Extract test sequences data
    with open(test_seqs_input, 'rb') as handle:
        test_seqs = pickle.load(handle)

    # Prepare data for check_guide_stats:
    # Need the design, the aligned target sequences to test against for each design
    with alive_bar(len(designs), spinner='fishes') as bar:
        for design in designs:
            matching_ref_idx_seqs = test_seqs[design.ref_idx]
            # Set correct u conservation
            design.u_conservation_test = matching_ref_idx_seqs[1]
            # Set what targets have this IGS id
            # Recall each inner dictionary is of form IGS_id: [guides, targets, perc cov, perc on target, true perc cov]
            design_id = design.igs + str(design.ref_idx)
            design.test_targets = matching_ref_idx_seqs[0][design_id][1]
            # Set igs stats
            design.perc_cov_test = matching_ref_idx_seqs[0][design_id][2]
            design.perc_on_target_test = matching_ref_idx_seqs[0][design_id][3]
            design.true_perc_cov_test = matching_ref_idx_seqs[0][design_id][4]
            # Extract the relevant guide sequences
            # (same ref_idx, and either guides with same IGS OR all guides if IGS is flexible)
            if not flexible_igs:
                igses_to_test = design.igs + str(design.ref_idx)
            else:
                igses_to_test = matching_ref_idx_seqs[0].keys()

            guides = []
            for igs_id in igses_to_test:
                guides.extend(matching_ref_idx_seqs[0][igs_id][0])
            design.guides_to_use = guides
            bar()

    # save into a separate file
    pickle_file_name = designs_input[:-7] + '_vs_' + test_seqs_input.split('.')[0].split('/')[-1] + '.coupled'
    with open(pickle_file_name, 'wb') as handle:
        pickle.dump(designs, handle)
    return pickle_file_name


def ribo_checker(coupled_designs_and_test_folder: str, number_of_workers: int, n_limit: int = 0):
    """
    This is the parallelization script I need to work on.
    """
    # Check if we have anything to test
    try:
        analysis_files = [file for file in os.listdir(coupled_designs_and_test_folder) if file[-8:] == '.coupled']
    except NotADirectoryError:
        if coupled_designs_and_test_folder[-8:] == '.coupled':
            analysis_files = [coupled_designs_and_test_folder]
        coupled_designs_and_test_folder = coupled_designs_and_test_folder.split('/')[:-1]

    if len(analysis_files) == 0:
        print('Please make sure to couple designs with the appropriate test sequences')
        return

    # Extract test sequences data
    to_test = []
    for file in analysis_files:
        file_name = file.split('.')[0]
        with open(file, 'rb') as handle:
            to_test_temp = pickle.load(handle)
        for design in to_test_temp:
            to_test.append((design, file_name))

    try:
        with open(file_name + '_checkpoint.txt', 'r') as d:
            last_design_analyzed = int(d.read())
    except FileNotFoundError:
        last_design_analyzed = 0

    # testing if this works first
    for current_design, (design_to_test, file_name) in enumerate(to_test):
        if current_design < last_design_analyzed:
            continue
        result = compare_to_test(design_to_test, n_limit=n_limit, test_dataset_name=file_name)
        result_dict = result.to_dict(all_data=False)
        with open(file_name + '_results.txt', 'a') as d:
            d.write(json.dumps(result_dict))

        with open(file_name + '_checkpoint.txt', 'w') as d:
            d.write(str(current_design))

        last_design_analyzed = current_design + 1

    return


def compare_to_test(coupled_design: RibozymeDesign, n_limit, test_dataset_name):
    def set_no_targets(ribo_design_attr: RibozymeDesign, n_limit_attr, test_dataset_name_attr):
        # No hits in the background! filter out those that have too many Ns
        # Scores are set at 0 based on the assumption that if there is no catalytic U there is no splicing
        # activity
        guide_len = len(ribo_design_attr.guide)
        # remove designs with too many Ns
        if ribo_design_attr.guide.count('N') / guide_len <= n_limit_attr:
            # Follow the same naming as the outputs of replace_ambiguity.
            ribo_design_attr.true_perc_cov_test = 0
            ribo_design_attr.perc_cov_test = 0
            ribo_design_attr.perc_on_target_test = 0
            ribo_design_attr.update_to_test(test_score_attr=0, name_of_test_dataset_attr=test_dataset_name_attr)
            return ribo_design_attr
        else:
            # # free up memory
            # del ribo_design_attr
            return None

    # See if we can skip computation all together because if no conserved u, no conserved IGS
    if not coupled_design.u_conservation_test:
        out = set_no_targets(coupled_design, n_limit, test_dataset_name)
        return out

    # Find average guide score
    guides_to_optimize = np.array([str(guide) for guide in coupled_design.guides_to_use])
    fn_pairwise = np.vectorize(pairwise_comparison, otypes=[RibozymeDesign],
                               excluded=['seq_b', 'score_type', 'only_consensus'])

    scores = fn_pairwise(guides_to_optimize, seq_b=coupled_design.guide, score_type=coupled_design.score_type,
                         only_consensus=True)
    coupled_design.update_to_test(test_score_attr=scores.mean(), name_of_test_dataset_attr=test_dataset_name)

    return coupled_design


# def test_ribo_design(design: str, test_folders: list[str], ref_seq_folder: str, igs_len: int = 5,
#                      score_type: str = 'naive', file_out: bool = False, folder_to_save: str = '', taxonomy=''):
#     print('Testing ribozyme design!')
#     start = time.perf_counter()
#     # first, read in all our sequences
#     targets_to_test_list = []
#     for file in test_folders:
#         targets_to_test_list.append(read_silva_fasta(in_file=file))
#     ref_name_and_seq = read_fasta(ref_seq_folder)[0]
#     igs = design[-igs_len:]
#     design_guide = Seq(design[:-igs_len - 1]).upper().back_transcribe()  # recall there's a G between the igs and guide
#
#     time1 = time.perf_counter()
#     print(f'Now finding catalytic sites and aligning to reference {ref_name_and_seq[0].replace("_", " ")}...')
#     # Align the targets to test to our reference sequence
#     names_guides_etc = []
#     for i, targets_to_test in enumerate(targets_to_test_list):
#         print(f'Pre-processing dataset {i} with {len(targets_to_test)} sequences to test')
#         with alive_bar(unknown='fish', spinner='fishes') as bar:
#             fn_align_to_ref = np.vectorize(align_to_ref, otypes=[TargetSeq],
#                                            excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])
#             fn_align_to_ref(targets_to_test, ref_name_and_seq=ref_name_and_seq, igs_length=igs_len,
#                             guide_length=len(design_guide), min_length=len(design_guide))
#             bar()
#         time2 = time.perf_counter()
#         round_convert_time(start=time1, end=time2, round_to=4,
#                            task_timed=f'finding catalytic sites and indexing target sequences to reference')
#
#         # get all guides and the name of their targets that have the correct igs: (og_idx, ref_idx, igs, guide)
#         time1 = time.perf_counter()
#         print('Finding matching IGSes and corresponding guides...')
#         names_temp = defaultdict(set)
#         guides_temp = defaultdict(list)
#         all_names_temp = set()
#         with alive_bar(len(targets_to_test), spinner='fishes') as bar:
#             for target in targets_to_test:
#                 for putative_guide in target.cat_sites:
#                     if putative_guide[2] == igs:
#                         names_temp[putative_guide[1]].add(target.full_name)
#                         guides_temp[putative_guide[1]].append(putative_guide[3])
#                         all_names_temp.add(target.full_name)
#                 bar()
#         names_guides_etc.append((names_temp, guides_temp, all_names_temp))
#
#     num_of_seqs_with_igs_list = [len(names) for _, _, names in names_guides_etc]
#     time2 = time.perf_counter()
#     round_convert_time(start=time1, end=time2, round_to=4,
#                        task_timed=f'matching IGSes and guides.')
#
#     # Calculate coverage for each and make a ribozyme design for each
#     time1 = time.perf_counter()
#     print('Calculating coverage for putative locations targeted by guide...')
#
#     locations_and_scores_list = []
#     for (names, guides, all_names), num_of_seqs_with_igs, targets_to_test, test_folder in (
#             zip(names_guides_etc, num_of_seqs_with_igs_list, targets_to_test_list, test_folders)):
#         locations_and_scores = np.array([], dtype=RibozymeDesign)
#         with alive_bar(len(guides.items()), spinner='fishes') as bar:
#             for index, guide_list in guides.items():
#                 perc_cov = num_of_seqs_with_igs / targets_to_test.size
#                 true_perc_cov = len(guide_list) / targets_to_test.size
#                 per_on_target = perc_cov / perc_cov
#                 score = return_score_from_type(sequence_to_test=design_guide, score_type=score_type)
#                 locations_and_scores = np.append(locations_and_scores,
#                                                  RibozymeDesign(id_attr=f'{igs}{index}', guide_attr=Seq(design_guide),
#                                                                 score_attr=score,
#                                                                 score_type_attr=score_type,
#                                                                 perc_cov_test_attr=perc_cov,
#                                                                 perc_on_target_test_attr=per_on_target,
#                                                                 true_perc_cov_test_attr=true_perc_cov))
#                 bar()
#         locations_and_scores = ribo_checker(locations_and_scores, targets_to_test, len(ref_name_and_seq[1]),
#                                             guide_length=len(design_guide), flexible_igs=True, n_limit=1,
#                                             refine_selective=False, for_test=True, test_dataset_name=test_folder)
#         locations_and_scores_list.extend(locations_and_scores)
#
#     time2 = time.perf_counter()
#     round_convert_time(start=time1, end=time2, round_to=4,
#                        task_timed=f'scoring design against test sequences')
#
#     if file_out:
#         out_file = folder_to_save + '/tested_premade_design_'
#         write_output_file(designs=locations_and_scores_list, folder_path=out_file, taxonomy=taxonomy)
#
#     end = time.perf_counter()
#     round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
#     print('########################################################\n')
#     return locations_and_scores
#


# def check_guide_stats(ribo_design: RibozymeDesign, igs_seqs_to_locs: dict,
#                       aligned_target_sequences: np.ndarray[TargetSeq], flexible_igs: bool,
#                       conserved_u_sites: np.ndarray[np.ndarray[int]], guide_length: int, n_limit: int,
#                       refine_selective: bool = False, target_background: bool = False, for_test: bool = True,
#                       test_dataset_name: str = ''):
#     def set_no_targets(for_test_attr, ribo_design_attr, guide_len_attr, n_limit_attr, test_dataset_name_attr):
#         # No hits in the background! filter out those that have too many Ns
#         # Scores are set at 0 based on the assumption that if there is no catalytic U there is no splicing
#         # activity
#
#         # remove designs with too many Ns
#         if ribo_design_attr.guide.count('N') / guide_len_attr <= n_limit_attr:
#             # Follow the same naming as the outputs of replace_ambiguity.
#             if not for_test_attr:
#                 ribo_design_attr.true_perc_cov_background = 0
#                 ribo_design_attr.perc_cov_background = 0
#                 ribo_design_attr.perc_on_target_background = 0
#                 ribo_design_attr.update_to_background(background_score_attr=0, new_guide=ribo_design_attr.guide,
#                                                       new_score=ribo_design_attr.score, reset_guides=True)
#             else:
#                 ribo_design_attr.true_perc_cov_test = 0
#                 ribo_design_attr.perc_cov_test = 0
#                 ribo_design_attr.perc_on_target_test = 0
#                 ribo_design_attr.update_to_test(test_score_attr=0, name_of_test_dataset_attr=test_dataset_name_attr)
#             return ribo_design_attr
#         else:
#             # # free up memory
#             # del ribo_design_attr
#             return None
#
#     # Now find how many conserved IGSes there are
#     # ref_idx is in base 1 indexing, indexes for arrays are base 0
#
#     # See if we can skip computation all together because if no conserved u, no conserved IGS
#     if for_test:
#         if not ribo_design.u_conservation_test:
#             out = set_no_targets(for_test, ribo_design, guide_length, n_limit, test_dataset_name)
#             return out
#     else:
#         if not ribo_design.u_conservation_background:
#             out = set_no_targets(for_test, ribo_design, guide_length, n_limit, test_dataset_name)
#             return out
#     try:
#         matching_igses = np.array(igs_seqs_to_locs[ribo_design.igs])  # base 0 indexing
#         orgs_with_igs = np.unique(matching_igses[:, 0])
#         # Convert to base 0 indexing
#         orgs_with_igs_on_target = matching_igses[np.argwhere(matching_igses[:, 1] == ribo_design.ref_idx - 1)][:, 0, 0]
#         on_target_count = orgs_with_igs_on_target.size
#         true_perc_cov = on_target_count / aligned_target_sequences.size
#         perc_cov = orgs_with_igs.size / aligned_target_sequences.size
#         perc_on_target = true_perc_cov / perc_cov
#         # Fill this with empty TargetSeq objects to later make finding taxonomy easier
#         background_names = set(back_seq.full_name for back_seq in aligned_target_sequences[matching_igses[:, 0]])
#     except KeyError:
#         # If there are no IGSes on target, and we don't want to consider other matching catalytic sites, set no targets
#         if not flexible_igs:
#             out = set_no_targets(for_test, ribo_design, guide_length, n_limit, test_dataset_name)
#             return out
#         true_perc_cov = 0
#         perc_cov = 0
#         perc_on_target = 0
#         background_names = None
#
#     # If the igs is not flexible, only check those that have a perfect igs upstream of the cat_site
#     if not flexible_igs:
#         # Fill this with empty TargetSeq objects to later make finding taxonomy more easy
#         # comparing two base 1 indexes, no conversion needed, but this will hold base 0
#         all_indexes_of_target_data = [np.argwhere(
#             np.array([cat_site[1] for cat_site
#                       in aligned_target_sequences[target_number].cat_sites]) == ribo_design.ref_idx)[0, 0]
#                                       for target_number in orgs_with_igs_on_target]
#
#         guides_to_optimize = [str(aligned_target_sequences[target].cat_sites[index_of_target_data][3]) for
#                               target, index_of_target_data in
#                               zip(orgs_with_igs_on_target, all_indexes_of_target_data)]
#     # Otherwise, extract all info in conserved cat_sites regardless of igs. As long at the catalytic U is on
#     # target, we'll extract those data
#     else:
#         # Convert to base 0 indexing
#         orgs_with_u_on_target = conserved_u_sites[
#                                     np.where(conserved_u_sites[:, 1] == ribo_design.ref_idx - 1)][:, 0]
#         if not orgs_with_u_on_target.any():
#             out = set_no_targets(for_test, ribo_design, guide_length, n_limit, test_dataset_name)
#             return out
#         # if we want to keep only u targets
#         # background_names = set(back_seq.full_name for back_seq in aligned_target_sequences[orgs_with_u_on_target])
#         # comparing two base 1 indexes, no conversion needed, but this will hold base 0
#         # I'm not sure why the following line won't work if I don't convert the list in np.argwhere into an array
#         try:
#             all_indexes_of_target_data = [np.argwhere(
#                 np.array([cat_site[1] for cat_site
#                           in aligned_target_sequences[target_number].cat_sites]) == ribo_design.ref_idx)[0, 0]
#                                           for target_number in orgs_with_u_on_target]
#         except:
#             print('uh oh')
#         guides_to_optimize = np.array([str(aligned_target_sequences[target].cat_sites[index_of_target_data][3]) for
#                                        target, index_of_target_data in
#                                        zip(orgs_with_u_on_target, all_indexes_of_target_data)])
#
#     # If necessary, optimize further:
#     if ribo_design.score < 1 and refine_selective and not for_test:
#         best_score = ribo_design.score
#
#         # Do a pairwise analysis: keep decreasing the ambiguity of our sequence
#         for other_guide in guides_to_optimize:
#             if best_score == 1:
#                 break
#
#             ribo_design.anti_guide = other_guide
#             new_seq_score, new_seq = replace_ambiguity(sequence_to_fix=ribo_design, update_score=False,
#                                                        target_background=target_background, n_limit=1)
#
#             if new_seq_score > best_score:
#                 best_guide = new_seq
#                 best_score = new_seq_score
#                 ribo_design.guide = best_guide
#                 ribo_design.score = best_score
#     # Find average guide score
#     fn_pairwise = np.vectorize(pairwise_comparison, otypes=[RibozymeDesign],
#                                excluded=['seq_b', 'score_type', 'only_consensus'])
#
#     scores = fn_pairwise(guides_to_optimize, seq_b=ribo_design.guide, score_type=ribo_design.score_type,
#                          only_consensus=True)
#     if not for_test:
#         ribo_design.true_perc_cov_background = true_perc_cov
#         ribo_design.perc_cov_background = perc_cov
#         ribo_design.perc_on_target_background = perc_on_target
#         ribo_design.background_targets = background_names
#         ribo_design.update_to_background(background_score_attr=scores.mean(), new_guide=ribo_design.guide,
#                                          new_score=ribo_design.score, reset_guides=True)
#     else:
#         ribo_design.true_perc_cov_test = true_perc_cov
#         ribo_design.perc_cov_test = perc_cov
#         ribo_design.perc_on_target_test = perc_on_target
#         ribo_design.test_targets = background_names
#         ribo_design.update_to_test(test_score_attr=scores.mean(), name_of_test_dataset_attr=test_dataset_name)
#
#     return ribo_design


def replace_ambiguity(sequence_to_fix: RibozymeDesign, target_background: bool = False,
                      n_limit: float = 1, update_score: bool = True):
    # Here is a dictionary with all the antinucleotides for each IUPAC ambiguity code
    anti_iupac_nucleotides = {'A': 'B', 'C': 'D', 'G': 'H', 'T': 'V', 'M': 'K', 'R': 'Y', 'W': 'S', 'S': 'W', 'Y': 'R',
                              'K': 'M', 'V': 'T', 'H': 'G', 'D': 'C', 'B': 'A', 'N': 'N', '-': 'N'}
    iupac_ambiguity_nucleotides = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'M': 'AC', 'R': 'AG', 'W': 'AT', 'S': 'CG',
                                   'Y': 'CT', 'K': 'GT', 'V': 'ACG', 'H': 'ACT', 'D': 'AGT', 'B': 'CGT', 'N': '',
                                   '-': ''}
    inv_iupac_ambiguity_nucleotides = {v: k for k, v in iupac_ambiguity_nucleotides.items()}

    new_seq = ''
    n_number = 0

    if len(sequence_to_fix.guide) > len(sequence_to_fix.anti_guide):
        # Sometimes when using a minlen that is different than guide_len things get a little funky in the anti_guide
        # at the ends of sequences. This makes it so we can still score these between each other.
        pad = '-' * (len(sequence_to_fix.guide) - len(sequence_to_fix.anti_guide))
        sequence_to_fix.anti_guide += pad

    for i, base in enumerate(sequence_to_fix.guide):
        if base not in ['A', 'T', 'C', 'G']:
            # Get the base at that same position in the bad sequence
            background_base = sequence_to_fix.anti_guide[i]

            # Easiest case: if there is total ambiguity, just get the opposite base from the bad sequence
            if base == 'N':
                if not target_background:
                    # Replace that base here with the antinucleotide
                    new_base = anti_iupac_nucleotides[background_base]
                else:
                    new_base = background_base

            # Harder case: try to remove some ambiguity by narrowing down the potential opposite bases
            else:
                background_extended_bases = iupac_ambiguity_nucleotides[background_base]
                good_extended_bases = iupac_ambiguity_nucleotides[base]
                if not target_background:
                    for remove_this_from_set in background_extended_bases:
                        # the set that only includes the bases in the good base
                        good_extended_bases = good_extended_bases.replace(remove_this_from_set, '')
                else:
                    # the set that only includes bases in both background and good bases
                    good_extended_bases = str(set(good_extended_bases) & set(background_extended_bases)).strip(
                        '{\'').strip('\'}')
                # Now check that we are actually reducing ambiguity. We don't want to increase ambiguity!
                if good_extended_bases and len(good_extended_bases) < len(iupac_ambiguity_nucleotides[base]):
                    new_base = inv_iupac_ambiguity_nucleotides[good_extended_bases]
                else:
                    new_base = base
        else:
            new_base = base
        if new_base == 'N':
            n_number += 1
            if n_number / len(sequence_to_fix.guide) > n_limit:
                return

        new_seq += new_base

    if update_score:
        # Now do a pairwise sequence with the target sequence to see the delta score:
        new_seq_score, pairwise_score = pairwise_comparison(seq_a=new_seq, seq_b=sequence_to_fix.anti_guide,
                                                            score_type=sequence_to_fix.score_type)
        sequence_to_fix.update_to_background(background_score_attr=pairwise_score, new_guide=Seq(new_seq),
                                             new_score=new_seq_score, reset_guides=True)
        return
    else:
        new_seq_score = return_score_from_type(sequence_to_test=new_seq, score_type=sequence_to_fix.score_type)
        return new_seq_score, new_seq


def pairwise_comparison(seq_a: Seq | str, seq_b: Seq | str, score_type: str = 'naive', only_consensus: bool = False):
    if len(seq_a) != len(seq_b):
        pad = '-' * abs(len(seq_a) - len(seq_b))
        if len(seq_a) > len(seq_b):
            seq_b += pad
        else:
            seq_a += pad

    pairwise_comparison_consensus = Bio.motifs.create([seq_a, seq_b],
                                                      alphabet='GATCRYWSMKHBVDN-').degenerate_consensus
    pairwise_score = return_score_from_type(sequence_to_test=pairwise_comparison_consensus, score_type=score_type)

    if only_consensus:
        return pairwise_score
    else:
        seq_a_score = return_score_from_type(sequence_to_test=seq_a, score_type=score_type)
        return seq_a_score, pairwise_score


def write_output_file(designs: np.ndarray[RibozymeDesign] | list[RibozymeDesign], folder_path: str, taxonomy='',
                      all_data: bool = False):
    designs_df = pd.DataFrame.from_records([item.to_dict(taxonomy=taxonomy, all_data=all_data) for item in designs])
    if os.path.exists(f'{folder_path}.csv'):
        os.remove(f'{folder_path}.csv')
    designs_df.to_csv(f'{folder_path}.csv', index=False)
    return


def make_graphs(control_designs: np.ndarray[RibozymeDesign] | list | str,
                universal_designs: list[np.ndarray[RibozymeDesign]] | list[str],
                selective_designs: list[np.ndarray[RibozymeDesign]] | list[str],
                ref_seq_designs: np.ndarray[RibozymeDesign] | list | str,
                var_regs: list[tuple[int, int]], save_fig: bool = False, file_loc: str = None, file_type: str = 'png',
                taxonomy: str = '', data_file: str = '', save_file_loc: str = '', test_folders: [str] = None):
    if not test_folders:
        print('Please provide test data')
        return
    # Make data into a pandas dataframe
    if not data_file:
        if type(control_designs) is str:
            all_data_df = pd.read_csv(control_designs)
            all_data_df.index = ['control'] * len(all_data_df.values)
        else:
            all_data_df = pd.DataFrame.from_records([item.to_dict(taxonomy='') for item in control_designs],
                                                    index='control' * len(control_designs))

        if type(ref_seq_designs) is str:
            ref_seq_designs_df = pd.read_csv(ref_seq_designs)
            ref_seq_designs_df.index = ['reference_seq'] * len(ref_seq_designs_df.values)
        else:
            ref_seq_designs_df = pd.DataFrame.from_records([item.to_dict(taxonomy='') for item in ref_seq_designs],
                                                           index='reference_seq' * len(ref_seq_designs))
        all_data_df = pd.concat([all_data_df, ref_seq_designs_df])

        for i, sub_list in enumerate(universal_designs):
            if type(sub_list) is str:
                sub_list_df = pd.read_csv(sub_list)
                sub_list_df.index = [f'universal_{i}'] * len(sub_list_df)
            else:
                sub_list_df = pd.DataFrame.from_records([item.to_dict(taxonomy='') for item in sub_list],
                                                        index=f'universal_{i}' * len(sub_list))
            all_data_df = pd.concat([all_data_df, sub_list_df])

        for i, sub_list in enumerate(selective_designs):
            if type(sub_list) is str:
                sub_list_df = pd.read_csv(sub_list)
                sub_list_df.index = [f'selective_{i}'] * len(sub_list_df)
            else:
                sub_list_df = pd.DataFrame.from_records([item.to_dict(taxonomy='') for item in sub_list],
                                                        index=f'selective_{i}' * len(sub_list))
            all_data_df = pd.concat([all_data_df, sub_list_df])
        #
        #
        # all_data_array = np.hstack([control_designs] + universal_designs + selective_designs + [ref_seq_designs])
        # index_names = ['control'] * len(control_designs) + \
        #               [f'universal_{i}' for i, data in enumerate(universal_designs) for _ in data] + \
        #               [f'selective_{i}' for i, data in enumerate(selective_designs) for _ in data] + \
        #               ['reference_seq'] * len(ref_seq_designs)
        #
        # all_data_df = pd.DataFrame.from_records([item.to_dict(taxonomy='') for item in all_data_array],
        #                                         index=index_names)
        if os.path.exists(file_loc):
            os.remove(file_loc)
        all_data_df.to_csv(file_loc)
    else:
        all_data_df = pd.read_csv(data_file)

    # Extract top scores for each category - this way it will count for us!
    all_data_df.rename(columns={'Unnamed: 0': 'tested_targets'}, inplace=True)
    all_data_df.index.name = 'Index'
    labels = list(set(all_data_df['tested_targets']))
    labels = [label for label in labels if type(label) is str]
    universal_labels = [name for name in labels if name[0] == 'u']
    selective_labels = [name for name in labels if name[0] == 's']
    labels.sort()
    control_designs_df = all_data_df.loc[all_data_df['tested_targets'] == 'control']
    control_designs_df_1376 = control_designs_df[control_designs_df['reference_idx'] == 1376]
    universal_designs_df = all_data_df.loc[all_data_df['tested_targets'].isin(universal_labels)]
    num_universal = len(universal_labels)
    selective_designs_df = all_data_df.loc[all_data_df['tested_targets'].isin(selective_labels)]
    num_selective = len(selective_labels)
    ref_seq_designs_df = all_data_df.loc[all_data_df['tested_targets'] == 'reference_seq']

    target_names = ['reference', 'bacterial', 'archaeal', 'eukaryotic', 'all kingdoms']

    # get top scores in each: 1 control, one of each group of selective and universal
    top_score_control = control_designs_df_1376[control_designs_df_1376[
        'name_of_test_dataset'].str.contains('Bac')].nlargest(1, 'composite_test_score')
    file_name_checks = ['Bac', 'Euk', 'Arc', 'All']
    for name in file_name_checks[1:]:
        top_score_temp = control_designs_df_1376[control_designs_df_1376[
            'name_of_test_dataset'].str.contains(name)].nlargest(1, 'composite_test_score')
        top_score_control = pd.concat([top_score_control, top_score_temp])

    top_scores_universal = universal_designs_df.loc[universal_designs_df['tested_targets'] ==
                                                    universal_labels[0]].nlargest(1, 'composite_test_score')
    for label in universal_labels[1:]:
        score_temp = universal_designs_df[universal_designs_df['tested_targets'] ==
                                          label].nlargest(1, 'composite_test_score')
        top_scores_universal = pd.concat([top_scores_universal, score_temp])

    top_scores_selective = selective_designs_df.loc[selective_designs_df['tested_targets'] ==
                                                    selective_labels[0]].nlargest(1, 'composite_test_score')
    for label in selective_labels[1:]:
        score_temp = selective_designs_df.loc[selective_designs_df['tested_targets'] ==
                                              label].nlargest(1, 'composite_test_score')
        top_scores_selective = pd.concat([top_scores_selective, score_temp])

    # top_test_scores = pd.concat([top_score_control, top_scores_universal, top_scores_selective])

    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.8, 16 * 0.8)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # colors = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']  # viridis
    colors = ['#0D365E', '#3F6D54', '#9B892D', '#F9A281', '#FACDFB']  # batlow
    to_analyze = ['reference_seq', 'universal_0', 'universal_1', 'universal_2', 'universal_3']

    max_vals = all_data_df.max(numeric_only=True)

    # Fig 1: MG1655
    reference_subset = ref_seq_designs_df[ref_seq_designs_df['name_of_test_dataset'].str.contains('Bac')]
    bacterial_subset = universal_designs_df[universal_designs_df['name_of_test_dataset'].str.contains('Bac')]
    bacterial_subset = bacterial_subset.loc[bacterial_subset['tested_targets'] == 'universal_0']
    top_score_control_subset = top_score_control[top_score_control['name_of_test_dataset'].str.contains('Bacteria')]
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
        sns.scatterplot(x=xvar, y=yvar, hue=hue_var, data=graph, ax=joint_ax[0], alpha=0.5, legend=False)
        jointplot_fig.axes[0].scatter(x=top_score_control_subset[xvar], y=top_score_control_subset[yvar],
                                      alpha=1, c='#fde725', label='Control')
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
        jointplot_fig.suptitle(f'Designs from {name} target sequences against bacterial test datasets')

        if save_fig:
            plt.savefig(fname=f'{save_file_loc}/figure_1_{name}.{file_type}', format=file_type)
        plt.show()

    # Fig 2: : Assessing universal design quality. 2a) IGS true percent coverage vs. guide score of bacterial designs
    # evaluated against datasets of different kingdoms. 2b) 16s rRNA location of all designs along the reference
    # E. coli MG1655 16s rRNA sequence.
    for i, target_name in enumerate(target_names):
        universal_subset = all_data_df.loc[all_data_df['tested_targets'] == to_analyze[i]]

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
            txt = label.get_text().split('/')[-1][0:3]
            new_label = [name for name in legend_names if txt.casefold() in name.casefold()]
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
    return


def give_taxonomy(silva_name: str | TargetSeq, level: str):
    if type(silva_name) is TargetSeq:
        return silva_name.give_taxonomy(level)
    level_names = {'Domain': 0, 'Phylum': 1, 'Class': 2, 'Order': 3,
                   'Family': 4, 'Genus': 5, 'Species': 6, 'Taxon': 7}
    if level not in level_names:
        print('Taxonomic level not found')
    # gives us the name at a certain taxonomic level
    putative_names = silva_name.split(',')
    out_names = []
    for name in putative_names:
        all_levels = name.split(';')
        try:
            if level == 'Domain':
                here_is_taxonomy = all_levels[0].split(' ')[-1]
            else:
                here_is_taxonomy = all_levels[level_names[level]]
            out_names.append(here_is_taxonomy)
        except IndexError:
            out_names.append('unknown')
            continue
    if out_names:
        counts_and_names = np.asarray(np.unique(out_names, return_counts=True)).T
        return counts_and_names.tolist()
    else:
        return None


def plot_variable_regions(ax, var_regs):
    for V in var_regs:
        ax.axvspan(V[0], V[1], facecolor='g', alpha=0.2)
    return


def muscle_msa_routine(sequences_to_align, name_of_file: str, muscle_exe_name: str, msa_fast: bool,
                       score_type: str, guide_len: int):
    # First create a dummy file with all the guides we want to optimize
    with open(f'to_align_{name_of_file}.fasta', 'w') as f:
        for i, line in enumerate(sequences_to_align):
            f.write('>seq' + str(i) + '\n' + str(line) + '\n')
    if msa_fast:
        subprocess.check_output([muscle_exe_name, '-super5', f'to_align_{name_of_file}.fasta', '-output',
                                 f'aln_{name_of_file}.afa'], stderr=subprocess.DEVNULL)
    else:
        subprocess.check_output([muscle_exe_name, '-align', f'to_align_{name_of_file}.fasta', '-output',
                                 f'aln_{name_of_file}.afa'], stderr=subprocess.DEVNULL)

    msa = AlignIO.read(open(f'aln_{name_of_file}.afa'), 'fasta')

    # delete files:
    if os.path.exists(f'aln_{name_of_file}.afa'):
        os.remove(f'aln_{name_of_file}.afa')
    if os.path.exists(f'to_align_{name_of_file}.fasta'):
        os.remove(f'to_align_{name_of_file}.fasta')

    # Now retrieve the alignment
    summary_align = AlignInfo.SummaryInfo(msa)
    alignments = [record.seq for record in summary_align.alignment]

    opti_seq = Bio.motifs.create(alignments, alphabet='GATCRYWSMKHBVDN-').degenerate_consensus.strip('-')
    # Make sure our sequence is our desired length
    truncated_seq = opti_seq[-guide_len:]

    score = return_score_from_type(sequence_to_test=truncated_seq, score_type=score_type)

    return truncated_seq, score


def return_score_from_type(sequence_to_test: Seq | str, score_type: str):
    if score_type == 'naive':
        # Naively tell me what the percent identity is: 1- ambiguity codons/ length
        score = get_naive_score(sequence_to_test)

    elif score_type == 'weighted':
        score = get_weighted_score(sequence_to_test)

    elif score_type == 'directional':
        score = get_directional_score(sequence_to_test)

    else:
        print('Score type not recognized. Defaulting to naive scoring')
        score = get_naive_score(sequence_to_test)
    return score


def compare_batches(input_designs: np.ndarray[RibozymeDesign]):
    # Extract batch data
    data_vals = {design.id: [] for design in input_designs}
    data_seqs = {design.id: [] for design in input_designs}
    data_guides = {design.id: [] for design in input_designs}
    for design in input_designs:
        key = design.id
        data_guides[key].append(design.guide)
        for guide_seq, val in design.guide_batches:
            data_vals[key].append(val)
            data_seqs[key].append(guide_seq)

    # Check how many sequences are identical
    set_seqs = {}
    for id, seqs in data_seqs.items():
        test = set(seqs)
        if len(test) > 1:
            set_seqs[id] = test
    # Check how many scores are identical
    set_scores = {}
    for id, scores in data_vals.items():
        test = set(scores)
        if len(test) > 1:
            set_scores[id] = test
    # Check if the obtained guide sequence is identical
    set_guides = {}
    for id, guides in data_guides.items():
        test = set(guides)
        if len(test) > 1:
            set_guides[id] = test

    # plot score distributions per id

    # plot sequence map of generated sequences per id

    # plot sequence map of final guide per id

    return
