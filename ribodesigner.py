import glob
import json
import multiprocessing as mp
import os
import pickle
import subprocess
import time
from collections import defaultdict
from math import log
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqUtils
from Bio.Align import AlignInfo, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from alive_progress import alive_bar
from icecream import ic
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

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

    def __init__(self, id_attr: str = '', guides_to_use_attr: list[Seq] = None, targets_attr: set = None,
                 guide_attr: Seq = '', score_attr: float = None, score_type_attr: str = '', perc_cov_attr: float = None,
                 perc_on_target_attr: float = None, true_perc_cov_attr: float = None,
                 background_tm_nn_attr: float = None,
                 background_score_attr: float = None, perc_cov_background_attr: float = None,
                 perc_on_target_background_attr: float = None, true_perc_cov_background_attr: float = None,
                 background_guides_attr: list[Seq] = None, anti_guide_attr: Seq = '', anti_guide_score_attr: int = None,
                 background_targets_attr: set = None, igs_attr: str = '', ref_idx_attr: int = None,
                 u_consv_background_attr: float = None, tested_design_attr: bool = False,
                 perc_cov_test_attr: float = None,
                 perc_on_target_test_attr: float = None, true_perc_cov_test_attr: float = None,
                 dict_initialize: dict = None):
        if dict_initialize:
            # First, do most data
            self.id = dict_initialize['id']
            self.igs = dict_initialize['igs']
            self.ref_idx = dict_initialize['reference_idx']
            self.optimized_to_targets = dict_initialize['optimized_to_targets']
            self.optimized_to_background = dict_initialize['optimized_to_background']
            self.tested = dict_initialize['tested']
            self.tested_design = dict_initialize['tested_design']
            self.guide = dict_initialize['guide']
            self.number_of_targets = dict_initialize['num_of_targets']
            self.score_type = dict_initialize['score_type']
            self.score = dict_initialize['score']
            self.perc_cov = dict_initialize['U_IGS_coverage']
            self.perc_on_target = dict_initialize['U_IGS_on target']
            self.true_perc_cov = dict_initialize['true_U_IGS_cov']
            self.composite_score = dict_initialize['composite_score']
            self.number_of_targets_background = dict_initialize['num_of_targets_background']
            self.u_conservation_background = dict_initialize['u_conservation_background']
            self.background_score = dict_initialize['background_score']
            self.perc_cov_background = dict_initialize['U_IGS_coverage_background']
            self.perc_on_target_background = dict_initialize['U_IGS_on target_background']
            self.true_perc_cov_background = dict_initialize['true_U_IGS_cov_background']
            self.composite_background_score = dict_initialize['composite_background_score']
            self.delta_igs_vs_background = dict_initialize['delta_igs_vs_background']
            self.delta_guide_vs_background = dict_initialize['delta_guide_vs_background']
            self.delta_vs_background = dict_initialize['delta_vs_background']
            self.name_of_test_dataset = dict_initialize['name_of_test_dataset']
            self.number_of_targets_test = dict_initialize['num_of_targets_test']
            self.u_conservation_test = dict_initialize['u_conservation_test']
            self.test_score = dict_initialize['test_score']
            self.test_tm_nn = dict_initialize['tm_nn_vs_test']
            self.perc_cov_test = dict_initialize['U_IGS_coverage_test']
            self.perc_on_target_test = dict_initialize['U_IGS_on target_test']
            self.true_perc_cov_test = dict_initialize['true_U_IGS_cov_test']
            self.composite_test_score = dict_initialize['composite_test_score']
            self.delta_igs_vs_test = dict_initialize['delta_igs_vs_test']
            self.delta_guide_vs_test = dict_initialize['delta_guide_vs_test']
            self.delta_vs_test = dict_initialize['delta_vs_test']
            # Then, targets will be a bit messy. Should fix this later, right now I just need it to work
            self.targets = zip(dict_initialize['target_Domain'], dict_initialize['target_Phylum'],
                               dict_initialize['target_Class'], dict_initialize['target_Order'],
                               dict_initialize['target_Family'], dict_initialize['target_Genus'],
                               dict_initialize['target_Species'], dict_initialize['target_Taxon'])
            self.test_targets = zip(dict_initialize['target_Domain_test'], dict_initialize['target_Phylum_test'],
                                    dict_initialize['target_Class_test'], dict_initialize['target_Order_test'],
                                    dict_initialize['target_Family_test'], dict_initialize['target_Genus_test'],
                                    dict_initialize['target_Species_test'], dict_initialize['target_Taxon_test'])
            self.background_targets = zip(dict_initialize['target_Domain_background'],
                                          dict_initialize['target_Phylum_background'],
                                          dict_initialize['target_Class_background'],
                                          dict_initialize['target_Order_background'],
                                          dict_initialize['target_Family_background'],
                                          dict_initialize['target_Genus_background'],
                                          dict_initialize['target_Species_background'],
                                          dict_initialize['target_Taxon_background'])
        else:
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
            self.background_tm_nn = background_tm_nn_attr
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
            self.test_tm_nn = None
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
        if type(taxonomy) is list:
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
            'U_IGS_coverage': self.perc_cov,
            'U_IGS_on target': self.perc_on_target,
            'true_U_IGS_cov': self.true_perc_cov,
            'composite_score': self.composite_score,
            'num_of_targets_background': self.number_of_targets_background,
            'u_conservation_background': self.u_conservation_background,
            'background_score': self.background_score,
            # 'tm_nn_vs_background': self.background_tm_nn,
            'U_IGS_coverage_background': self.perc_cov_background,
            'U_IGS_on target_background': self.perc_on_target_background,
            'true_U_IGS_cov_background': self.true_perc_cov_background,
            'composite_background_score': self.composite_background_score,
            'delta_igs_vs_background': self.delta_igs_vs_background,
            'delta_guide_vs_background': self.delta_guide_vs_background,
            'delta_vs_background': self.delta_vs_background,
            'name_of_test_dataset': self.name_of_test_dataset,
            'num_of_targets_test': self.number_of_targets_test,
            'u_conservation_test': self.u_conservation_test,
            'test_score': self.test_score,
            'tm_nn_vs_test': self.test_tm_nn,
            'U_IGS_coverage_test': self.perc_cov_test,
            'U_IGS_on target_test': self.perc_on_target_test,
            'true_U_IGS_cov_test': self.true_perc_cov_test,
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
            return
        if not perc_cov_attr and true_perc_cov_attr and perc_on_target_attr:
            self.perc_cov = true_perc_cov_attr / perc_on_target_attr

        if not perc_on_target_attr and true_perc_cov_attr and perc_cov_attr:
            self.perc_on_target = true_perc_cov_attr / perc_cov_attr

        if not true_perc_cov_attr and perc_on_target_attr and perc_cov_attr:
            self.true_perc_cov = perc_on_target_attr * perc_cov_attr
        return

    def calc_background_percent_coverages(self, perc_cov_background_attr: float = None,
                                          true_perc_cov_background_attr: float = None,
                                          perc_on_target_background_attr: float = None):
        self.perc_cov_background = perc_cov_background_attr
        self.true_perc_cov_background = true_perc_cov_background_attr
        self.perc_on_target_background = perc_on_target_background_attr

        if not perc_on_target_background_attr and true_perc_cov_background_attr and perc_on_target_background_attr:
            self.perc_cov_background = true_perc_cov_background_attr / perc_on_target_background_attr

        if not perc_on_target_background_attr and true_perc_cov_background_attr and perc_cov_background_attr:
            self.perc_on_target_background = true_perc_cov_background_attr / perc_cov_background_attr

        if not true_perc_cov_background_attr and perc_on_target_background_attr and perc_cov_background_attr:
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
                             reset_guides: bool = False, tm_nn_attr: float = None):
        self.optimized_to_background = True
        self.background_score = background_score_attr
        self.score = new_score
        self.background_tm_nn = tm_nn_attr
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

    def update_to_test(self, test_score_attr: float, name_of_test_dataset_attr: str, test_tm_nn_attr: float = None):
        self.tested = True
        self.name_of_test_dataset = name_of_test_dataset_attr
        self.test_score = test_score_attr
        self.test_tm_nn = test_tm_nn_attr
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
            # Score on targets,U_IGS cov background,U_IGS on target background,true U_IGS cov background,Composite score background,
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
                 guide_length: int = 50, min_length: int = 35, ref_sequence_file=None, min_true_cov: float = 0.7,
                 fileout: bool = False, folder_to_save: str = '',
                 score_type: str = 'weighted', msa_fast: bool = False, gaps_allowed: bool = True,
                 percent_of_target_seqs_used: float = 1.0, get_consensus_batches: bool = True,
                 random_guide_sample_size: int = 10, store_batch_results: bool = False):
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
    #                      random_sample_size=random_guide_sample_size, random_sample=get_consensus_batches)
    #     ic()
    # time1 = time.perf_counter()
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        # vectorize function. For the multiprocessing module I made it return the updated RibozymeDesign

        fn_msa = np.vectorize(optimize_designs, otypes=[RibozymeDesign],
                              excluded=['score_type', 'guide_len', 'msa_fast', 'gaps_allowed', 'compare_to_background',
                                        'random_sample_size', 'random_sample', 'store_batch_results'])
        # prepare data for starmap
        in_data = [(sub_array, score_type, guide_length, msa_fast, gaps_allowed, False,
                    random_guide_sample_size, get_consensus_batches, store_batch_results)
                   for sub_array in np.array_split(to_optimize, process_num)]

        # multithread
        with mp.Pool(processes=process_num) as pool:
            out_data = pool.starmap(fn_msa, in_data)
        optimized_seqs = np.concatenate([sub_array for sub_array in out_data])
        bar()

    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='generating optimized designs')

    pickle_file_name = target_sequences_folder.split('.')[0].split('/')[-1].split('\\')[-1] + '_universal'
    print('Now pickling output file...')
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        with open(os.path.normpath(f'{folder_to_save}/designs_{pickle_file_name}.pickle'), 'wb') as handle:
            pickle.dump(optimized_seqs, handle)
            bar()

    if fileout:
        if not os.path.exists(folder_to_save):
            os.mkdir(folder_to_save)
        out_file = os.path.normpath(folder_to_save + '/designs_' + pickle_file_name)
        write_output_file(designs=optimized_seqs, folder_path=out_file, all_data=store_batch_results)

    end = time.perf_counter()
    round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    print('########################################################\n')
    return os.path.normpath(f'{folder_to_save}/designs_{pickle_file_name}.pickle')


def prepare_test_seqs(test_folder, ref_sequence_file, guide_length, igs_length, min_length, folder_to_save,
                      graph_results, var_regs, graph_file_type, get_consensus_batches, batch_num, score_type, msa_fast,
                      remove_x_dupes_in_graph, lim):
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
    test_seqs_dict, igs_ids_counts, counts_of_ref_idx = filter_igs_candidates(aligned_targets=test_names_and_seqs,
                                                                              min_true_cov=0, igs_len=igs_length,
                                                                              return_test_seqs=True)
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='finding putative ribozyme sites')

    time1 = time.perf_counter()
    guide_scores = {}

    if get_consensus_batches:
        print('Getting consensus sequences for test locations...')
        fn_muscle = np.vectorize(muscle_msa_routine, otypes=[str, float], excluded=['muscle_exe_name', 'msa_fast',
                                                                                         'score_type', 'guide_len'])
        with alive_bar(len(test_seqs_dict), spinner='fishes') as bar:
            for idx, loc_data in test_seqs_dict.items():
                for idx_id, idx_id_data in loc_data[0].items():
                    # divide number of guide sequences into batches of around batch_num guides
                    guides_to_use = idx_id_data[0]
                    # If there is only one guide, just keep it and continue
                    if len(guides_to_use) == 1:
                        guides = guides_to_use
                        score = 1
                    elif (len(guides_to_use) / batch_num) >= 2:
                        random_indices = np.random.permutation(np.arange(len(guides_to_use)))

                        batch_indices = np.array_split(random_indices, len(guides_to_use) / batch_num)
                        batches = np.array([(guides_to_use[index] for index in indices) for indices in batch_indices])

                        # Prepare to run msa
                        names = np.array([f'{idx_id}_{i}' for i in range(len(batches))])
                        # Should return an array of tuples - (truncated_seq, score)
                        seqs_and_scores = fn_muscle(batches, names, muscle_exe_name='muscle5', msa_fast=msa_fast,
                                                    score_type=score_type, guide_len=guide_length)
                        guides = [*seqs_and_scores[0]]
                        # Get the average score
                        score = sum(seqs_and_scores[1])/ len(seqs_and_scores[1])

                    else:
                        best_guide, best_score = muscle_msa_routine(sequences_to_align=guides_to_use, name_of_file=idx_id,
                                                           muscle_exe_name='muscle5', msa_fast=msa_fast,
                                                           score_type=score_type, guide_len=guide_length)
                        guides = [best_guide]
                        score = best_score
                    test_seqs_dict[idx][0][idx_id][0] = guides
                    # test_seqs_dict[idx][0][idx_id].append(score)
                    guide_scores[idx_id] = score
                    if len(test_seqs_dict[idx][0][idx_id]) > 6:
                        print('uh oh')
                bar()
        time2 = time.perf_counter()
        round_convert_time(start=time1, end=time2, round_to=4, task_timed='getting consensus sequences')

    # prepare data for graphing
    if not os.path.exists(folder_to_save):
        os.mkdir(folder_to_save)

    title = test_folder.split('.')[0].split('/')[-1].split('\\')[-1]
    save_file_name = os.path.normpath(f'{folder_to_save}/test_sequences_{title}')
    if graph_results:
        # Graph IGS data
        num_of_seqs = len(test_names_and_seqs)
        u_conservation_list = []
        igs_true_perc_cov_list = []
        ref_idxes_list = []
        guide_scores_list = []
        in_var_reg = []
        for igs_id, igs_id_count in igs_ids_counts.items():
            ref_idx = int(igs_id[igs_length:])
            ref_idxes_list.append(ref_idx)
            u_conservation_list.append(counts_of_ref_idx[ref_idx] / num_of_seqs)
            igs_true_perc_cov_list.append(igs_id_count / num_of_seqs)
            guide_scores_list.append(guide_scores[igs_id])
            # for l, r in var_regs:
            #     if l <= ref_idx <= r:
            #         in_var_reg_temp = True
            #         break
            #     else:
            #         in_var_reg_temp = False
            # in_var_reg.append(in_var_reg_temp)
        make_test_seqs_graph(title, x_data=u_conservation_list, xlabel='U coverage',
                             y_data=igs_true_perc_cov_list, ylabel='IGS true coverage', loc_data=ref_idxes_list,
                             var_regs=var_regs, save_file_name=save_file_name, file_type=graph_file_type, alpha=0.3,
                             dataset_len=num_of_seqs)
        # Graph guide score data
        make_guide_score_plot(xdata=igs_true_perc_cov_list, xlabel='IGS true coverage', ydata=guide_scores_list,
                              ylabel='Average guide score', loc_data=ref_idxes_list, var_regs=var_regs, lim=lim,
                              save_file_name=save_file_name, bins_wanted = 100, file_type=graph_file_type, save_fig=True)


    if os.path.exists(f'{save_file_name}.pickle'):
        os.remove(f'{save_file_name}.pickle')
    # remove .fasta from any files for pickle naming
    with open(f'{save_file_name}.pickle', 'wb') as handle:
        pickle.dump(test_seqs_dict, handle)

    end = time.perf_counter()
    round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    print('########################################################\n')
    return f'{save_file_name}.pickle'


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


def find_cat_sites(target_sequence: TargetSeq, igs_length: int = 5, guide_length: int = 50, min_length: int = 35):
    """
    Finds all instances of a U or T in a set of sequences and creates ribozyme designs for these sites.

    :param target_sequence:
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


def adjust_var_regs(known_seq_file: str, known_var_regs: list | int, unknown_seq_file: str):
    known_seq = read_fasta(known_seq_file)[0][1]
    unknown_seq = read_fasta(unknown_seq_file)[0][1]
    aligner = PairwiseAligner(mode='global')
    alignments = aligner.align(unknown_seq, known_seq)[0]

    # seq_a is the test sequence, seq_b is the reference sequence
    unknown_seq_aligned, known_seq_aligned = alignments

    if type(known_var_regs) == list:
        unknown_var_regs = []
        for idxs in known_var_regs:
            temp_idx_for_var_reg = []
            for unaligned_current_idx in idxs:
                # extract the position that is known
                current_spot = unaligned_current_idx
                known_current_idx = len(known_seq_aligned[:current_spot].replace('-', ''))
                while known_current_idx != unaligned_current_idx:
                    current_spot += 1
                    known_current_idx = len(known_seq_aligned[:current_spot].replace('-', ''))

                # now use this info to find what the corresponding unknown sequence index is
                unknown_current_idx = len(unknown_seq_aligned[:current_spot].replace('-', ''))
                temp_idx_for_var_reg.append(unknown_current_idx)
            unknown_var_regs.append((temp_idx_for_var_reg[0], temp_idx_for_var_reg[1]))
        return unknown_var_regs
    else:
        current_spot = known_var_regs
        known_current_idx = len(known_seq_aligned[:current_spot].replace('-', ''))
        while known_current_idx != known_var_regs:
            current_spot += 1
            known_current_idx = len(known_seq_aligned[:current_spot].replace('-', ''))

        # now use this info to find what the corresponding unknown sequence index is
        unknown_current_idx = len(unknown_seq_aligned[:current_spot].replace('-', ''))
        return unknown_current_idx


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

    counts_of_igs = defaultdict(lambda: 0)
    counts_of_ref_idx = defaultdict(lambda: 0)
    igs_ids_counts = defaultdict(lambda: 0)

    # Extract all the IGS id numbers - that's the IGS sequence and the reference index number
    total_targets = aligned_targets.size
    with alive_bar(total_targets, spinner='fishes') as bar:
        for seq in aligned_targets:
            igses = set()
            ref_ids = set()
            igs_ids = set()
            for item in seq.cat_sites:
                igses.add(f'{item[2]}')
                ref_ids.add(item[1])
                igs_ids.add(f'{item[2]}{item[1]}')

            for igs in igses:
                counts_of_igs[igs] += 1
            for ref_id in ref_ids:
                counts_of_ref_idx[ref_id] += 1
            for igs_id in igs_ids:
                igs_ids_counts[igs_id] += 1
            bar()

    if return_test_seqs:
        # IGS_id: [guides, targets, perc cov, perc on target, true perc cov]
        igs_over_min_true_cov = {igs: [[], set(), counts_of_igs[igs[:igs_len]] / total_targets,
                                       counts / counts_of_igs[igs[:igs_len]], counts / total_targets, counts]
                                 for igs, counts in igs_ids_counts.items()}
        # Here each item in the list is an IGS id (IGS + reference index), a guide, and the location
        # of the target sequence
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
            output_dict[ref_id] = (output_dict[ref_id], count_of_ref / len(aligned_targets), count_of_ref)
        return output_dict, igs_ids_counts, counts_of_ref_idx
    else:
        # this gives us a dictionary where the ID is matched to the true percent coverage
        # Measure true percent coverage of these and keep IGSes that are at least the minimum
        # true percent coverage needed
        igs_over_min_true_cov = {igs: [[], set(), counts_of_igs[igs[:igs_len]] / total_targets,
                                       counts / counts_of_igs[igs[:igs_len]], counts / total_targets]
                                 for igs, counts in igs_ids_counts.items()
                                 if counts / total_targets >= min_true_cov}

        # Here each item in the list is an IGS id (IGS + reference index), a guide, and the location of the
        # target sequence in our initial aligned_targets array
        igs_subsets = [(f'{item[2]}{item[1]}', item[3], i) for i, seq in enumerate(aligned_targets) for item in
                       seq.cat_sites]

        print(f'{len(igs_over_min_true_cov)} putative designs found.')

        print('Filtering IGSes that meet our min true U_IGS coverage...')
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
                                               perc_cov_attr=item[2], perc_on_target_attr=item[3],
                                               true_perc_cov_attr=item[4])
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
            new_seq_score, new_seq = replace_ambiguity(sequence_to_fix=to_optimize, target_background=True,
                                                       opti_len=guide_len, n_limit=1, update_score=False)

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


def get_naive_score(opti_seq, opti_len):
    # This version penalizes sequences that are too short
    score = 1 - (str(opti_seq).count('N') + str(opti_seq).count('-')) / opti_len
    return score


def get_weighted_score(opti_seq, opti_len):
    # This version penalizes sequences that are too short
    # prepare scoring matrix a(x)
    a = give_scoring_dict()

    prelim_score = 0
    for i, x in enumerate(str(opti_seq)):
        prelim_score += a[x]
    score = prelim_score / opti_len

    return score


# noinspection DuplicatedCode
def get_directional_score(opti_seq, opti_len):
    # This version penalizes sequences that are too short
    # Note that the end of the sequence is closer to the catalytic U!
    a = give_scoring_dict()

    # score based on extended identity and position of the bases
    ideal_score = 0
    current_score = 0

    for i, x in enumerate(str(opti_seq)):
        w = (i + 1)  # Calculate weight based on proximity to IGS
        ideal_score += w
        current_score += w * a[x]
    prelim_score = current_score / ideal_score

    # Penalize sequences below the minimum length
    score = prelim_score * (len(opti_seq)) / opti_len

    return score


def couple_designs_to_test_seqs(designs_input: str | list[str], test_seqs_input: str, file_to_save: str = '',
                                flexible_igs: bool = False, igs_len: int = 5, ref_idx_u_loc: int = 0,
                                score_type: str = 'weighted'):
    # For now, you must tell the program the location of the U for random sequences to analyze.
    # First, make a results folder, checkpoint, and coupled folder if they do not exist
    subfile = os.path.normpath(file_to_save + '/coupled')
    if not os.path.exists(subfile):
        os.mkdir(subfile)

    # Check if we've pre_processed designs correctly
    if (designs_input[-7:] != '.pickle') | (test_seqs_input[-7:] != '.pickle'):
        if designs_input[-igs_len - 1] not in ['G', 'g']:
            print('Sequence does not seem to have a U in the indicated igs length. Please adjust.')
            return
        igs = str(Seq(designs_input[-igs_len:]).upper().back_transcribe())
        guide = Seq(designs_input[:-igs_len - 1]).upper().back_transcribe()
        igs_id = igs + str(ref_idx_u_loc)
        # Can get the score at least
        score = return_score_from_type(sequence_to_test=guide, score_type=score_type, opti_len=len(guide))
        designs = [RibozymeDesign(id_attr=igs_id, guide_attr=guide, igs_attr=igs, ref_idx_attr=ref_idx_u_loc,
                                  score_type_attr=score_type, score_attr=score)]

        # Rename designs_input to correctly name file later on
        designs_input = os.path.normpath(subfile + '/' + igs_id + '.pickle')
    else:
        # Extract design sequence data
        with open(designs_input, 'rb') as handle:
            designs = pickle.load(handle)
    # Extract test sequences data
    with open(test_seqs_input, 'rb') as handle:
        test_seqs = pickle.load(handle)

    # Prepare data for check_guide_stats:
    # Need the design, the aligned target sequences to test against for each design
    coupled_datasets_name = designs_input.split('.')[0].split('/')[-1].split('\\')[-1] + '_vs_' + \
                            test_seqs_input.split('.')[0].split('/')[-1].split('\\')[-1]
    pickle_file_name = os.path.normpath(
        subfile + f'/{len(designs)}_designs_' + coupled_datasets_name + '.coupled').replace('_designs_designs_',
                                                                                            '_designs_')
    # Check if this file already exists
    if os.path.exists(pickle_file_name) and os.path.exists(os.path.normpath(subfile + '/big_checkpoint.txt')):
        print('Coupled file and checklist already exist! Skipping coupling...')
        return
    print(f'Coupling designs for {coupled_datasets_name}')
    with alive_bar(len(designs), spinner='fishes') as bar:
        for design in designs:
            matching_ref_idx_seqs = test_seqs[design.ref_idx]
            # If there are no test sequences with that index, set stats to zero and check the next one
            if len(matching_ref_idx_seqs) == 0:
                design.u_conservation_test = 0
                design.perc_cov_test = 0
                design.perc_on_target_test = 0
                design.true_perc_cov_test = 0
                design.number_of_targets_test = 0
                bar()
                continue
            # Set correct u conservation
            design.u_conservation_test = matching_ref_idx_seqs[1]

            # Set what targets have this IGS id
            # Recall each inner dictionary is of form IGS_id: [guides, targets, perc cov, perc on target, true perc cov]
            design_id = design.igs + str(design.ref_idx)
            try:
                design.test_targets = matching_ref_idx_seqs[0][design_id][1]
                # Set igs stats
                design.perc_cov_test = matching_ref_idx_seqs[0][design_id][2]
                design.perc_on_target_test = matching_ref_idx_seqs[0][design_id][3]
                design.true_perc_cov_test = matching_ref_idx_seqs[0][design_id][4]
                design.number_of_targets_test = matching_ref_idx_seqs[0][design_id][5]
            except KeyError:  # if there is not any matching IGSes
                design.perc_cov_test = 0
                design.perc_on_target_test = 0
                design.true_perc_cov_test = 0
                # recall that number of targets is how many targets have that IGS!
                design.number_of_targets_test = 0
                if not flexible_igs:
                    bar()
                    continue
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

    # save into a separate folder. Each pickle file will have 50 designs in it
    print('Now saving...')
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        with open(pickle_file_name, 'wb') as handle:
            pickle.dump(designs, handle)
        bar()

    print('Updating big checkpoint file...')
    file_name_to_check = os.path.normpath(subfile + '/big_checkpoint.txt')
    if not os.path.exists(file_name_to_check):
        with open(file_name_to_check, 'a') as d:
            for i in range(len(designs)):
                d.write(f'{i}\t{pickle_file_name}\t{i}\n')
    else:
        with open(file_name_to_check, 'r') as d:
            last_idx = int(d.readlines()[-1].split('\t')[0])
        with open(file_name_to_check, 'a') as d:
            for i in range(len(designs)):
                d.write(f'{i + last_idx + 1}\t{pickle_file_name}\t{i}\n')
    print('Done!')
    return pickle_file_name


def ribo_checker(coupled_folder: str, number_of_workers: int, worker_number: int, n_limit: int = 0,
                 display_bar: bool = True, opti_len: int = 50, get_tm_nn: bool = False):
    """
    This is the parallelization script I need to work on.
    """
    worker_number = int(worker_number)
    number_of_workers = int(number_of_workers)

    work_done_file = os.path.normpath(f'{coupled_folder}/work_done_{worker_number}.txt')
    designs_skipped = os.path.normpath(f'{coupled_folder}/skipped_designs{worker_number}.txt')
    results_folder = os.path.normpath(f'{coupled_folder}/results')

    # Check if we have anything to test
    analysis_files = [file for file in os.listdir(coupled_folder) if file.endswith('.coupled')]

    # Generate input results folder if it does not exist yet
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

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

    # Decide then how much work is allocated to each worker by dividing these evenly
    # Dictionary {file_name: [(big index, small index)]}
    this_worker_worklist = np.array_split(np.array(work_to_do_list, dtype=tuple), number_of_workers)[worker_number]

    files_to_open_dict = defaultdict(lambda: [])
    for big_idx, file_name, small_idx in this_worker_worklist:
        files_to_open_dict[file_name].append((int(big_idx), int(small_idx)))

    # Print how much work is left to do
    print(f'Total work:{total_work}\nWork to be done: {work_to_do}\nWork completed: {work_completed} '
          f'({round(work_completed / total_work * 100, 2)}% done)\n'
          f'Work for worker {worker_number}: {this_worker_worklist.shape[0]}\n')

    print('Opening files and extracting designs...')
    designs_to_test = []
    # Now only open the files relevant to this current worker and only extract the designs needed
    if display_bar:
        with alive_bar(len(files_to_open_dict.keys()), spinner='fishes') as bar:
            for file_name in files_to_open_dict.keys():
                with open(file_name, 'rb') as handle:
                    to_test_temp = pickle.load(handle)
                for big_idx, small_idx in files_to_open_dict[file_name]:
                    # Extract designs matching small index, save as (design, big index, file_name, small index)
                    designs_to_test.append((to_test_temp[small_idx], big_idx, file_name))
                bar()
    else:
        for file_name in files_to_open_dict.keys():
            with open(file_name, 'rb') as handle:
                to_test_temp = pickle.load(handle)
            for big_idx, small_idx in files_to_open_dict[file_name]:
                # Extract designs matching small index, save as (design, big index, file_name, small index)
                designs_to_test.append((to_test_temp[small_idx], big_idx, file_name))

    # Analyze designs
    print('Now analyzing designs...')
    if display_bar:
        with alive_bar(len(designs_to_test), spinner='fishes') as bar:
            for (design_to_test, big_idx, file_name) in designs_to_test:
                result = compare_to_test(design_to_test, n_limit=n_limit, test_dataset_name=file_name,
                                         guide_len=opti_len, get_tm_nn=get_tm_nn)
                naming_for_file = file_name.split('/')[-1].split('\\')[-1].split('.')[0] + f'_worker_{worker_number}'

                # If our result does not meet n_limit requirements, skip it
                if not result:
                    with open(work_done_file, 'a') as d:
                        d.write(str(big_idx) + '\n')
                    with open(designs_skipped, 'a') as d:
                        d.write(str(result.id) + '\n')
                    bar()
                    continue

                result_dict = result.to_dict(all_data=False)
                with open(os.path.normpath(f'{results_folder}/{naming_for_file}_results.txt'), 'a') as d:
                    d.write(json.dumps(result_dict) + '\n')

                with open(work_done_file, 'a') as d:
                    d.write(str(big_idx) + '\n')
                bar()
    else:
        for (design_to_test, big_idx, file_name) in designs_to_test:
            result = compare_to_test(design_to_test, n_limit=n_limit, test_dataset_name=file_name, guide_len=opti_len,
                                     get_tm_nn=get_tm_nn)
            naming_for_file = file_name.split('/')[-1].split('\\')[-1].split('.')[0] + f'_worker_{worker_number}'

            # If our result does not meet n_limit requirements, skip it
            if not result:
                with open(work_done_file, 'a') as d:
                    d.write(str(big_idx) + '\n')
                with open(designs_skipped, 'a') as d:
                    d.write(str(result.id) + '\n')
                continue

            result_dict = result.to_dict(all_data=False)
            with open(os.path.normpath(f'{results_folder}/{naming_for_file}_results.txt'), 'a') as d:
                d.write(json.dumps(result_dict) + '\n')

            with open(work_done_file, 'a') as d:
                d.write(str(big_idx) + '\n')
    return


def select_designs(tested_to_targets_path: list[str], designs_required: int, results_folder: str,
                   design_type: str = 'universal', igs_min: float = 0.7, guide_min: float = 0.7,
                   choose_var_reg_site: bool = False, start_idx: int = 1295, end_idx: int = 1435,
                   save_pickle: bool = True, file_extra_text: str = '', add_overhangs: bool = False,
                   igs_overhang: str = 'GGTCTCttttt', guide_overhang: str = 'gtcgtGAGACC'):
    """
    By default find things in E. coli conserved region between V8-V9.default is E. coli from:
    Chakravorty, S., Helb, D., Burday, M., Connell, N. & Alland, D. A detailed analysis of 16S ribosomal RNA gene
    segments for the diagnosis of pathogenic bacteria. J Microbiol Methods 69, 330-339 (2007).
    :return:
    """

    if design_type not in ['universal', 'selective']:
        print('Please indicate whether the designs are either universal or selective')
        return -1

    # Import tested designs
    # This function is being repurposed from the graph making functions pay no mind to how many nested lists there are
    with alive_bar(spinner='fishes') as bar:
        target_seqs_df = import_data_to_df(tested_to_targets_path, None)
        bar()

    # Sort by highest IGS difference then by highest guide difference
    if design_type == 'universal':
        igs_var = 'true_U_IGS_cov_test'
        guide_var = 'test_score'
    else:
        igs_var = 'delta_igs_vs_test'
        guide_var = 'delta_guide_vs_test'

    # Remove anything with an igs true coverage below the minimum asked
    filtered_df = target_seqs_df[target_seqs_df[igs_var] >= igs_min]

    # Remove anything with a guide score below the minimum asked
    filtered_df = filtered_df[filtered_df[guide_var] >= guide_min]

    # Remove those not in included variable regions if asked
    if choose_var_reg_site:
        filtered_df = filtered_df[(filtered_df['reference_idx'] <= end_idx) &
                                  (filtered_df['reference_idx'] >= start_idx)]
    # sort to select
    filtered_df.sort_values(by=[igs_var, guide_var], ascending=False)

    # Save designs into file
    print(f'{filtered_df.shape[0]} out of {target_seqs_df.shape[0]} designs meet given parameters! '
          f'Selecting top {designs_required}...\n')
    file_name = os.path.normpath(f'{results_folder}/{design_type}_top_{designs_required}_designs{file_extra_text}.csv')
    if os.path.exists(file_name):
        os.remove(file_name)
    filtered_df.head(designs_required).to_csv(file_name, index=False)
    # Go ahead and save all the designs above parameters if asked
    if save_pickle:
        file_name_pickle = file_name.split('.')[0] + '.pickle'
        filtered_df.to_pickle(file_name_pickle)
    print('File saved!')

    # Now get designs to order them easily
    igses = filtered_df.head(designs_required)['igs'].tolist()
    guides = filtered_df.head(designs_required)['guide'].tolist()
    to_order = []
    for igs, guide in zip(igses, guides):
        design = guide + 'g' + igs
        if add_overhangs:
            design = igs_overhang + design + guide_overhang
        to_order.append(design)
    return to_order, filtered_df.head(designs_required)


def compare_to_test(coupled_design: RibozymeDesign, n_limit, test_dataset_name, guide_len, get_tm_nn):
    def set_no_targets(ribo_design_attr: RibozymeDesign, n_limit_attr: int, test_dataset_name_attr: str,
                       guide_len_attr: int):
        # No hits in the background! filter out those that have too many Ns
        # Score is set at 0 based if there is no u conservation or guide sequences

        # remove designs with too many Ns
        if ribo_design_attr.guide.count('N') / guide_len_attr <= n_limit_attr:
            # Follow the same naming as the outputs of replace_ambiguity.
            ribo_design_attr.update_to_test(test_score_attr=0, name_of_test_dataset_attr=test_dataset_name_attr)
            return ribo_design_attr
        else:
            return None

    # See if we can skip computation all together because if no conserved u, no conserved IGS
    if coupled_design.u_conservation_test == 0 or len(coupled_design.guides_to_use) == 0:
        out = set_no_targets(coupled_design, n_limit, test_dataset_name, guide_len)
        return out

    # Find average guide score
    guides_to_optimize = np.array([str(guide) for guide in coupled_design.guides_to_use])
    fn_pairwise = np.vectorize(pairwise_comparison, otypes=[RibozymeDesign],
                               excluded=['seq_b', 'score_type', 'only_consensus', 'opti_len', 'get_tm_nn'])

    scores = fn_pairwise(guides_to_optimize, seq_b=coupled_design.guide, score_type=coupled_design.score_type,
                         only_consensus=True, opti_len=guide_len, get_tm_nn=get_tm_nn)
    if get_tm_nn:
        score = np.array([val[0] for val in scores]).mean()
        tm = np.array([val[1] for val in scores]).mean()
        coupled_design.update_to_test(test_score_attr=score, name_of_test_dataset_attr=test_dataset_name,
                                      test_tm_nn_attr=tm)
    else:
        coupled_design.update_to_test(test_score_attr=scores.mean(), name_of_test_dataset_attr=test_dataset_name)

    return coupled_design


def replace_ambiguity(sequence_to_fix: RibozymeDesign, opti_len: int, target_background: bool = False,
                      n_limit: float = 1, update_score: bool = True, get_tm_nn: bool = False):
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
        if get_tm_nn:
            new_seq_score, pairwise_score, tm = pairwise_comparison(seq_a=new_seq, seq_b=sequence_to_fix.anti_guide,
                                                                    score_type=sequence_to_fix.score_type,
                                                                    opti_len=opti_len, get_tm_nn=get_tm_nn)
            sequence_to_fix.update_to_background(background_score_attr=pairwise_score, new_guide=Seq(new_seq),
                                                 new_score=new_seq_score, reset_guides=True, tm_nn_attr=tm)
        else:
            new_seq_score, pairwise_score = pairwise_comparison(seq_a=new_seq, seq_b=sequence_to_fix.anti_guide,
                                                                score_type=sequence_to_fix.score_type,
                                                                opti_len=opti_len)
            sequence_to_fix.update_to_background(background_score_attr=pairwise_score, new_guide=Seq(new_seq),
                                                 new_score=new_seq_score, reset_guides=True)
        return
    else:
        new_seq_score = return_score_from_type(sequence_to_test=new_seq, score_type=sequence_to_fix.score_type,
                                               opti_len=opti_len)
        return new_seq_score, new_seq


def pairwise_comparison(seq_a: Seq | str, seq_b: Seq | str, score_type: str = 'naive', only_consensus: bool = False,
                        opti_len: int = None, priori: list = None, get_tm_nn: bool = False):
    # If optilen is not given, then will assume the optimum length is the length of the longest sequence

    if priori is None:
        priori = [0.25, 0.25, 0.25, 0.25]
    if len(seq_a) != len(seq_b):
        pad = '-' * abs(len(seq_a) - len(seq_b))
        if len(seq_a) > len(seq_b):
            seq_b += pad
        else:
            seq_a += pad

    mat = words2countmatrix([seq_a, seq_b], priori=priori)
    pairwise_comparison_consensus = consensus(mat, priori=priori)

    if get_tm_nn:
        # Getting an estimated Tm for the consensus sequence
        tm = get_tm_gc(pairwise_comparison_consensus, strict=False, mismatch=True, mismatch_base='-')

    if not opti_len:
        opti_len = len(pairwise_comparison_consensus)

    pairwise_score = return_score_from_type(sequence_to_test=pairwise_comparison_consensus, score_type=score_type,
                                            opti_len=opti_len)

    if only_consensus:
        if get_tm_nn:
            return pairwise_score, tm
        else:
            return pairwise_score
    else:
        seq_a_score = return_score_from_type(sequence_to_test=seq_a, score_type=score_type, opti_len=opti_len)
        if get_tm_nn:
            return seq_a_score, pairwise_score, tm
        else:
            return seq_a_score, pairwise_score


def write_output_file(designs: np.ndarray[RibozymeDesign] | list[RibozymeDesign], folder_path: str, taxonomy='',
                      all_data: bool = False):
    designs_df = pd.DataFrame.from_records([item.to_dict(taxonomy=taxonomy, all_data=all_data) for item in designs])
    if os.path.exists(f'{folder_path}.csv'):
        os.remove(f'{folder_path}.csv')
    designs_df.to_csv(f'{folder_path}.csv', index=False)
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


def muscle_msa_routine(sequences_to_align, name_of_file: str, muscle_exe_name: str, msa_fast: bool,
                       score_type: str, guide_len: int, priori=None):
    # First create a dummy file with all the guides we want to optimize
    if priori is None:
        priori = [0.25, 0.25, 0.25, 0.25]
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

    mat = words2countmatrix(alignments, priori=priori)
    opti_seq = consensus(mat, priori=priori)

    # Make sure our sequence is our desired length
    truncated_seq = opti_seq[-guide_len:]

    score = return_score_from_type(sequence_to_test=truncated_seq, score_type=score_type, opti_len=guide_len)

    return truncated_seq, score


def return_score_from_type(sequence_to_test: Seq | str, score_type: str, opti_len: int):
    if score_type == 'naive':
        # Naively tell me what the percent identity is: 1- ambiguity codons/ length
        score = get_naive_score(sequence_to_test, opti_len)

    elif score_type == 'weighted':
        score = get_weighted_score(sequence_to_test, opti_len)

    elif score_type == 'directional':
        score = get_directional_score(sequence_to_test, opti_len)

    else:
        print('Score type not recognized. Defaulting to naive scoring')
        score = get_naive_score(sequence_to_test, opti_len)
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


def give_scoring_dict():
    score_dict = {'A': 1, 'T': 1, 'G': 1, 'C': 1, 'U': 1, 'R': 0.5, 'Y': 0.5, 'M': 0.5, 'K': 0.5, 'S': 0.5, 'W': 0.5,
                  'H': 1 / 3, 'B': 1 / 3, 'D': 1 / 3, 'V': 1 / 3, 'N': 0, '-': 0}
    return score_dict


def combine_data(folder_path):
    # Get all files in the folder that end in .txt
    files = [os.path.normpath(f'{folder_path}/{f}') for f in os.listdir(folder_path)
             if f.endswith('.txt')]
    worker_file_names = defaultdict(lambda: [])

    for file_temp in files:
        file_name = file_temp.split('/')[-1].split('\\')[-1].split('_worker')[0]
        worker_file_names[os.path.normpath(f'{folder_path}/combined/{file_name}.txt')].append(file_temp)

    # Append the data from each file in a batch to combined name file
    file_to_test = os.path.normpath(f'{folder_path}/combined')
    if not os.path.exists(file_to_test):
        os.mkdir(file_to_test)
    else:
        for file in os.listdir(file_to_test):
            os.remove(os.path.normpath(f'{folder_path}/combined/{file}'))

    for new_file_name, items_to_write in worker_file_names.items():
        for item in items_to_write:
            with open(item, 'r') as r:
                designs = r.readlines()
            with open(new_file_name, 'a') as d:
                for line in designs:
                    d.write(line)
    return


"""
The functions below are from
    https://github.com/rsa-tools/rsat-code/blob/68048c6135224d5007abb388e2f64837ff57bc72/python-scripts/lib/matrix.py
Thank you so so so much for making this! For my reference here is the citation: 
https://github.com/rsa-tools/rsat-code/blob/master/CITATION.cff
"""


# Changed some variable names
def words2countmatrix(words, priori: list = None):
    """
    Convert a list of words to a simple count matrix.
    """
    if not priori:
        priori = [0.25, 0.25, 0.25, 0.25]
    w = len(words[0])
    m = [[0.0] * 4 for i in range(w)]
    n = len(words)
    for i in range(w):
        for s in range(n):
            try:
                letter = words[s][i]
            except IndexError:
                letter = 'N'
            if letter == 'A':
                m[i][0] = m[i][0] + 1.0
            elif letter == 'C':
                m[i][1] = m[i][1] + 1.0
            elif letter == 'G':
                m[i][2] = m[i][2] + 1.0
            elif letter == 'T':
                m[i][3] = m[i][3] + 1.0
            elif letter == 'N':
                m[i][0] = m[i][0] + priori[0]
                m[i][1] = m[i][1] + priori[1]
                m[i][2] = m[i][2] + priori[2]
                m[i][3] = m[i][3] + priori[3]
            # Below are my own edits
            # adding some ambiguity: in case the input sequence has an ambiguous base!
            elif letter == 'R':
                m[i][0] = m[i][0] + priori[0]/(priori[0] + priori[2])
                m[i][2] = m[i][2] + priori[2]/(priori[0] + priori[2])
            elif letter == 'Y':
                m[i][1] = m[i][1] + priori[1]/(priori[1] + priori[3])
                m[i][3] = m[i][3] + priori[3]/(priori[1] + priori[3])
            elif letter == 'W':
                m[i][0] = m[i][0] + priori[0]/(priori[0] + priori[3])
                m[i][3] = m[i][3] + priori[3]/(priori[0] + priori[3])
            elif letter == 'S':
                m[i][1] = m[i][1] + priori[1]/(priori[1] + priori[2])
                m[i][2] = m[i][2] + priori[2]/(priori[1] + priori[2])
            elif letter == 'M':
                m[i][0] = m[i][0] + priori[0]/(priori[0] + priori[1])
                m[i][1] = m[i][1] + priori[1]/(priori[0] + priori[1])
            elif letter == 'K':
                m[i][2] = m[i][2] + priori[2]/(priori[2] + priori[3])
                m[i][3] = m[i][3] + priori[3]/(priori[2] + priori[3])
            elif letter == 'H':
                m[i][0] = m[i][0] + priori[0]/(priori[0] + priori[1] + priori[3])
                m[i][1] = m[i][1] + priori[1]/(priori[0] + priori[1] + priori[3])
                m[i][3] = m[i][3] + priori[3]/(priori[0] + priori[1] + priori[3])
            elif letter == 'B':
                m[i][1] = m[i][1] + priori[1]/(priori[1] + priori[2] + priori[3])
                m[i][2] = m[i][2] + priori[2]/(priori[1] + priori[2] + priori[3])
                m[i][3] = m[i][3] + priori[3]/(priori[1] + priori[2] + priori[3])
            elif letter == 'V':
                m[i][0] = m[i][0] + priori[0]/(priori[0] + priori[1] + priori[2])
                m[i][1] = m[i][1] + priori[1]/(priori[0] + priori[1] + priori[2])
                m[i][2] = m[i][2] + priori[2]/(priori[0] + priori[1] + priori[2])
            elif letter == 'D':
                m[i][0] = m[i][0] + priori[0]/(priori[0] + priori[2] + priori[3])
                m[i][2] = m[i][2] + priori[2]/(priori[0] + priori[2] + priori[3])
                m[i][3] = m[i][3] + priori[3]/(priori[0] + priori[2] + priori[3])
    return m


'''
A                       (Adenine) 
C                       (Cytosine)
G                       (Guanine)
T                       (Thymine)
R       = A or G        (puRines)
Y       = C or T        (pYrimidines)
W       = A or T        (Weak hydrogen bonding)
S       = G or C        (Strong hydrogen bonding)
M       = A or C        (aMino group at common position)
K       = G or T        (Keto group at common position)
H       = A, C or T     (not G)
B       = G, C or T     (not A)
V       = G, A, C       (not T)
D       = G, A or T     (not C)
N       = G, A, C or T  (aNy)
'''

# Modified to include a gap
CONSENSUS = {'A': 'A',
             'C': 'C',
             'G': 'G',
             'T': 'T',
             'AG': 'R',
             'CT': 'Y',
             'AT': 'W',
             'CG': 'S',
             'AC': 'M',
             'GT': 'K',
             'ACT': 'H',
             'CGT': 'B',
             'ACG': 'V',
             'AGT': 'D',
             'ACGT': 'N',
             '': '-'
             }


# Changed some variable names
def consensus(matrix, priori: list = None, mask=False):
    w = len(matrix)

    if not priori:
        priori = [0.25, 0.25, 0.25, 0.25]

    str_list = []
    for i in range(w):
        letter = ''
        for b in range(4):
            if matrix[i][b] != 0.0 and log(matrix[i][b] / priori[b]) >= 0:
                letter = letter + 'ACGT'[b]
        if mask and (len(letter) == 4 or len(letter) == 3 or len(letter) == 2 or len(letter) == 0):
            letter = 'N'
        else:
            letter = CONSENSUS[letter]

        # if len(letter) != 1:
        #    letter = '[' + letter + ']'

        str_list.append(letter)
    return ''.join(str_list)


"""
The functions below is from biopython, modified to fit my needs!
"""


def get_tm_gc(seq, strict=True, valueset=7, userset=None, mismatch=True, mismatch_base='X'):
    """Return the Tm using empirical formulas based on GC content.

    General format: Tm = A + B(%GC) - C/N + salt correction - D(%mismatch)

    A, B, C, D: empirical constants, N: primer length
    D (amount of decrease in Tm per % mismatch) is often 1, but sometimes other
    values have been used (0.6-1.5). Use 'X' to indicate the mismatch position
    in the sequence. Note that this mismatch correction is a rough estimate.

    Arguments:
     - valueset: A few often cited variants are included:

        1. Tm = 69.3 + 0.41(%GC) - 650/N
           (Marmur & Doty 1962, J Mol Biol 5: 109-118; Chester & Marshak 1993),
           Anal Biochem 209: 284-290)
        2. Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch
           'QuikChange' formula. Recommended (by the manufacturer) for the
           design of primers for QuikChange mutagenesis.
        3. Tm = 81.5 + 0.41(%GC) - 675/N + 16.6 x log[Na+]
           (Marmur & Doty 1962, J Mol Biol 5: 109-118; Schildkraut & Lifson
           1965, Biopolymers 3: 195-208)
        4. Tm = 81.5 + 0.41(%GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x
           [Na+])) - %mismatch
           (Wetmur 1991, Crit Rev Biochem Mol Biol 126: 227-259). This is the
           standard formula in approximative mode of MELTING 4.3.
        5. Tm = 78 + 0.7(%GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+]))
           - %mismatch
           (Wetmur 1991, Crit Rev Biochem Mol Biol 126: 227-259). For RNA.
        6. Tm = 67 + 0.8(%GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+]))
           - %mismatch
           (Wetmur 1991, Crit Rev Biochem Mol Biol 126: 227-259). For RNA/DNA
           hybrids.
        7. Tm = 81.5 + 0.41(%GC) - 600/N + 16.6 x log[Na+]
           Used by Primer3Plus to calculate the product Tm. Default set.
        8. Tm = 77.1 + 0.41(%GC) - 528/N + 11.7 x log[Na+]
           (von Ahsen et al. 2001, Clin Chem 47: 1956-1961). Recommended 'as a
           tradeoff between accuracy and ease of use'.

     - userset: Tuple of four values for A, B, C, and D. Usersets override
       valuesets.
     - mismatch: If 'True' (default) every 'X' in the sequence is counted as
       mismatch.

    """
    if strict and any(x in seq for x in "KMNRYBVDH"):
        raise ValueError(
            "ambiguous bases B, D, H, K, M, N, R, V, Y not allowed when 'strict=True'"
        )

    # Ambiguous bases: add 0.5, 0.67 or 0.33% depending on G+C probability:
    percent_gc = SeqUtils.gc_fraction(seq, "weighted") * 100

    # gc_fraction counts X as 0.5
    if mismatch:
        percent_gc -= seq.count(mismatch_base) * 50.0 / len(seq)

    if userset:
        a, b, c, d = userset
    else:
        if valueset == 1:
            a, b, c, d = (69.3, 0.41, 650, 1)
        if valueset == 2:
            a, b, c, d = (81.5, 0.41, 675, 1)
        if valueset == 3:
            a, b, c, d = (81.5, 0.41, 675, 1)
        if valueset == 4:
            a, b, c, d = (81.5, 0.41, 500, 1)
        if valueset == 5:
            a, b, c, d = (78.0, 0.7, 500, 1)
        if valueset == 6:
            a, b, c, d = (67.0, 0.8, 500, 1)
        if valueset == 7:
            a, b, c, d = (81.5, 0.41, 600, 1)
        if valueset == 8:
            a, b, c, d = (77.1, 0.41, 528, 1)
    if valueset > 8:
        raise ValueError("allowed values for parameter 'valueset' are 0-8.")

    melting_temp = a + b * percent_gc - c / len(seq)

    if mismatch:
        melting_temp -= d * (seq.count(mismatch_base) * 100.0 / len(seq))
    return melting_temp


def extract_info(results_file_path: str, dataset: str):
    with open(results_file_path, 'r') as read_file:
        list_of_designs = read_file.readlines()
    # Now extact data here:
    column_types = {'id': str, 'igs': str, 'reference_idx': int, 'optimized_to_targets': bool,
                    'optimized_to_background': bool, 'tested': bool, 'tested_design': bool, 'guide': str,
                    'num_of_targets': int, 'score_type': str, 'score': float, 'U_IGS_coverage': float, 'U_IGS_on target': float,
                    'true_U_IGS_cov': float, 'composite_score': float, 'num_of_targets_background': int,
                    'u_conservation_background': float, 'background_score': float, 'U_IGS_coverage_background': float,
                    'U_IGS_on target_background': float, 'true_U_IGS_cov_background': float,
                    'composite_background_score': float,
                    'delta_igs_vs_background': float, 'delta_guide_vs_background': float, 'delta_vs_background': float,
                    'name_of_test_dataset': str, 'num_of_targets_test': int, 'u_conservation_test': float,
                    'test_score': float, 'tm_nn_vs_test': float, 'U_IGS_coverage_test': float, 'U_IGS_on target_test': float,
                    'true_U_IGS_cov_test': float, 'composite_test_score': float, 'delta_igs_vs_test': float,
                    'delta_guide_vs_test': float, 'delta_vs_test': float, 'target_Domain': str,
                    'target_Domain_background': str, 'target_Domain_test': str, 'target_Phylum': str,
                    'target_Phylum_background': str, 'target_Phylum_test': str, 'target_Class': str,
                    'target_Class_background': str, 'target_Class_test': str, 'target_Order': str,
                    'target_Order_background': str, 'target_Order_test': str, 'target_Family': str,
                    'target_Family_background': str, 'target_Family_test': str, 'target_Genus': str,
                    'target_Genus_background': str, 'target_Genus_test': str, 'target_Species': str,
                    'target_Species_background': str, 'target_Species_test': str, 'target_Taxon': str,
                    'target_Taxon_background': str, 'target_Taxon_test': str}

    designs = []
    for design in list_of_designs:
        design_dict = {}
        individual_attributes = design.strip(' {').split(', "')
        for key_val in individual_attributes:
            key, val = key_val.split('": ')
            key = key.strip('"')
            try:
                val_typed = column_types[key](val.strip('"\''))
            except ValueError:
                if column_types[key] == int or column_types[key] == float:
                    val_typed = column_types[key](0)
                else:
                    print(column_types[key])
            design_dict[key] = val_typed
        designs.append(design_dict)
    if dataset is not None:
        design_df = pd.DataFrame.from_records(designs, index=[dataset] * len(designs)).convert_dtypes()
    else:
        design_df = pd.DataFrame.from_records(designs, index=dataset).convert_dtypes()
    return design_df

def import_data_to_df(designs_path: list[str], name: str = 'designs', check_integrity: bool = True):
    dfs = []
    all_dupes_check = defaultdict(lambda: [])
    all_without_dupes = defaultdict(lambda: [])

    for i in range(0, len(designs_path)):
        temp = extract_info(designs_path[i], name)
        file_name = designs_path[i].split('/')[-1].split('\\')[-1].split('_worker')[0]
        if check_integrity:
            dupes_check = len(temp.index)
            all_dupes_check[file_name].append(dupes_check)
            without_dupes = len(temp.drop_duplicates().index)
            all_without_dupes[file_name].append(without_dupes)
            if dupes_check != without_dupes:
                print(f'\n{name} {file_name} has duplicate outputs!')
            dfs.append(temp)

    for key, val in all_without_dupes.items():
        all_vals = sum(val)
        design_amt = int(key.split('_')[0])
        if design_amt != all_vals:
            print(f'\n{name} {key} designs file is not done cooking! Missing {design_amt - all_vals} designs. '
                  f'Run again to complete analysis.')

    designs_df = pd.concat(dfs)
    return designs_df


def make_guide_score_plot(xdata: list, xlabel: str, ydata: list, ylabel: str, loc_data: list,
                          var_regs: list, save_file_name: str, bins_wanted: int = 100, file_type='png', save_fig=True,
                          add_control=False, control_x_data=None, control_y_data=None, control_loc_data=None, lim=1580,
                          dual_plot: bool = True):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (10.193, 6.932)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    if dual_plot:
        jointplot_fig = plt.figure()
        gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=14)
        joint_ax = {
            0: jointplot_fig.add_subplot(gridspec[0:6, 6:14]),
            1: jointplot_fig.add_subplot(gridspec[0:6, 0:6]),
        }
        # Plot scatter and kde plots
        # below is for individual log norms
        good_guides = 0
        for x, y in zip(xdata, ydata):
            if x > 0.75 and y > 0.75:
                good_guides += 1

        jointplot_fig.axes[0].axhline(y=0.75, color='orange', linestyle='-')
        jointplot_fig.axes[0].axvline(x=0.75, color='orange', linestyle='-')
        hist_data = jointplot_fig.axes[0].hist2d(xdata, ydata, cmap='viridis', bins=bins_wanted, norm='log')
        jointplot_fig.colorbar(hist_data[-1], ax=jointplot_fig.axes[0])

        jointplot_fig.axes[0].annotate(f'n = {good_guides}\ncandidate\ndesigns', xy=(0.75, 0.15), xycoords='axes fraction')

        jointplot_fig.axes[0].set_xlabel(xlabel)
        jointplot_fig.axes[0].set_ylabel(ylabel)
        # Set variable regions for location plots
        plot_variable_regions(joint_ax[1], var_regs)
        # Plot test data for each testing condition
        sns.scatterplot(x=loc_data, y=ydata, ax=joint_ax[1], alpha=0.3, legend=False, linewidth=0, size=0.5,
                        color='#000000')

        if add_control:
            sns.scatterplot(x=control_loc_data, y=control_y_data, ax=joint_ax[1], legend=False, alpha=1, c='#43a2ca',
                            edgecolors = 'black', linewidth = 0.8, marker = '^', sizes = [150] * len(control_loc_data))
            sns.scatterplot(x=control_x_data, y=control_y_data, ax=joint_ax[0], legend=False, alpha=1, c='#43a2ca',
                            edgecolors='black', linewidth=0.8, marker='^', sizes=[150] * len(control_loc_data))

        jointplot_fig.axes[1].set_xlabel('16s rRNA sequence position on reference sequence')
        jointplot_fig.axes[1].set_ylabel(ylabel)
        jointplot_fig.axes[0].set(xlim=[0, 1.02], ylim=[0, 1.02], ylabel=None)
        # recall e coli ref seq length is 1542, so 1580 should be plenty of space!
        jointplot_fig.axes[1].set(xlim=[-0.1, lim])
        jointplot_fig.axes[1].set(xlabel='Reference SSU rRNA index')
        jointplot_fig.axes[1].sharey(jointplot_fig.axes[0])
        jointplot_fig.axes[0].tick_params(labelleft=False)
    else:
        jointplot_fig = plt.figure()
        gridspec = jointplot_fig.add_gridspec(nrows=6, ncols=6)
        joint_ax = {
            0: jointplot_fig.add_subplot(gridspec[0:6, 0:6]),
        }
        # Plot scatter and kde plots
        # below is for individual log norms
        good_guides = 0
        for x, y in zip(xdata, ydata):
            if x > 0.75 and y > 0.75:
                good_guides += 1

        jointplot_fig.axes[0].axhline(y=0.75, color='orange', linestyle='-')
        jointplot_fig.axes[0].axvline(x=0.75, color='orange', linestyle='-')
        hist_data = jointplot_fig.axes[0].hist2d(xdata, ydata, cmap='viridis', bins=bins_wanted, norm='log')
        jointplot_fig.colorbar(hist_data[-1], ax=jointplot_fig.axes[0])

        jointplot_fig.axes[0].annotate(f'n = {good_guides}\ncandidate\ndesigns', xy=(0.75, 0.15),
                                       xycoords='axes fraction')

        jointplot_fig.axes[0].set_xlabel(xlabel)
        jointplot_fig.axes[0].set_ylabel(ylabel)

        if add_control:
            sns.scatterplot(x=control_x_data, y=control_y_data, ax=joint_ax[0], legend=False, alpha=1, c='#43a2ca',
                            edgecolors='black', linewidth=0.8, marker='^', sizes=[150] * len(control_loc_data))

        jointplot_fig.axes[0].set_ylabel(ylabel)
        jointplot_fig.axes[0].set(xlim=[0, 1.02], ylim=[0, 1.02])

    name = save_file_name.split('/')[-1].split('\\')[-1].replace('_', ' ').replace('scores', 'guide scores')
    jointplot_fig.suptitle(name.replace('vs test', 'tested to target'))
    plt.tight_layout()
    if save_fig:
        save_file_name = os.path.normpath(f'{save_file_name}_guide_scores.{file_type}')
        plt.savefig(fname=save_file_name, format=file_type)
    plt.show()

    return


def plot_variable_regions(ax, var_regs, invert=False):
    if invert:
        for V in var_regs:
            ax.axhspan(V[0], V[1], facecolor='g', alpha=0.2)
    else:
        for V in var_regs:
            ax.axvspan(V[0], V[1], facecolor='g', alpha=0.2)
    return


def make_test_seqs_graph(title: str, x_data: list, xlabel: str, y_data: list, ylabel: str,
                         loc_data: list[int], var_regs: list, save_file_name: str, alpha=0.5,
                         file_type: str = 'png', dataset_len: int = None, add_three_panel: bool = True,
                         add_control=False, control_x_data=None, control_y_data=None, control_loc_data=None):
    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.6, 18 * 0.6)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # Prepare axes for first figure
    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=7, ncols=14)
    joint_ax = {
        0: jointplot_fig.add_subplot(gridspec[1:7, 0:7]),
        1: jointplot_fig.add_subplot(gridspec[0:1, 0:7]),
        2: jointplot_fig.add_subplot(gridspec[1:7, 7:14]),
        3: jointplot_fig.add_subplot(gridspec[0:1, 7:14])
    }
    # Prepare hue data: color based on yes no var regs
    yes_no_var_reg = ['Conserved region'] * len(loc_data)
    for i, loc in enumerate(loc_data):
        for range_min, range_max in var_regs:
            if range_min <= loc <= range_max:
                yes_no_var_reg[i] = 'Variable region'
                break
    # Plot scatter and kde plots
    # color based on location:
    yes_no_palette = {'Conserved region': '#000000', 'Variable region': 'g'}
    sns.scatterplot(x=x_data, y=y_data, linewidth=0, alpha=alpha, hue=loc_data, ax=joint_ax[0], palette='viridis')
    norm = plt.Normalize(min(loc_data), max(loc_data))
    colorbar_data = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    colorbar_data.set_array([])
    joint_ax[0].get_legend().remove()
    jointplot_fig.colorbar(colorbar_data, ax=joint_ax[0], pad=0.1, orientation='horizontal',
                           label='Location along E. coli 16s')
    sns.scatterplot(x=x_data, y=y_data, linewidth=0, alpha=alpha, hue=yes_no_var_reg, ax=joint_ax[2],
                    palette=yes_no_palette)
    cmap = ListedColormap(colors=['#000000', 'g'])
    bounds = [0, 0.5, 1]
    norm = BoundaryNorm(bounds, cmap.N)
    colorbar_data = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    colorbar_data.set_array([])
    joint_ax[2].get_legend().remove()

    formatter = plt.FixedFormatter(['Conserved region', 'Variable region'])
    plt.colorbar(colorbar_data, ticks=plt.FixedLocator([0.25, 0.75]), format=formatter, ax=joint_ax[2], pad=0.1,
                 orientation='horizontal', label='Variable region location')

    sns.kdeplot(x=x_data, ax=joint_ax[1], fill=True, common_norm=True, alpha=.3, legend=False, color='#000000', cut=0,
                clip=[0.1, 1])
    joint_ax[1].label_outer()
    # sns.kdeplot(y=y_data, ax=joint_ax[2], fill=True, common_norm=True, alpha=.3, legend=False, color='#000000', cut=0)
    # jointplot_fig.axes[1].tick_params(axis='both', labelleft=False)
    sns.kdeplot(x=x_data, ax=joint_ax[3], hue=yes_no_var_reg, fill=True, common_norm=True, alpha=.3, legend=False,
                palette=yes_no_palette, cut=0, clip=[0.1, 1])
    if add_control:
        jointplot_fig.axes[0].scatter(x=control_x_data, y=control_y_data, alpha=1, c='#43a2ca', edgecolors='black',
                                      linewidth=0.8, marker='^', sizes=[150] * len(control_loc_data))
        jointplot_fig.axes[2].scatter(x=control_x_data, y=control_y_data, alpha=1, c='#43a2ca', edgecolors='black',
                                      linewidth=0.8, marker='^', sizes=[150] * len(control_loc_data))
    # sns.kdeplot(y=y_data, ax=joint_ax[5], hue=yes_no_var_reg, fill=True, common_norm=True, alpha=.3, legend=False,
    #             palette=yes_no_palette, cut=0)
    jointplot_fig.axes[0].set(xlabel=xlabel, ylabel=ylabel, xlim=[0, 1], ylim=[0, 1])
    jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    joint_ax[2].label_outer()
    joint_ax[3].label_outer()
    jointplot_fig.axes[1].tick_params(axis='both', labelleft=False, labelbottom=False, left=False)
    jointplot_fig.axes[2].set(xlabel=xlabel, xlim=[0, 1], ylim=[0, 1])
    jointplot_fig.axes[3].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[3].tick_params(axis='both', labelleft=False, labelbottom=False, left=False)
    jointplot_fig.axes[2].tick_params(axis='y', labelleft=False)
    if dataset_len is not None:
        jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} distributions '
                                                         f'n = {dataset_len}')
    else:
        jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} distributions')
    plt.tight_layout()
    plt.savefig(fname=f'{save_file_name}_var_regs.{file_type}', format=file_type)
    plt.show()

    # Prepare axes for second figure
    plot_three_panel_graph(var_regs, loc_data, x_data, alpha, y_data, xlabel, ylabel, dataset_len, title,
                           save_file_name, file_type, add_control=add_control, control_x_data=control_x_data,
                           control_y_data=control_y_data, control_loc_data=control_loc_data)

    return



def plot_three_panel_graph(var_regs, loc_data, x_data, alpha, y_data, xlabel, ylabel, dataset_len, title,
                           save_file_name, file_type, add_control: bool = False, control_loc_data=None,
                           control_x_data=None, control_y_data=None):
    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (30 * 0.6, 18 * 0.6)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')
    # Prepare axes for figure
    jointplot_fig = plt.figure()
    gridspec = jointplot_fig.add_gridspec(nrows=9, ncols=7)
    joint_ax = {
        0: jointplot_fig.add_subplot(gridspec[0:3, 0:7]),
        1: jointplot_fig.add_subplot(gridspec[3:6, 0:7]),
        2: jointplot_fig.add_subplot(gridspec[6:9, 0:7])
    }
    # Set variable regions for location plots
    plot_variable_regions(joint_ax[0], var_regs)
    plot_variable_regions(joint_ax[1], var_regs)
    plot_variable_regions(joint_ax[2], var_regs)
    # Plot test data for each testing condition
    unique_locs = {loc: [set(), 0] for loc in set(loc_data)}
    for loc, xval in zip(loc_data, x_data):
        unique_locs[loc][0].add(xval)
        unique_locs[loc][1] += 1
    to_plot_x = []
    to_plot_loc = []
    to_plot_nums_of_vals = []
    for loc, (val, num_of_vals) in unique_locs.items():
        to_plot_x.append(*val)
        to_plot_loc.append(loc)
        to_plot_nums_of_vals.append(num_of_vals)
    if len(to_plot_loc) != len(to_plot_x):
        print('x data contains non-unique values and cannot be plotted by reducing duplicates')
        sns.scatterplot(x=loc_data, y=x_data, linewidth=0, size=0.5, color='#000000',
                        ax=jointplot_fig.axes[0], alpha=alpha / 2, legend=False)
    else:
        sns.scatterplot(x=to_plot_loc, y=to_plot_x, linewidth=0, size=0.5, color='#000000',
                        ax=jointplot_fig.axes[0], alpha=1, legend=False)
        jointplot_fig.axes[0].set_title(f'{len(to_plot_loc)} unique U sites', loc='left')
        jointplot_fig.axes[0].label_outer()
    jointplot_fig.axes[1].bar(to_plot_loc, to_plot_nums_of_vals, color='#000000', edgecolor='#000000')
    # jointplot_fig.axes[1].set_title(f'{len(loc_data)} unique U-IGS sites', loc='left')
    jointplot_fig.axes[1].label_outer()
    sns.scatterplot(x=loc_data, y=y_data, linewidth=0, size=0.5, color='#000000',
                    ax=jointplot_fig.axes[2], alpha=alpha / 2, legend=False)
    if add_control:
        jointplot_fig.axes[0].scatter(x=control_loc_data, y=control_x_data, alpha=1, c='#43a2ca', edgecolors='black',
                                      linewidth=0.8, marker='^', sizes=[150] * len(control_loc_data))
        jointplot_fig.axes[2].scatter(x=control_loc_data, y=control_y_data, alpha=1, c='#43a2ca', edgecolors='black',
                                      linewidth=0.8, marker='^', sizes=[150] * len(control_loc_data))
    jointplot_fig.axes[2].set_title(f'{len(loc_data)} unique U-IGS sites', loc='left')
    # Set graph settings for pretti graphing
    jointplot_fig.axes[0].set_ylabel(xlabel)
    jointplot_fig.axes[1].set_ylabel('Number of U-IGS')
    jointplot_fig.axes[2].set_ylabel(ylabel)
    max_loc = max(loc_data)
    jointplot_fig.axes[0].set(xlabel=None, xlim=[-0.1, max_loc + 20], ylim=[-0.1, 1.1])
    jointplot_fig.axes[2].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[2].sharey(jointplot_fig.axes[0])
    jointplot_fig.axes[2].set(xlabel=None)
    jointplot_fig.axes[2].set(xlabel='Reference 16s rRNA index')
    jointplot_fig.axes[1].sharex(jointplot_fig.axes[0])
    jointplot_fig.axes[0].tick_params(labelbottom=False)
    if dataset_len is not None:
        jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} distributions '
                                                         f'n = {dataset_len}')
    else:
        jointplot_fig.suptitle(title.replace('_', ' ') + f' {xlabel} and {ylabel} distributions')
    plt.tight_layout()
    plt.savefig(fname=save_file_name + '.' + file_type, format=file_type)
    plt.show()
    return


def check_checkpoint_file(coupled_folder):
    print('Checking that there are files to analyze or if checkpoint file is corrupted...')
    # Check if we have anything to test
    analysis_files = [file for file in os.listdir(coupled_folder) if file.endswith('.coupled')]

    if len(analysis_files) == 0:
        print('Please make sure to couple designs with the appropriate test sequences')
        return -1

    # check the amount of work by summing the number of designs on each file in the coupled folder
    lengths = [int(name.split('/')[-1].split('\\')[-1].split('_designs')[0]) for name in analysis_files]
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
        if total_work == work_completed:
            print('All work to do has been done, enjoy analyzing your data!')
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
                         var_regs=None, graph_type: str = 'png', get_consensus_batches_test: bool = True,
                         get_consensus_batches_designs: bool = False,
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
                                  get_consensus_batches=get_consensus_batches_test, batch_num=batch_num,
                                  score_type=score_type, msa_fast=msa_fast,
                                  remove_x_dupes_in_graph=remove_x_dupes_in_graph, lim=var_regs_lim)
        else:
            print(f'{title} exists already! Moving on...')
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
                             gaps_allowed=gaps_allowed, get_consensus_batches=get_consensus_batches_designs,
                             random_guide_sample_size=random_guide_sample_size)
        else:
            print(f'{target_title} exists already! Moving on...')

        for test_outfile in all_test_file_names:
            _ = couple_designs_to_test_seqs(designs_input=target_save_file_name, test_seqs_input=test_outfile,
                                            flexible_igs=flexible_igs, file_to_save=out_path)
    return all_target_file_names, all_test_file_names
