import glob
import os
import time
import numpy as np
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy.random import default_rng
from Bio.Seq import Seq
import Bio.motifs
from collections import defaultdict
import subprocess
from Bio.Align import AlignInfo, PairwiseAligner
from Bio import AlignIO
from math import exp, log
import seaborn as sns
import matplotlib.pyplot as plt
import figure_plotting_functions_v2 as fpf
import math
from alive_progress import alive_bar
import matplotlib.patches as mpatches
from cmcrameri import cm


# Expanded version of SilvaSequence
class TargetSeq:
    # Important: TargetSeq id must be formatted like a Silva sequence if we want the give_taxonomy function to work!
    def __init__(self, id_attr: str = '', full_name_attr: str = '', seq_attr: Seq = '', cat_sites_attr: list = []):
        self.id = id_attr
        self.full_name = full_name_attr
        self.seq = seq_attr
        self.cat_sites = cat_sites_attr

    def __str__(self):
        return f'{self.id}({self.seq})'

    def __repr__(self):
        return f'{self.id}, cat_sites={self.cat_sites}'

    def give_taxonomy(self, level: str):
        level_names = {'Domain': 0, 'Phylum': 1, 'Class': 2, 'Order': 3,
                       'Family': 4, 'Genus': 5, 'Species': 6, 'Taxon': 7}
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

    def __init__(self, id_attr: str, guides_to_use_attr: list[Seq] = None, targets_attr: set() = None,
                 guide_attr: Seq = '', score_attr: int = None, score_type_attr: str = '', perc_cov_attr: float = None,
                 perc_on_target_attr: float = None, true_perc_cov_attr: float = None,
                 background_score_attr: float = None, perc_cov_background_attr: float = None,
                 perc_on_target_background_attr: float = None, true_perc_cov_background_attr: float = None,
                 background_guides_attr: list[Seq] = None, anti_guide_attr: Seq = '', anti_guide_score_attr: int = None,
                 background_targets_attr: set = None, igs_attr: str = '', ref_idx_attr: int = None,
                 u_consv_background_attr: float = None, tested_design_attr: bool = False):
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
        self.perc_cov = perc_cov_attr
        self.perc_on_target = perc_on_target_attr
        self.true_perc_cov = true_perc_cov_attr
        if not true_perc_cov_attr or perc_cov_attr or perc_on_target_attr:
            self.calc_percent_coverages(perc_cov_attr=perc_cov_attr, perc_on_target_attr=perc_on_target_attr,
                                        true_perc_cov_attr=true_perc_cov_attr)

        # If we initialized using an optimized design, set these attributes now.
        if self.score:
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
        self.perc_cov_background = perc_cov_background_attr
        self.perc_on_target_background = perc_on_target_background_attr
        self.true_perc_cov_background = true_perc_cov_background_attr
        self.calc_background_percent_coverages(perc_cov_background_attr=perc_cov_background_attr,
                                               perc_on_target_background_attr=perc_on_target_background_attr,
                                               true_perc_cov_background_attr=true_perc_cov_background_attr)

        if self.true_perc_cov_background and background_score_attr:
            self.composite_background_score = self.true_perc_cov_background * background_score_attr
        else:
            self.composite_background_score = None

        if self.composite_score and self.composite_background_score:
            self.delta_vs_background = self.composite_score - self.composite_background_score
        elif self.composite_score:
            self.delta_vs_background = self.composite_score
        else:
            self.delta_vs_background = 0

        if self.composite_background_score:
            self.optimized_to_background = True
        else:
            self.optimized_to_background = False

    def to_dict(self, taxonomy: str = ''):
        if taxonomy:
            if self.targets:
                target_input = str(give_taxonomy(str(self.targets),
                                                 level=taxonomy)).replace(',', ';').replace('\'', '')
            else:
                target_input = ''
            if self.background_targets:
                background_target_input = str(give_taxonomy(str(self.background_targets),
                                                            level=taxonomy)).replace(',', ';').replace('\'', '')
            else:
                background_target_input = ''
        else:
            target_input = str(self.targets).replace(',', '|')
            background_target_input = str(self.background_targets).replace(',', '|')
        return {
            'id': self.id,
            'igs': self.igs,
            'reference_idx': self.ref_idx,
            'optimized_to_targets': self.optimized_to_targets,
            'optimized_to_background': self.optimized_to_background,
            'tested_design': self.tested_design,
            'guide': str(self.guide),
            'guides_to_use': self.guides_to_use,
            'targets': target_input,
            'num_of_targets': self.number_of_targets,
            'score_type': self.score_type,
            'score': self.score,
            '%_coverage': self.perc_cov,
            '%_on target': self.perc_on_target,
            'true_%_cov': self.true_perc_cov,
            'composite_score': self.composite_score,
            'background_targets': background_target_input,
            'num_of_targets_background': self.number_of_targets_background,
            'background_guides': self.background_guides,
            'u_conservation_background': self.u_conservation_background,
            'anti_guide': str(self.anti_guide),
            'anti_guide_score': self.anti_guide_score,
            'background_score': self.background_score,
            '%_coverage_background': self.perc_cov_background,
            '%_on target_background': self.perc_on_target_background,
            'true_%_cov_background': self.true_perc_cov_background,
            'composite_background_score': self.composite_background_score,
            'delta_vs_background': self.delta_vs_background
        }

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
        if perc_cov_attr and true_perc_cov_attr and perc_on_target_attr:
            self.perc_cov = perc_cov_attr
            self.true_perc_cov = true_perc_cov_attr
            self.perc_on_target = perc_on_target_attr
        else:
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
        if perc_on_target_background_attr is not None and true_perc_cov_background_attr is not None and \
                perc_on_target_background_attr is not None:
            self.perc_cov_background = perc_cov_background_attr
            self.true_perc_cov_background = true_perc_cov_background_attr
            self.perc_on_target_background = perc_on_target_background_attr
        else:
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
        self.delta_vs_background = self.composite_score - self.composite_background_score

    def print_attributes(self, taxonomy='Order'):
        if self.optimized_to_targets and self.optimized_to_background:
            # first_row = 'IGS,Reference index,Number of species targeted,Score,Score type,% cov,% on target,true % cov,' \
            #             'Composite score,background score,% cov background,% on target background,' \
            #             'true % cov background,Composite score background,Delta composite score vs background,' \
            #             'Optimized guide,Optimized guide + G + IGS\n'
            text = f'{self.igs},{self.ref_idx},{self.number_of_targets},{self.score},{self.score_type},{self.perc_cov},' \
                   f'{self.perc_on_target},{self.true_perc_cov},{self.composite_score},{self.background_score},' \
                   f'{self.perc_cov_background},{self.perc_on_target_background},{self.true_perc_cov_background},' \
                   f'{self.composite_background_score},{self.delta_vs_background},{self.guide},' \
                   f'{self.guide}G{self.igs}\n'
        elif self.optimized_to_targets:
            # first_row = 'IGS,Reference index,Number of species targeted,Score,Score type,% cov,% on target,
            # True % cov,Composite score,Optimized guide,Optimized guide + G + IGS\n'
            text = f'{self.igs},{self.ref_idx},{self.number_of_targets},{self.score},{self.score_type},{self.perc_cov},' \
                   f'{self.perc_on_target},{self.true_perc_cov},{self.composite_score},' \
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
                 background_sequences_folder: str = '', min_delta: float = 0, optimize_seq: bool = True,
                 min_true_cov: float = 0.7, identity_thresh: float = 0.7, fileout: bool = False,
                 folder_to_save: str = '', score_type: str = 'quantitative', msa_fast: bool = False,
                 keep_single_targets: bool = False, gaps_allowed: bool = True, percent_of_target_seqs_used: float = 1.0,
                 percent_of_background_seqs_used: float = 1, n_limit: int = 0):
    """Generates ribozyme designs to target a set of sequences.
    :param percent_of_background_seqs_used: In case background data is very large we can get a random sample of the
    sequences used without replacement
    :param percent_of_target_seqs_used: In case target data is very large we can get a random sample of the sequences used
    without replacement
    :param gaps_allowed:
    :param min_delta: when making targeted designs, disregard any that have a composite score less than this
    :param barcode_seq_file: file containing the desired insertion barcode sequence (5' -> 3')
    :param ribobody_file: file containing the ribozyme sequence to use as the body template (5' -> 3')
    :param target_sequences_folder: folder containing all the sequences you want to design ribozymes for. Can be either
    a folder of fasta files or a single fasta file.
    :param igs_length: how long you want your IGS to be in base pairs. Default is 5 bp.
    :param guide_length: how long you want your guide sequence to be in base pairs. Default is 5 bp.
    :param min_length: minimum guide sequence length from 3' end. Must be smaller or equal to guide_length.
    Default is 35 bp. Ex: if you want your guide sequence to bind to at least 35 nt at the 3' end of the target
    sequence, set min_length = 35.
    :param ref_sequence_file:
    :param targeted:
    :param background_sequences_folder: folder containing all the sequences you do NOT want to design ribozymes for.
    Can be either a folder of fasta files or a single fasta file.
    :param optimize_seq: recommended. Uses MUSCLE multiple sequence alignment to generate optimized designs for several
    target sequences.
    :param min_true_cov: minimum percentage of targets you want to hit at a conserved location with a single optimized
    design. Default is 0.7 (70% of targets).
    :param identity_thresh: How much sequence identity do you want to use for the MSA. Only applicable for designs made
    without considering IUPAC ambiguity codes.
    :param fileout: whether we want a csv file output or not. Default is False.
    :param folder_to_save: the path where the folder we will save our outputs in if fileout = True
    :param score_type:
    :param msa_fast: whether to use super5 MUSCLE MSA or just regular MUSCLE MSA. Recommended for large datasets (over
    300 sequences) for faster data processing.
    :param keep_single_targets: whether we want to get designs that only hit one target sequence. This is overriden by
    min_true_cov - i.e. if min_true_cov is higher than the equivalent of one sequence then the program will not keep
    single targets.
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
    with alive_bar(unknown='fish') as bar:
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
    time1 = time.perf_counter()
    with alive_bar(unknown='fish') as bar:
        fn_msa = np.vectorize(optimize_designs, otypes=[RibozymeDesign],
                              excluded=['score_type', 'thresh', 'msa_fast', 'gaps_allowed'])
        fn_msa(to_optimize, score_type=score_type, thresh=identity_thresh, msa_fast=msa_fast, gaps_allowed=gaps_allowed)
        bar()
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='generating optimized designs')

    if fileout:
        write_output_file(designs=to_optimize, folder_path=folder_to_save)

    background_names_and_seqs = read_silva_fasta(in_file=background_sequences_folder)

    print(f'Found {background_names_and_seqs.size} total background sequences to analyze.')

    if percent_of_background_seqs_used < 1:
        background_names_and_seqs = np.random.choice(
            background_names_and_seqs, size=round(background_names_and_seqs.size * percent_of_target_seqs_used),
            replace=False)
        print(f'Randomly sampling {target_names_and_seqs.size} sequences to analyze.\n')

    # first find the IGSes and locations of the background sequences. We do not want to hit these.
    # Also align these to the reference sequence
    with alive_bar(unknown='fish') as bar:
        fn_align_to_ref(background_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length,
                        guide_length=guide_length, min_length=min_length)
        time2 = time.perf_counter()
        bar()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'finding catalytic sites and indexing background sequences to reference')

    print('Now applying designed ribozymes with background sequences and getting statistics...')
    opti_target_seqs = ribo_checker(to_optimize, background_names_and_seqs, len(ref_name_and_seq[1]),
                                    identity_thresh=identity_thresh, guide_length=guide_length,
                                    score_type=score_type, gaps_allowed=gaps_allowed, msa_fast=msa_fast,
                                    flexible_igs=True, n_limit=n_limit, target_background=not selective)

    if fileout:
        write_output_file(designs=opti_target_seqs, folder_path=folder_to_save)

    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4, task_timed='comparing designs against background '
                                                                      'sequences')

    end = time.perf_counter()
    round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    print('########################################################\n')
    return opti_target_seqs


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


def read_silva_fasta(in_file: str, file_type: str = 'fasta', exclude=[], include=[],
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
    min = hours % 1 * 60
    # Sec
    sec = min % 1 * 60

    if hours > 1:
        text = f'Time taken {task_timed}: {int(hours)} hrs, {int(min)} min, {round(sec, round_to)} sec\n'
    elif min > 1:
        text = f'Time taken {task_timed}: {int(min)} min, {round(sec, round_to)} sec\n'
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
    idx = [i for i, ltr in enumerate(target_sequence.seq) if ltr == 'T']
    # remove indexes that are < guide_length or > len(sequ) - igs_length (must have enough residues to attach to)
    idx_new = [res for res in idx if igs_length <= res < (len(target_sequence.seq) - min_length)]

    if not idx_new:
        print(f'No viable catalytic sites in {target_sequence.id}')
        target_sequence.cat_sites = None
    else:
        data = []
        for i in idx_new:
            # generate complementary guide sequence guide_length residues *downstream* of U site
            guide = target_sequence.seq[i + 1:i + guide_length + 1].reverse_complement()
            # generate complementary IGS sequence igs_length bases long *upstream* of U site
            igs = target_sequence.seq[i - igs_length:i].reverse_complement()
            # Store as (idx, igs, guide)
            data.append((i + 1, igs, guide))  # we're adding 1 to the idx because of 0 based indexing
        target_sequence.cat_sites = data


def align_to_ref(target_sequence: TargetSeq, ref_name_and_seq, igs_length: int = 5, guide_length: int = 50,
                 min_length: int = 35):
    """
    :param ref_name_and_seq: TargetSequence object
    :return:
    """
    if not target_sequence.cat_sites:
        find_cat_sites(target_sequence, igs_length, guide_length, min_length)
        if not target_sequence.cat_sites:
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

    # obtain index of new U
    idx_seq_a = [i for i, ltr in enumerate(seq_a) if ltr == 'T']

    data = []
    current_spot = 0
    for og_idx, igs, guide in target_sequence.cat_sites:
        seq_a_idx = len(seq_a[:idx_seq_a[current_spot]].replace('-', '')) + 1
        while seq_a_idx != og_idx:
            current_spot += 1
            seq_a_idx = len(seq_a[:idx_seq_a[current_spot]].replace('-', '')) + 1
        # find what index that is based on the reference sequence
        ref_string = seq_b[:idx_seq_a[current_spot]]
        ref_idx = len(ref_string.replace('-', '')) + 1  # turns zero based indexing to ref_seq numbering

        data.append((og_idx, ref_idx, igs, guide))

    target_sequence.cat_sites = data
    return target_sequence


def filter_igs_candidates(aligned_targets: np.ndarray[TargetSeq], min_true_cov: float = 0, igs_len: float = 5) -> None:
    """
    returns a dictionary of possible designs that meet our needs: checks each target sequence and filters the IGSes that
    work for the most targets. Returns a dictionary for each IGSid with the list of guides it needs to optimize, as
    well as the true percent coverage of these and the organism number that IGS exists in (this is the index in
    aligned_targets and we can use this later to calculate percent coverage and get taxonomy information of which
    sequences a particular target can hit.
    :param aligned_targets:
    :param min_true_cov:
    :return:
    """

    print('Finding IGS sequences...')

    all_igs_ids = []
    counts_of_igs = defaultdict(int)

    # Extract all the IGS id numbers - that's the IGS sequence and the reference index number
    with alive_bar(len(aligned_targets)) as bar:
        for seq in aligned_targets:
            igses = set()
            for item in seq.cat_sites:
                all_igs_ids.append(f'{item[2]}{item[1]}')
                igses.add(f'{item[2]}')

            for igs in igses:
                if counts_of_igs[igs]:
                    counts_of_igs[igs] += 1
                else:
                    counts_of_igs[igs] = 1
            bar()

    # Measure true percent coverage of these and keep IGSes that are at least the minimum true percent coverage needed
    igs_ids_counted, igs_ids_counts = np.unique(all_igs_ids, return_counts=True)

    # this gives us a dictionary where the ID is matched to the true percent coverage
    igs_over_min_true_cov = {igs: [[], set(), counts / aligned_targets.size] for igs, counts
                             in zip(igs_ids_counted, igs_ids_counts)
                             if counts / aligned_targets.size >= min_true_cov}

    # Here each item in the list is an IGS id (IGS + reference index), a guide, and the location of the target sequence
    # in our initial aligned_targets array
    igs_subsets = [(f'{item[2]}{item[1]}', item[3], i) for i, seq in enumerate(aligned_targets) for item in
                   seq.cat_sites]

    # get all the guides that have good IGSes
    print('Filtering IGSes that meet our min true % coverage...')
    with alive_bar(len(igs_subsets)) as bar:
        for igs_id, guide, target_num in igs_subsets:
            if igs_id in igs_over_min_true_cov:
                # Add guide to list
                igs_over_min_true_cov[igs_id][0].append(guide)
                # Add organism number to set. We can calculate the percent coverage with this number later on.
                igs_over_min_true_cov[igs_id][1].add(aligned_targets[target_num].full_name)
            bar()

    print(f'{len(igs_over_min_true_cov)} putative designs found.')
    # Now make an array of all of the putative designs for later use.
    to_optimize = np.array([RibozymeDesign(id_attr=igs_id, guides_to_use_attr=item[0], targets_attr=item[1],
                                           true_perc_cov_attr=item[2],
                                           perc_cov_attr=counts_of_igs[igs_id[:igs_len]] / aligned_targets.size)
                            for igs_id, item in igs_over_min_true_cov.items()])
    return to_optimize


def optimize_designs(to_optimize: RibozymeDesign, score_type: str, thresh: float = 0.7, guide_len: int = 50,
                     msa_fast: bool = True, gaps_allowed: bool = True, compare_to_background: bool = False):
    """
    Takes in a RibozymeDesign with guide list assigned and uses MUSCLE msa to optimize the guides.
    when compare_to_background is set to true, will do an MSA on background seqs and will update the anti-guide sequence

    """
    # First create a dummy file with all the guides we want to optimize
    if not compare_to_background:
        with open(f'to_align_{to_optimize.id}.fasta', 'w') as f:
            for line in range(len(to_optimize.guides_to_use)):
                f.write('>seq' + str(line) + '\n' + str(to_optimize.guides_to_use[line]) + '\n')
    else:
        with open(f'to_align_{to_optimize.id}.fasta', 'w') as f:
            for line in range(len(to_optimize.background_guides)):
                f.write('>seq' + str(line) + '\n' + str(to_optimize.background_guides[line]) + '\n')

    muscle_exe = 'muscle5'

    if msa_fast:
        subprocess.check_output([muscle_exe, '-super5', f'to_align_{to_optimize.id}.fasta', '-output',
                                 f'aln_{to_optimize.id}.afa'], stderr=subprocess.DEVNULL)
    else:
        subprocess.check_output([muscle_exe, '-align', f'to_align_{to_optimize.id}.fasta', '-output',
                                 f'aln_{to_optimize.id}.afa'], stderr=subprocess.DEVNULL)

    msa = AlignIO.read(open(f'aln_{to_optimize.id}.afa'), 'fasta')

    # delete files:
    if os.path.exists(f'aln_{to_optimize.id}.afa'):
        os.remove(f'aln_{to_optimize.id}.afa')
    if os.path.exists(f'to_align_{to_optimize.id}.fasta'):
        os.remove(f'to_align_{to_optimize.id}.fasta')

    # Now retrieve the alignment
    summary_align = AlignInfo.SummaryInfo(msa)
    alignments = [record.seq for record in summary_align.alignment]

    # thanks to https://github.com/mdehoon for telling me about .degenerate_consensus!!!
    if score_type != 'quantitative':
        if gaps_allowed:
            opti_seq = Bio.motifs.create(alignments, alphabet='GATCRYWSMKHBVDN-').degenerate_consensus.strip('-')
        else:
            opti_seq = Bio.motifs.create(alignments, alphabet='GATCRYWSMKHBVDN-').degenerate_consensus.strip(
                '-').replace('-', 'N')
        # Make sure our sequence is our desired length
        truncated_seq = opti_seq[-guide_len:]

        if score_type == 'naive':
            # Naively tell me what the percent identity is: 1- ambiguity codons/ length
            score = get_naive_score(truncated_seq)

        elif score_type == 'weighted':
            score = get_weighted_score(truncated_seq)

        elif score_type == 'directional':
            score = get_directional_score(truncated_seq)

        else:
            print('Score type not recognized. Defaulting to naive scoring')
            score = get_naive_score(truncated_seq)

    else:
        if gaps_allowed:
            left_seq = summary_align.gap_consensus(threshold=thresh, ambiguous='N')[-guide_len:]
        else:
            left_seq = summary_align.dumb_consensus(threshold=thresh, ambiguous='N')[-guide_len:]
        score, opti_seq = get_quantitative_score(left_seq, alignments, count_gaps=gaps_allowed)
        truncated_seq = opti_seq[-guide_len:]

    # Now update our initial design
    if not compare_to_background:
        to_optimize.update_after_optimizing(score_attr=score, guide_attr=truncated_seq, score_type_attr=score_type,
                                            reset_guides=True)
    else:
        # In case the sequence is using min_len that is different than the guide length
        to_optimize.anti_guide = opti_seq[-len(to_optimize.guide):]
        to_optimize.anti_guide_score = score


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


def ribo_checker(designs: np.ndarray[RibozymeDesign], aligned_target_sequences: np.ndarray[TargetSeq], ref_seq_len: int,
                 identity_thresh: float, guide_length: int, gaps_allowed: bool = True,
                 msa_fast: bool = False, flexible_igs: bool = False, target_background: bool = False,
                 n_limit: int = 0, score_type: str = 'naive'):
    """
    Checks generated designs against a set of sequences to either find how well they align to a given dataset.
    Will return a set of sequences that match for each design as well as a score showing how well they matched.
    :param target_background:
    :param designs:
    :param aligned_target_sequences:
    :param ref_seq_len:
    :param identity_thresh:
    :param guide_length:
    :param score_type:
    :param gaps_allowed:
    :param flexible_igs:
    :param msa_fast:
    :param n_limit: what percentage of Ns are acceptable?
    :return:
    """

    # This should update designs and return uracil_sites_targets
    uracil_sites_targets = check_cat_site_conservation(designs, aligned_target_sequences, ref_seq_len)

    print('Done! Checking IGS conservation...')
    igs_sites_targets, conserved_u_sites = check_igs_conservation(uracil_sites_targets, aligned_target_sequences)

    print('Done! Now checking finding guide stats for all U sites...')

    to_optimize, opti_seqs, designs_below_limit = check_guide_stats(designs, igs_sites_targets,
                                                                    aligned_target_sequences, flexible_igs,
                                                                    conserved_u_sites, guide_length, n_limit)

    print('Found guide stats. Now comparing designs to associated binding sites in background...')
    # Do a MSA. Remember that we are passing the references to objects, hopefully, if I coded this correctly.
    with alive_bar(unknown='fish') as bar:
        fn_msa = np.vectorize(optimize_designs, otypes=[RibozymeDesign],
                              excluded=['score_type', 'thresh', 'msa_fast', 'guide_len', 'gaps_allowed',
                                        'compare_to_background'])
        fn_msa(to_optimize, score_type=score_type, thresh=identity_thresh, msa_fast=msa_fast, guide_len=guide_length,
               gaps_allowed=gaps_allowed, compare_to_background=True)
        bar()

    to_reduce_ambiguity = np.append(opti_seqs, to_optimize)

    # Now reduce ambiguity and score
    print('\nNow scoring ribozyme designs against background sequences...')
    start = time.perf_counter()
    with alive_bar(unknown='fish') as bar:
        # Keep in mind that this will only update designs that are above the limit of ambiguity codons N given!
        fn_ambiguity = np.vectorize(replace_ambiguity, otypes=[RibozymeDesign],
                                    excluded=['target_background', 'n_limit', 'score_params'])
        fn_ambiguity(to_reduce_ambiguity, target_background=target_background,
                     n_limit=n_limit, score_params=[])

        less_ambiguous_designs_above_limit = np.append(
            designs_below_limit, [design for design in to_reduce_ambiguity if design.optimized_to_background])
        bar()

    round_convert_time(start=start, end=time.perf_counter(), round_to=4,
                       task_timed='scoring against background sequences')

    return less_ambiguous_designs_above_limit


def test_ribo_design(design: str, target_folder: str, ref_seq_folder: str, igs_len: int = 5, score_type: str = 'naive',
                     thresh: float = 0.7, msa_fast: bool = False, gaps_allowed: bool = False, file_out: bool = False,
                     folder_to_save: str = '', taxonomy='Order'):
    print('Testing ribozyme design!')
    start = time.perf_counter()
    # first, read in all our sequences
    targets_to_test = read_silva_fasta(in_file=target_folder)
    print(f'Found {targets_to_test.size} total target sequences to analyze.')
    ref_name_and_seq = read_fasta(ref_seq_folder)[0]
    igs = design[-igs_len:]
    design_guide = Seq(design[:-igs_len - 1]).upper().back_transcribe()  # recall there's a G between the igs and guide

    time1 = time.perf_counter()
    print(f'Now finding catalytic sites and aligning to reference {ref_name_and_seq[0].replace("_", " ")}...')
    # Align the targets to test to our reference sequence
    with alive_bar(unknown='fish') as bar:
        fn_align_to_ref = np.vectorize(align_to_ref, otypes=[TargetSeq],
                                       excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])
        fn_align_to_ref(targets_to_test, ref_name_and_seq=ref_name_and_seq, igs_length=igs_len,
                        guide_length=len(design_guide), min_length=len(design_guide))
        bar()
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'finding catalytic sites and indexing target sequences to reference')

    # get all guides and the name of their targets that have the correct igs: (og_idx, ref_idx, igs, guide)
    time1 = time.perf_counter()
    print('Finding matching IGSes and corresponding guides...')
    names = defaultdict(set)
    guides = defaultdict(list)
    all_names = set()
    with alive_bar(len(targets_to_test)) as bar:
        for target in targets_to_test:
            for putative_guide in target.cat_sites:
                if putative_guide[2] == igs:
                    names[putative_guide[1]].add(target.full_name)
                    guides[putative_guide[1]].append(putative_guide[3])
                    all_names.add(target.full_name)
            bar()
    num_of_seqs_with_igs = len(all_names)
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'matching IGSes and guides.')

    # Calculate coverage for each and make a ribozyme design for each
    time1 = time.perf_counter()
    print('Calculating coverage for putative locations targeted by guide...')
    locations_and_scores = np.array([])
    with alive_bar(len(guides.items())) as bar:
        for index, guide_list in guides.items():
            perc_cov = num_of_seqs_with_igs / targets_to_test.size
            true_perc_cov = len(guide_list) / targets_to_test.size
            per_on_target = perc_cov / perc_cov
            locations_and_scores = np.append(locations_and_scores,
                                             RibozymeDesign(id_attr=f'{igs}{index}', guide_attr=design_guide,
                                                            background_targets_attr=names[index],
                                                            background_guides_attr=guide_list,
                                                            score_type_attr=score_type,
                                                            perc_cov_background_attr=perc_cov,
                                                            perc_on_target_background_attr=per_on_target,
                                                            true_perc_cov_background_attr=true_perc_cov,
                                                            tested_design_attr=True))
            bar()

    fn_msa = np.vectorize(optimize_designs, otypes=[RibozymeDesign],
                          excluded=['score_type', 'thresh', 'msa_fast', 'guide_len', 'gaps_allowed',
                                    'compare_to_background'])
    fn_msa(locations_and_scores, score_type=score_type, thresh=thresh, msa_fast=msa_fast, guide_len=len(design_guide),
           gaps_allowed=gaps_allowed, compare_to_background=True)

    # Now score our design against these sequences
    with alive_bar(len(locations_and_scores)) as bar:
        for to_score in locations_and_scores:
            pairwise_comparison_consensus = Bio.motifs.create([to_score.guide, to_score.anti_guide],
                                                              alphabet='GATCRYWSMKHBVDN-').degenerate_consensus
            design_score, background_score = pairwise_comparison(to_score.guide, pairwise_comparison_consensus)
            to_score.background_guides = None
            to_score.score = design_score
            to_score.background_score = background_score
            to_score.composite_background_score = background_score * to_score.true_perc_cov_background
            bar()

    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'scoring design against background')

    if file_out:
        write_output_file(designs=locations_and_scores, folder_path=folder_to_save, taxonomy=taxonomy)

    end = time.perf_counter()
    round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    print('########################################################\n')
    return locations_and_scores


def check_cat_site_conservation(designs, aligned_target_sequences, ref_seq_len):
    designed_idxs = np.array(list(set(ribodesign.ref_idx - 1 for ribodesign in designs)))  # convert to base 0 indexing
    # Adding to maximum size because sometimes the alignment for the very end of a sequence does strange things
    max_len = max(designed_idxs.max(), ref_seq_len)
    uracil_sites_targets = np.zeros((aligned_target_sequences.size, max_len + 1))
    uracil_sites_targets[:, designed_idxs] = 1

    print('Finding U sites...')
    with alive_bar(len(aligned_target_sequences)) as bar:
        for i, target_sequence in enumerate(aligned_target_sequences):  # base 1 indexing
            temp_uracil_indexes = np.array([cat_site[1] - 1 for cat_site in
                                            target_sequence.cat_sites])  # convert to base 0 indexing to match array indexing locations
            # everything that is a 2 is a U that appears in both the designs and the target sequence,
            # a 1 is a U in the design
            uracil_sites_targets[i, temp_uracil_indexes] = uracil_sites_targets[
                                                               i, temp_uracil_indexes] * 2  # base 0 index
            bar()

    num_of_targets = aligned_target_sequences.size
    print('Found U sites. Checking how many are conserved...')
    with alive_bar(len(designed_idxs)) as bar:
        for ribodesign_num, i in enumerate(designed_idxs):
            u_values, u_counts = np.unique(uracil_sites_targets[:, i], return_counts=True)
            try:
                percentage_of_cat_sites_conserved = u_counts[np.where(u_values == 2)[0][0]] / num_of_targets
            except:
                percentage_of_cat_sites_conserved = 0
            designs[ribodesign_num].u_conservation_background = percentage_of_cat_sites_conserved
            bar()
    return uracil_sites_targets


def check_igs_conservation(uracil_sites_targets, aligned_target_sequences):
    # What is the conservation of IGSes at a particular position for each design?
    # Keep only the sites that are conserved for later.
    igs_sites_targets = np.array(np.clip(uracil_sites_targets, 0, 1), copy=True, dtype=str)  # base 0 indexing
    # Keeps the indexes of where Us are conserved, so base 0 indexing
    conserved_u_sites = np.argwhere(uracil_sites_targets == 2)

    # convert to base 0
    cat_sites_background = {(i, cat_site[1] - 1): str(cat_site[2]) for i, design in
                            enumerate(aligned_target_sequences) for cat_site in design.cat_sites}

    target_test = np.array(np.clip(uracil_sites_targets, 0, 1), copy=True, dtype=str)
    with alive_bar(len(conserved_u_sites)) as bar:
        for row, col in conserved_u_sites:
            target_test[row, col] = cat_sites_background[(row, col)]
            bar()

    return igs_sites_targets, conserved_u_sites


def check_guide_stats(designs, igs_sites_targets, aligned_target_sequences, flexible_igs, conserved_u_sites,
                      guide_length, n_limit):
    # Now find how many conserved IGSes there are!
    to_optimize = []
    designs_below_limit = []
    opti_seqs = []

    # ref_idx is in base 1 indexing, indexes for arrays are base 0
    with alive_bar(len(designs)) as bar:
        for ribo_design in designs:  # base 1 indexing
            check = np.argwhere(igs_sites_targets == ribo_design.igs)  # base 0 indexing
            # Fill this with empty TargetSeq objects to later make finding taxonomy more easy
            background_names = set(back_seq.full_name for back_seq in aligned_target_sequences[check[:, 0]])
            ribo_design.background_targets = background_names
            orgs_with_igs_on_target = check[np.argwhere(check[:, 1] == ribo_design.ref_idx - 1)][:, 0,
                                      0]  # convert to base 0 indexing
            on_target_count = orgs_with_igs_on_target.size
            if on_target_count > 0:
                true_perc_cov = on_target_count / aligned_target_sequences.size
                orgs_with_igs = np.unique(check[:, 0])
                perc_cov = orgs_with_igs.size / aligned_target_sequences.size
                perc_on_target = true_perc_cov / perc_cov
                ribo_design.calc_background_percent_coverages(perc_cov_background_attr=perc_cov,
                                                              true_perc_cov_background_attr=true_perc_cov,
                                                              perc_on_target_background_attr=perc_on_target)
                # If the igs is not flexible, only check those that have a perfect igs upstream of the cat_site
                if not flexible_igs:
                    # comparing two base 1 indexes, no conversion needed, but this will hold base 0
                    all_indexes_of_target_data = [np.argwhere(
                        np.array([cat_site[1] for cat_site
                                  in aligned_target_sequences[target_number].cat_sites]) == ribo_design.ref_idx)[0, 0]
                                                  for target_number in orgs_with_igs_on_target]

                    guides_to_optimize = [str(aligned_target_sequences[target].cat_sites[index_of_target_data][3]) for
                                          target, index_of_target_data in
                                          zip(orgs_with_igs_on_target, all_indexes_of_target_data)]
                # Otherwise, extract all info in conserved cat_sites regardless of igs. As long at the catalytic U is on
                # target, we'll extract those data
                else:
                    orgs_with_u_on_target = conserved_u_sites[
                                                np.where(conserved_u_sites[:, 1] == ribo_design.ref_idx - 1)][:,
                                            0]  # convert to base 0
                    # comparing two base 1 indexes, no conversion needed, but this will hold base 0
                    # I'm not sure why the following line won't work if I don't convert the list in np.argwhere into an array
                    all_indexes_of_target_data = [np.argwhere(
                        np.array([cat_site[1] for cat_site
                                  in aligned_target_sequences[target_number].cat_sites]) == ribo_design.ref_idx)[0, 0]
                                                  for target_number in orgs_with_u_on_target]
                    guides_to_optimize = [str(aligned_target_sequences[target].cat_sites[index_of_target_data][3]) for
                                          target, index_of_target_data in
                                          zip(orgs_with_u_on_target, all_indexes_of_target_data)]

                if len(guides_to_optimize) > 1:
                    # Pass the reference to this design to further optimize to the background
                    ribo_design.background_guides = guides_to_optimize
                    to_optimize.append(ribo_design)
                else:
                    # No need to optimize to the background with a full MSA but can still remove ambiguity
                    ribo_design.anti_guide = guides_to_optimize[0]
                    opti_seqs.append(ribo_design)
            else:
                # No hits in the background! filter out those that have too many Ns
                # Scores are set at 0 based on the assumption that if there is no catalytic U there is no splicing
                # activity
                ribo_design.calc_background_percent_coverages(perc_cov_background_attr=0,
                                                              true_perc_cov_background_attr=0,
                                                              perc_on_target_background_attr=0)
                # remove designs with too many Ns
                if ribo_design.guide.count('N') / guide_length <= n_limit:
                    # Follow the same naming as the outputs of replace_ambiguity.
                    ribo_design.update_to_background(background_score_attr=0, new_guide=ribo_design.guide,
                                                     new_score=ribo_design.score, reset_guides=True)
                    designs_below_limit.append(ribo_design)
            bar()

    return to_optimize, opti_seqs, designs_below_limit


def replace_ambiguity(sequence_to_fix: RibozymeDesign, target_background: bool = False,
                      n_limit: float = 1, score_params: list = []):
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

    # Now do a pairwise sequence with the target sequence to see the delta score:
    pairwise_comparison_consensus = Bio.motifs.create([new_seq, sequence_to_fix.anti_guide],
                                                      alphabet='GATCRYWSMKHBVDN-').degenerate_consensus

    new_seq_score, pairwise_score = pairwise_comparison(seq_a=new_seq, seq_b=pairwise_comparison_consensus,
                                                        score_type=sequence_to_fix.score_type,
                                                        score_params=score_params)

    sequence_to_fix.update_to_background(background_score_attr=pairwise_score, new_guide=new_seq,
                                         new_score=new_seq_score, reset_guides=True)


def pairwise_comparison(seq_a: Seq, seq_b: Seq, score_type: str = 'naive', score_params: list = []):
    # seq_b is pairwise
    if score_type == 'quantitative':
        if score_params:
            a_params = [seq_a].extend(score_params)
            seq_a_score = map(get_quantitative_score, a_params)
            b_params = [seq_b].extend(score_params)
            seq_b_score = map(get_quantitative_score, b_params)
        else:
            print('Cannot do quantitative score without a list of alignments. Performing weighted score instead...')
            seq_a_score = get_weighted_score(seq_a)
            seq_b_score = get_weighted_score(seq_b)
    elif score_type == 'weighted':
        seq_a_score = get_weighted_score(seq_a)
        seq_b_score = get_weighted_score(seq_b)
    elif score_type == 'directional':
        seq_a_score = get_directional_score(seq_a)
        seq_b_score = get_directional_score(seq_b)
    elif score_type == 'naive':
        seq_a_score = get_naive_score(seq_a)
        seq_b_score = get_naive_score(seq_b)
    else:
        print(f'Score type {score_type} is not supported. Please choose either weighted, quantitative, '
              f'or directional.')
        return -1

    return seq_a_score, seq_b_score


def write_output_file(designs: np.ndarray[RibozymeDesign], folder_path: str, taxonomy='Order'):
    if designs[0].optimized_to_targets and designs[0].optimized_to_background:
        if os.path.exists(f'{folder_path}/Targeted designs against background above threshold.csv'):
            os.remove(f'{folder_path}/Targeted designs against background above threshold.csv')
        first_row = 'IGS,Reference index,Number of species targeted,Score,Score_type,% cov,% on target,true % cov,' \
                    'Composite score,background score,% cov background,% on target background,' \
                    'true % cov background,Composite score background,Delta composite score vs background,' \
                    'Optimized guide,Optimized guide + G + IGS\n'
        with open(f'{folder_path}/Targeted designs against background above threshold.csv', 'w') as f:
            f.write(first_row)
            for ribo_design in designs:
                f.write(ribo_design.print_attributes())
        return

    elif designs[0].optimized_to_targets:
        if os.path.exists(f'{folder_path}/All targeted designs.csv'):
            os.remove(f'{folder_path}/All targeted designs.csv')
        first_row = 'IGS,Reference index,Number of species targeted,Score,Score type,% cov,% on target,True % cov,' \
                    'Composite score,Optimized guide,Optimized guide + G + IGS\n'
        with open(f'{folder_path}/All targeted designs.csv', 'w') as f:
            f.write(first_row)
            for ribo_design in designs:
                f.write(ribo_design.print_attributes())
        return
    elif designs[0].tested_design:
        if os.path.exists(f'{folder_path}/Sequence tested on background.csv'):
            os.remove(f'{folder_path}/Sequence tested on background.csv')
        first_row = 'IGS,Reference index,Number of tested species targeted,Putative targets,Score type,Design score,' \
                    'Score on targets,% cov background,% on target background,true % cov background,' \
                    'Composite score background,Optimized guide,Optimized guide + G + IGS,Ideal guide,' \
                    'Ideal guide + G + IGS\n'
        with open(f'{folder_path}/Sequence tested on background.csv', 'w') as f:
            f.write(first_row)
            for ribo_design in designs:
                f.write(ribo_design.print_attributes(taxonomy=taxonomy))
    else:
        print('Please optimize to targets or optimize to targets then to background first.')
        return


def make_graphs(control_designs: np.ndarray[RibozymeDesign], universal_designs: list[np.ndarray[RibozymeDesign]],
                selective_designs: list[np.ndarray[RibozymeDesign]], var_regs: list[int], save_fig: bool = False,
                file_loc: str = None, file_type: str = 'png', taxonomy: str = 'Order', data_file: str = ''):
    # Make data into a pandas dataframe
    if not data_file:
        all_data_array = np.hstack([control_designs] + universal_designs + selective_designs)
        index_names = ['control'] * len(control_designs) + \
                      [f'universal_{i}' for i, data in enumerate(universal_designs) for _ in data] + \
                      [f'selective_{i}' for i, data in enumerate(selective_designs) for _ in data]

        all_data_df = pd.DataFrame.from_records([item.to_dict(taxonomy=taxonomy) for item in all_data_array],
                                                index=index_names)
        if os.path.exists(file_loc):
            os.remove(file_loc)
        all_data_df.to_csv(file_loc)
    else:
        all_data_df = pd.read_csv(data_file, index_col=0)

    # Extract top scores for each category
    labels = list(set(all_data_df.index))
    labels.sort()
    control_designs_df = all_data_df.loc['control']
    universal_designs_df = all_data_df.loc[[name for name in labels if name[0] == 'u']]
    selective_designs_df = all_data_df.loc[[name for name in labels if name[0] == 's']]
    background_designs_df = selective_designs_df.append(control_designs_df)

    # get top scores in each: 1 control, one of each group of selective and universal
    try:
        top_score_control = control_designs_df.sort_values('composite_background_score', ascending=False)
    except:
        top_score_control = control_designs_df

    top_scores_universal = universal_designs_df.sort_values('composite_score', ascending=False).groupby(
        universal_designs_df.index).head(2)
    top_scores_selective = selective_designs_df.sort_values('composite_background_score', ascending=False).groupby(
        selective_designs_df.index).head(2)

    top_background_scores = top_scores_selective.append(top_score_control)

    # Set plot parameters
    custom_params = {"axes.spines.right": False, "axes.spines.top": False, 'figure.figsize': (20, 16)}
    sns.set_theme(context='talk', style="ticks", rc=custom_params, palette='viridis')

    # The data we will need is:
    # 1a: average U conservation in background, igs conservation, delta composite score, counts === work on this later
    # 1b: guide score, igs true % coverage === work on this later
    # 2a: average U conservation in background, igs conservation, guide score, ref_idx === work on this later
    # 2b: composite score, ref_idx === work on this later
    # 3: composite scores === work on this later
    # 4: delta score for selective designs, ref_idx === work on this later
    # 5a and 5b: orders, counts for each order

    # Graph 1a, 1b: old plots, showing the distribution of scores for each design on background:
    # guide scores with flexible IGS, scores with rigid IGS

    # fig_1a = plt.subplots(sharex='all', layout='constrained')
    # axes_1a = fpf.generate_axes(fig_1a, '1a', design_dict.keys())
    # fpf.generate_summary_graphs(fig_1a, axes_1a, design_dict)
    #
    # fig_1b = plt.subplots(sharex='all', layout='constrained')
    # axes_1b = fpf.generate_axes(fig_1b, '1b', design_dict.keys())
    # fpf.generate_igs_vs_guide_graph(fig_1b, axes_1b, design_dict)

    # Graph 2: y-axis is the composite score, x-axis is the 16s rRNA gene, plot the universal, control, and each
    # selective for all designs in different panels (same data as above but order along gene)
    for selective_dataset in [name for name in labels if name[0] == 's']:
        fig2a, ax2a = plt.subplots(3, sharex='all')
        fig2a.suptitle(selective_dataset)
        fig_2a_titles = ['Average conservation of catalytic site', '% of reads with IGS', 'Guide score']
        fig_2a_columns = ['u_conservation_background', 'true_%_cov_background', 'background_score']

        for ax, name, col in zip(ax2a, fig_2a_titles, fig_2a_columns):
            fpf.plot_variable_regions(ax, var_regs)
            sns.scatterplot(x='reference_idx', y=col, data=selective_designs_df.loc[selective_dataset], ax=ax,
                            alpha=0.7)
            ax.set_title(name)
        plt.tight_layout()
        plt.show()

    fig2b, ax2b = plt.subplots(len(labels), sharex='all')
    for selective_dataset, ax in zip(labels, ax2b):
        fpf.plot_variable_regions(ax, var_regs)
        ax.set_title(selective_dataset)
        try:
            if selective_dataset[0] == 'u':
                sns.scatterplot(x='reference_idx', y='composite_score', data=all_data_df.loc[selective_dataset], ax=ax,
                                alpha=0.7)
            else:
                sns.scatterplot(x='reference_idx', y='composite_background_score',
                                data=all_data_df.loc[selective_dataset],
                                ax=ax, alpha=0.7)
        except:
            continue
    plt.tight_layout()
    plt.show()

    # Graph 3: violin plot, showing composite score distribution in all designs of a given category
    fig_3_data = pd.concat([universal_designs_df['composite_score'],
                            background_designs_df['composite_background_score'].rename(
                                {'composite_background_score': 'composite_score'})])

    sns.violinplot(x=fig_3_data.index, y=fig_3_data, inner=None, cut=0)
    plt.tight_layout()
    plt.show()

    # Graph 4: y-axis is the delta score, x-axis is the 16s rRNA gene, plot all selective designs in different panels
    fig4, ax4 = plt.subplots(3, sharex='all')
    # fpf.generate_axes(fig3, graph='4', titles=list(selective_designs_df.index))

    for ax, name in zip(ax4, [name for name in labels if name[0] == 's']):
        fpf.plot_variable_regions(ax, var_regs)
        sns.scatterplot(x='reference_idx', y='delta_vs_background', data=selective_designs_df.loc[name], ax=ax,
                        alpha=0.7)
        ax.set_title(name)
    plt.tight_layout()
    plt.show()

    # Graph 5a: bar graph, y-axis is fraction of total reads, x-axis is each condition: control (original design),
    # universal, and three selective designs one for a different order in wastewater. One design each, with the best
    # composite score. Maybe pick the best performing one with flexible and rigid IGS

    # Extract possible taxonomy labels
    all_targets = set()
    universal_targets = [targets.replace('[', '').replace(']]', '').split(']; ') for targets in
                         top_scores_universal['targets']]
    background_targets = [targets.replace('[', '').replace(']]', '').split(']; ') for targets in
                          top_background_scores['background_targets']]

    for items in (universal_targets + background_targets):
        for item in items:
            try:
                name, num = item.split('; ')
            except:
                continue
            all_targets.add(name)

    fig_5_data = pd.DataFrame(columns=list(all_targets))
    # universal
    for i, design_targets in enumerate(top_scores_universal['targets']):
        dict = {}
        all_targets = design_targets.replace('[', '').replace(']]', '').split(']; ')
        for item in all_targets:
            try:
                name, num = item.split('; ')
                dict[name] = int(num)
            except:
                continue
        fig_5_data.loc[top_scores_universal.index[i]] = dict
    for i, design_targets in enumerate(top_background_scores['background_targets']):
        dict = {}
        all_targets = design_targets.replace('[', '').replace(']]', '').split(']; ')
        for item in all_targets:
            try:
                name, num = item.split('; ')
                dict[name] = int(num)
            except:
                continue
        fig_5_data.loc[top_background_scores.index[i]] = dict
    filter_nans = fig_5_data.groupby(fig_5_data.index).sum()

    fig5, ax5 = plt.subplots(ncols=2, nrows=1, sharey='all', layout='constrained')
    # Figure 5a:
    # Get 100%
    totals = filter_nans.sum(axis=1)
    # Calculate fractions
    fractions = filter_nans.div(totals, axis=0)
    fractions.plot.barh(stacked=True, ax=ax5[0], legend=0)
    ax5[0].set_xlabel('Fraction')
    # Figure 5b
    filter_nans.plot.barh(stacked=True, ax=ax5[1], legend=0)
    ax5[1].set_xlabel('Counts')
    fig5.suptitle('Orders targeted of dataset')
    leg_handles, leg_labels = ax5[0].get_legend_handles_labels()
    fig5.legend(leg_handles, leg_labels, ncols=math.ceil(len(leg_labels) / 14), loc='lower center', fancybox=True,
                fontsize='x-small')
    # show the graph
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.27)
    plt.show()

    return


def give_taxonomy(silva_name: str, level: str):
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
            here_is_taxonomy = all_levels[level_names[level]]
        except IndexError:
            here_is_taxonomy = ''
        out_names.append(here_is_taxonomy)
    if out_names:
        counts_and_names = np.asarray(np.unique(out_names, return_counts=True)).T
        return counts_and_names.tolist()
    else:
        return None
