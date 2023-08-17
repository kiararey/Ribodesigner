import glob
import os
import random
import warnings
from Bio.Seq import Seq
import Bio.motifs
import re
import pandas as pd
import numpy as np
from collections import defaultdict
import subprocess
from Bio.Align import AlignInfo, MultipleSeqAlignment, PairwiseAligner
from Bio.Align.Applications import MuscleCommandline
from collections import Counter
from math import exp, log
from multiprocessing import Pool
import time
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from numpy.random import default_rng
from Bio import SeqIO, AlignIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import repeat


# Expanded version of SilvaSequence
class TargetSeq:
    # Important: TargetSeq id must be formatted like a Silva sequence if we want the give_taxonomy function to work!
    def __init__(self, id: str = '', taxonomy: str = '', seq: Seq = ''):
        self.id = id
        self.taxonomy = taxonomy
        self.seq = seq
        self.cat_sites = []

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
        all_levels = self.id.split(';')
        try:
            here_is_taxonomy = all_levels[level_names[level]]
        except IndexError:
            here_is_taxonomy = ''
        return here_is_taxonomy


class RibozymeDesign:
    def __init__(self, id: str, guides_to_use: list[Seq], targets: set, guide: Seq = '', score: int = None,
                 score_type: str = '', perc_cov: float = None, perc_on_target: float = None,
                 true_perc_cov: float = None, composite_background_score: float = None):
        self.id = id
        self.ref_idx = int(id[5:])
        self.guide = guide
        if not guide:
            self.guides_to_use = guides_to_use
        else:
            self.guides_to_use = None
        self.targets = targets
        self.igs = id[:5]
        self.score = score
        self.score_type = score_type
        self.perc_cov = perc_cov
        if not perc_on_target and perc_cov and true_perc_cov:
            self.perc_on_target = true_perc_cov / perc_cov
        else:
            self.perc_on_target = perc_on_target
        self.true_perc_cov = true_perc_cov
        if score:
            self.composite_score = true_perc_cov * score
        else:
            self.composite_score = None
        if guide and score and score_type:
            self.optimized_to_targets = True
        else:
            self.optimized_to_targets = False
        # For targeted ribozyme designs only:
        self.composite_background_score = composite_background_score
        if self.composite_score and self.composite_background_score:
            self.delta_vs_background = self.composite_score - composite_background_score
        else:
            self.delta_vs_background = None
        self.number_of_targets = len(targets)
        if composite_background_score:
            self.optimized_to_background = True
        else:
            self.optimized_to_background = False

    def __str__(self):
        # ID: igs + ref_idx, then also display the guide sequence
        return f'ID:{self.igs}{self.ref_idx}, Guide:{self.guide}'

    def __repr__(self):
        return f'{type(self).__name__}(IGS={self.igs}, ref_idx={self.ref_idx}, guide_sequence={self.guide}, ' \
               f'score={self.score}, score_type={self.score_type}, optimized_to_targets={self.optimized_to_targets}, ' \
               f'optimized_to_background={self.optimized_to_background}'



def ribodesigner(target_sequences_folder: str, barcode_seq_file: str, ribobody_file: str, igs_length: int = 5,
                 guide_length: int = 50, min_length: int = 35, ref_sequence_file=None, targeted: bool = False,
                 background_sequences_folder: str = '', min_delta: float = 0, optimize_seq: bool = True,
                 min_true_cov: float = 0.7, identity_thresh: float = 0.7, fileout: bool = False,
                 folder_to_save: str = '', score_type: str = 'quantitative', msa_fast: bool = False,
                 keep_single_targets: bool = False, gaps_allowed: bool = True, percent_of_target_seqs_used: float = 1.0,
                 percent_of_background_seqs_used: float = 1, seed_target: int = 1, seed_background: float = 1,
                 n_limit: int = 0):
    """Generates ribozyme designs to target a set of sequences.
    :param generate_summary: this will generate a summary graph of the target sequences showing the IGS true percent
    coverage distribution vs. the guide score distribution. Keep in mind that setting this to True will SIGNIFICANTLY
    increase computational load, because in essence the program would have to calculate guide scores for any possible
    location regardless of min_true_cov. Only recommended for small datasets, very fast computers, or computers with a
    lot of memory owned by people with a lot of patience.
    :param seed: seed to use for random sampling for reproducibility
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
    :param seed_target:
    :param seed_background:
    :param n_limit: What percentage of Ns are acceptable in the final design? This is applied when doing targeted
    designs only for the time being.
    """

    start = time.perf_counter()

    # Make the Ribozyme sequence by combining the main body and the barcode
    barcode_seq = back_transcribe_seq_file(barcode_seq_file)
    ribobody = back_transcribe_seq_file(ribobody_file)
    ribo_seq = ribobody + barcode_seq

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
        ref_name_and_seq = random.choice(target_names_and_seqs)
    else:
        ref_name_and_seq = read_fasta(ref_sequence_file)[0]

    print(f'Found {target_names_and_seqs.size} total target sequences to analyze.')

    if percent_of_target_seqs_used < 1:
        target_names_and_seqs = np.random.choice(target_names_and_seqs,
                                                 size=round(target_names_and_seqs.size * percent_of_target_seqs_used),
                                                 replace=False)
        print(f'Randomly sampling {target_names_and_seqs.size} sequences to analyze.\n')

    if keep_single_targets and min_true_cov > 1 / len(target_names_and_seqs):
        # If we only want to get more than a certain percentage, override keep_single_targets
        keep_single_targets = False

    time1 = time.perf_counter()
    fn = np.vectorize(align_to_ref, otypes=[TargetSeq],
                      excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])
    fn(target_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length, guide_length=guide_length,
       min_length=min_length)
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'finding catalytic sites and indexing sequences to reference '
                                  f'{ref_name_and_seq[0].replace("_", " ")}')

    time1 = time.perf_counter()
    if optimize_seq:
        # first, filter through possible designs and get the ones that meet our threshold
        to_optimize = filter_igs_candidates(target_names_and_seqs, min_true_cov)
        time2 = time.perf_counter()
        round_convert_time(start=time1, end=time2, round_to=4, task_timed='prepping sequences for optimization')

        # Now, optimize all of the possible guide sequences for each IGS:
        time1 = time.perf_counter()

    #     time1 = time.perf_counter()
    #     opti_seqs = optimize_sequences(to_optimize, identity_thresh, guide_length, ribo_seq, to_keep_single_targets,
    #                                    fileout=fileout, file=folder_to_save, score_type=score_type, msa_fast=msa_fast,
    #                                    gaps_allowed=gaps_allowed, min_true_cov=min_true_cov)
    #     time2 = time.perf_counter()
    #     round_convert_time(start=time1, end=time2, round_to=4, task_timed='generating optimized designs')
    #
    #     if targeted:
    #         background_names_and_seqs = read_fasta(background_sequences_folder)
    #         print(f'Found {len(background_names_and_seqs)} total background sequences to analyze.')
    #
    #         if percent_of_background_seqs_used < 1:
    #             array = np.array(background_names_and_seqs, dtype=tuple)
    #             rng = default_rng(seed=seed_background)
    #             background_names_and_seqs = rng.choice(
    #                 array, size=round(len(background_names_and_seqs) * percent_of_background_seqs_used), replace=False)
    #             print(f'Randomly sampling {len(background_names_and_seqs)} background sequences to analyze.\n')
    #         # first find the IGSes and locations of the background sequences. We do not want to hit these.
    #         background_sequences_data = find_cat_sites(background_names_and_seqs, igs_length,
    #                                                    guide_length, min_length)
    #
    #         time1 = time.perf_counter()
    #         aligned_background_sequences = align_to_ref(background_sequences_data, ref_name_and_seq)
    #         time2 = time.perf_counter()
    #         round_convert_time(start=time1, end=time2, round_to=4, task_timed='indexing background sequences')
    #
    #         print('Now applying designed ribozymes with background sequences and getting statistics...')
    #         average_conserved, conserved_igs_true_perc_coverage, delta_composite_scores, background_guide_scores, \
    #         igs_and_guides_comp_scores = \
    #             ribo_checker(opti_seqs, aligned_background_sequences, len(ref_name_and_seq[1]),
    #                          identity_thresh=identity_thresh, guide_length=guide_length, score_type=score_type,
    #                          gaps_allowed=gaps_allowed, msa_fast=msa_fast, flexible_igs=True, n_limit=n_limit)
    #
    #         opti_target_seqs = compare_targeted_sequences(opti_seqs, average_conserved, conserved_igs_true_perc_coverage,
    #                                                       background_guide_scores, delta_composite_scores,
    #                                                       igs_and_guides_comp_scores, min_delta=min_delta,
    #                                                       file_out=fileout, file=folder_to_save)
    #         time2 = time.perf_counter()
    #         round_convert_time(start=time1, end=time2, round_to=4, task_timed='comparing designs against background '
    #                                                                           'sequences')
    #
    #         end = time.perf_counter()
    #         round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    #
    #         return opti_target_seqs
    #
    #     end = time.perf_counter()
    #     round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    #     print('########################################################\n')
    #     return opti_seqs
    #
    # else:
    #     ranked_sorted_IGS = find_repeat_targets(aligned_seqs, ribo_seq, fileout=fileout, file=folder_to_save)
    #     time2 = time.perf_counter()
    #     print('All guide sequences generated.')
    #     round_convert_time(start=time1, end=time2, round_to=4, task_timed='generating designs')
    #     end = time.perf_counter()
    #     round_convert_time(start=start, end=end, round_to=4, task_timed='overall')
    #     print('########################################################\n')
    #     return ranked_sorted_IGS


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
                id = title.split(' ')[0]
                taxonomy = title[len(id) + 1].split(';')
                putative_sequence = TargetSeq(id, taxonomy, Seq(seq).upper().back_transcribe())
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


def filter_igs_candidates(aligned_targets: np.ndarray[TargetSeq], min_true_cov: float = 0) -> None:
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

    print('Finding repeat IGSes...')

    # Extract all the IGS id numbers - that's the IGS sequence and the reference index number
    all_igs_ids = [f'{item[2]}{item[1]}' for seq in aligned_targets for item in seq.cat_sites]
    # Measure true percent coverage of these and keep IGSes that are at least the minimum true percent coverage needed
    igs_ids_counted, igs_counts = np.unique(all_igs_ids, return_counts=True)

    if [igs for igs, counts in zip(igs_ids_counted, igs_counts) if counts > aligned_targets.size]:
        print('WARNING: Multiple IGSes per sequence at the same position detected. True percent coverage ')
        return None

    else:
        print('Test passed.')

        # this gives us a dictionary where the ID is matched to the true percent coverage
    igs_over_min_true_cov = {igs: [[], set(), counts / aligned_targets.size] for igs, counts
                             in zip(igs_ids_counted, igs_counts)
                             if counts / aligned_targets.size >= min_true_cov}

    # Here each item in the list is an IGS id (IGS + reference index), a guide, and the location of the target sequence
    # in our initial aligned_targets array
    igs_subsets = [(f'{item[2]}{item[1]}', item[3], i) for i, seq in enumerate(aligned_targets) for item in
                   seq.cat_sites]

    # get all the guides that have good IGSes
    for igs_id, guide, target_num in igs_subsets:
        if igs_id in igs_over_min_true_cov:
            # Add guide to list
            igs_over_min_true_cov[igs_id][0].append(guide)
            # Add organism number to set. We can calculate the percent coverage with this number later on.
            igs_over_min_true_cov[igs_id][1].add(aligned_targets[target_num].id)

    print(f'{len(igs_over_min_true_cov)} repeat IGSes found.')
    # Now make an array of all of the putative designs for later use.
    to_optimize = np.array([RibozymeDesign(id=igs_id, guides_to_use=item[0], targets=item[1], true_perc_cov=item[2],
                                           perc_cov=len(item[1])/ aligned_targets.size)
                            for igs_id, item in igs_over_min_true_cov.items()])
    return to_optimize

    def optimize_designs(to_optimize: np.ndarray[RibozymeDesign], score_type: str, msa_fast: bool = True)



