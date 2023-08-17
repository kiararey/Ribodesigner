import glob
import os
import random
import time
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy.random import default_rng
from Bio.Seq import Seq
import Bio.motifs
from collections import defaultdict
import subprocess
from Bio.Align import AlignInfo, MultipleSeqAlignment, PairwiseAligner
from Bio.Align.Applications import MuscleCommandline
from collections import Counter
from Bio import SeqIO, AlignIO
from math import exp, log
from multiprocessing import Pool
import seaborn as sns
import matplotlib.pyplot as plt
import scipy


# Expanded version of SilvaSequence
class TargetSeq:
    # Important: TargetSeq id must be formatted like a Silva sequence if we want the give_taxonomy function to work!
    def __init__(self, id: str = '', seq: Seq = ''):
        self.id = id
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
                 perc_on_target_background_attr: float = None, true_perc_cov_background_attr: float = None):
        self.id = id_attr
        self.igs = id_attr[:5]
        self.ref_idx = int(id_attr[5:])
        self.guide = guide_attr
        self.guides_to_use = guides_to_use_attr
        self.targets = targets_attr
        self.number_of_targets = len(targets_attr)
        self.score = score_attr
        self.score_type = score_type_attr

        # calculate percent coverage, percent on target, or true percent coverage if needed and possible.
        self.perc_cov = perc_cov_attr
        self.perc_on_target = perc_on_target_attr
        self.true_perc_cov = true_perc_cov_attr
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

        # For initialized targeted ribozyme designs only:
        self.background_score = background_score_attr
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

    def __str__(self):
        # ID: igs + ref_idx, then also display the guide sequence
        return f'ID:{self.igs}{self.ref_idx}, Guide:{self.guide}'

    def __repr__(self):
        return f'{type(self).__name__}(IGS={self.igs}, ref_idx={self.ref_idx}, guide_sequence={self.guide}, ' \
               f'score={self.score}, score_type={self.score_type}, optimized_to_targets={self.optimized_to_targets}, ' \
               f'optimized_to_background={self.optimized_to_background}'

    def calc_percent_coverages(self, perc_cov_attr: float = None, true_perc_cov_attr: float = None,
                               perc_on_target_attr: float = None):
        # calculate percent coverage, percent on target, or true percent coverage if needed and possible.
        if not perc_cov_attr and true_perc_cov_attr and perc_on_target_attr:
            self.perc_cov = true_perc_cov_attr / perc_on_target_attr

        if not perc_on_target_attr and true_perc_cov_attr and perc_cov_attr:
            self.perc_on_target = true_perc_cov_attr / perc_cov_attr

        if not true_perc_cov_attr and perc_on_target_attr and perc_cov_attr:
            self.true_perc_cov = perc_on_target_attr * perc_cov_attr

    def calc_background_percent_coverages(self, perc_cov_background_attr: float = None,
                                          true_perc_cov_background_attr: float = None,
                                          perc_on_target_background_attr: float = None):
        # calculate percent coverage, percent on target, or true percent coverage if needed and possible.
        if not perc_cov_background_attr and true_perc_cov_background_attr and perc_on_target_background_attr:
            self.perc_cov_background = true_perc_cov_background_attr / perc_on_target_background_attr

        if not perc_on_target_background_attr and true_perc_cov_background_attr and perc_cov_background_attr:
            self.perc_on_target_background = true_perc_cov_background_attr / perc_cov_background_attr

        if not true_perc_cov_background_attr and perc_on_target_background_attr and perc_cov_background_attr:
            self.true_perc_cov_background = perc_on_target_background_attr * perc_cov_background_attr


    def update_after_optimizing(self, score_attr, guide_attr, score_type_attr):
        self.score = score_attr
        self.guide = guide_attr
        # self.guides_to_use = None
        self.score_type = score_type_attr
        self.optimized_to_targets = True
        self.composite_score = self.true_perc_cov * score_attr

    def update_to_background(self, updated_guide_attr: Seq, background_score_attr: float,
                             perc_cov_background_attr: float, perc_on_target_background_attr: float,
                             true_perc_cov_background_attr: float):
        self.guide = updated_guide_attr
        self.optimized_to_background = True
        self.background_score = background_score_attr
        self.calc_background_percent_coverages(perc_cov_background_attr=perc_cov_background_attr,
                                               perc_on_target_background_attr=perc_on_target_background_attr,
                                               true_perc_cov_background_attr=true_perc_cov_background_attr)
        self.composite_background_score = self.true_perc_cov_background * background_score_attr
        self.delta_vs_background = self.composite_score - self.composite_background_score

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

    time1 = time.perf_counter()
    print(f'Now finding catalytic sites and aligning to reference {ref_name_and_seq[0].replace("_", " ")}')
    fn_align_to_ref = np.vectorize(align_to_ref, otypes=[TargetSeq],
                      excluded=['ref_name_and_seq', 'igs_length', 'guide_length', 'min_length'])
    fn_align_to_ref(target_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length,
                    guide_length=guide_length, min_length=min_length)
    time2 = time.perf_counter()
    round_convert_time(start=time1, end=time2, round_to=4,
                       task_timed=f'finding catalytic sites and indexing target sequences to reference')

    time1 = time.perf_counter()
    if optimize_seq:
        # first, filter through possible designs and get the ones that meet our threshold
        to_optimize = filter_igs_candidates(aligned_targets=target_names_and_seqs, min_true_cov=min_true_cov)
        time2 = time.perf_counter()
        round_convert_time(start=time1, end=time2, round_to=4, task_timed='finding putative ribozyme sites')

        # Now, optimize all of the possible guide sequences for each IGS:
        print('Now optimizing guide sequences...')
        time1 = time.perf_counter()
        fn_msa = np.vectorize(optimize_designs, otypes=[RibozymeDesign],
                              excluded=['score_type', 'thresh', 'msa_fast', 'gaps_allowed'])
        fn_msa(to_optimize, score_type=score_type, thresh=identity_thresh, msa_fast=msa_fast, gaps_allowed=gaps_allowed)
        time2 = time.perf_counter()
        round_convert_time(start=time1, end=time2, round_to=4, task_timed='generating optimized designs')

        if targeted:
            background_names_and_seqs = read_silva_fasta(in_file=background_sequences_folder)

            print(f'Found {background_names_and_seqs.size} total background sequences to analyze.')

            if percent_of_background_seqs_used < 1:
                background_names_and_seqs = np.random.choice(
                    background_names_and_seqs, size=round(background_names_and_seqs.size * percent_of_target_seqs_used),
                    replace=False)
                print(f'Randomly sampling {target_names_and_seqs.size} sequences to analyze.\n')

            # first find the IGSes and locations of the background sequences. We do not want to hit these.
            # Also align these to the reference sequence
            fn_align_to_ref(background_names_and_seqs, ref_name_and_seq=ref_name_and_seq, igs_length=igs_length,
                            guide_length=guide_length, min_length=min_length)
            time2 = time.perf_counter()
            round_convert_time(start=time1, end=time2, round_to=4,
                               task_timed=f'finding catalytic sites and indexing background sequences to reference')

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
                seq_id = title.split(' ')[0]
                putative_sequence = TargetSeq(seq_id, Seq(seq).upper().back_transcribe())
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
            igs_over_min_true_cov[igs_id][1].add(aligned_targets[target_num].reset_cat_sites_and_seq())

    print(f'{len(igs_over_min_true_cov)} putative designs found.')
    # Now make an array of all of the putative designs for later use.
    to_optimize = np.array([RibozymeDesign(id_attr=igs_id, guides_to_use_attr=item[0], targets_attr=item[1],
                                           true_perc_cov_attr=item[2],
                                           perc_cov_attr=len(item[1]) / aligned_targets.size)
                            for igs_id, item in igs_over_min_true_cov.items()])
    return to_optimize


def optimize_designs(to_optimize: RibozymeDesign, score_type: str, thresh: float = 0.7, msa_fast: bool = True,
                     gaps_allowed: bool = True):
    """
    Takes in a RibozymeDesign with guide list assigned and uses MUSCLE msa to optimize the guides.

    """
    # First create a dummy file with all the guides we want to optimize
    with open(f'to_align_{to_optimize.id}.fasta', 'w') as f:
        for line in range(to_optimize.number_of_targets):
            f.write('>seq' + str(line) + '\n' + str(to_optimize.guides_to_use[line]) + '\n')

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

        if score_type == 'naive':
            # Naively tell me what the percent identity is: 1- ambiguity codons/ length
            score = get_naive_score(opti_seq)

        elif score_type == 'weighted':
            score = get_weighted_score(opti_seq)

        elif score_type == 'directional':
            score = get_directional_score(opti_seq)

    else:
        if gaps_allowed:
            left_seq = summary_align.gap_consensus(threshold=thresh, ambiguous='N')
        else:
            left_seq = summary_align.dumb_consensus(threshold=thresh, ambiguous='N')
        score, opti_seq = get_quantitative_score(left_seq, alignments, count_gaps=gaps_allowed)

    # Now update our initial design
    to_optimize.update_after_optimizing(score_attr=score, guide_attr=opti_seq, score_type_attr=score_type)


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
         'H': 1/3, 'B': 1/3, 'D': 1/3, 'V': 1/3, 'N': 0, '-': 0}
    n = len(opti_seq)

    prelim_score = 0
    for i, x in enumerate(str(opti_seq)):
        prelim_score += a[x]
    score = prelim_score / n

    return score


def get_directional_score(opti_seq):
    # prepare scoring matrix a(x)
    a = {'A': 1, 'T': 1, 'G': 1, 'C': 1, 'U': 1, 'R': 0.5, 'Y': 0.5, 'M': 0.5, 'K': 0.5, 'S': 0.5, 'W': 0.5,
         'H': 1/3, 'B': 1/3, 'D': 1/3, 'V': 1/3, 'N': 0, '-': 0}
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


def ribo_checker(designs, aligned_target_sequences, ref_seq_len, identity_thresh: float, guide_length: int,
                 score_type: str = 'quantitative', gaps_allowed: bool = True, msa_fast: bool = False,
                 flexible_igs: bool = False, do_not_target_background: bool = True, n_limit: int = 0):
    """
    Checks generated designs against a set of sequences to either find how well they align to a given dataset.
    Will return a set of sequences that match for each design as well as a score showing how well they matched.
    :param designs:
    :param aligned_target_sequences:
    :param ref_seq_len:
    :param identity_thresh:
    :param guide_length:
    :param do_not_target_background: boolean set to True if you want to further optimize the designs by reducing ambiguity
    :param score_type:
    :param gaps_allowed:
    :param flexible_igs:
    :param msa_fast:
    :param n_limit: what percentage of Ns are acceptable?
    :return:
    """
    # extract all design IGSes and guides.
    # This will be a list: [igs, guide, ref_idx (0 base indexing), guide_id(igs+ref_idx)]
    igs_and_guides_designs = [(str(design[0]), design[8], design[1] - 1, str(design[0]) + str(design[1]))
                              for design in designs]
    igs_and_guides_initial_comp_scores = {str(design[0]) + str(design[1]): design[6] for design in designs}

    # Is there a U at each position for each design at each target sequence?
    designed_idxs = np.array(list({design[2] for design in igs_and_guides_designs}))  # base 0 indexing
    uracil_sites_targets = np.zeros((len(aligned_target_sequences), ref_seq_len))
    uracil_sites_targets[:, designed_idxs] = 1
    for i, target_data in enumerate(aligned_target_sequences):  # base 1 indexing
        temp_uracil_indexes = np.array([ref_idx - 1 for _, ref_idx in target_data[2][2]])  # convert to base 0 indexing
        # everything that is a 2 is a U that appears in both the designs and the target sequence,
        # a 1 is a U in the design
        uracil_sites_targets[i, temp_uracil_indexes] = uracil_sites_targets[i, temp_uracil_indexes] * 2  # base 0 index

    average_conserved = []
    num_of_targets = len(aligned_target_sequences)
    for i in designed_idxs:
        u_values, u_counts = np.unique(uracil_sites_targets[:, i], return_counts=True)
        try:
            percentage_of_cat_sites_conserved = u_counts[np.where(u_values == 2)[0][0]]/num_of_targets
        except:
            percentage_of_cat_sites_conserved = 0
        average_conserved.append(percentage_of_cat_sites_conserved)

    # # for histogram
    # u_values, u_counts = np.unique(uracil_sites_targets, return_counts=True)  # base 0 indexing
    # is_u_conserved = [[False, True], u_counts[1:]]
    # print(
    #     f'There are {is_u_conserved[1][1]} conserved Us and {is_u_conserved[1][0]} not conserved.')

    # What is the conservation of IGSes at a particular position for each design?
    # Keep only the sites that are conserved for later.
    igs_sites_targets = np.array(np.clip(uracil_sites_targets, 0, 1), copy=True, dtype=str)  # base 0 indexing
    conserved_u_sites = np.argwhere(uracil_sites_targets == 2)

    for target_number, ref_idx in conserved_u_sites:  # base 0 indexing
        # extract IGS there:
        index_of_target_data = np.argwhere(np.array(aligned_target_sequences[target_number]
                                                    [2][2])[:, 1] == ref_idx + 1)[0, 0]  # convert to base 1 indexing
        igs_sites_targets[target_number, ref_idx] = str(aligned_target_sequences[target_number][2][0]
                                                        [index_of_target_data])

    # Now find how many conserved IGSes there are!
    to_optimize = {}
    no_splices_in_background = []
    opti_seqs = []

    conserved_igs_true_perc_coverage = {}

    for igs, designed_guide, ref_idx, guide_id in igs_and_guides_designs:  # base 0 indexing
        check = np.argwhere(igs_sites_targets == igs)
        orgs_with_igs_on_target = check[np.argwhere(check[:, 1] == ref_idx)][:, 0, 0]
        on_target_count = len(orgs_with_igs_on_target)
        if on_target_count > 0:
            true_perc_cov = on_target_count / len(aligned_target_sequences)
            orgs_with_igs, counts = np.unique(check[:, 0], return_counts=True)
            perc_cov = len(orgs_with_igs) / len(aligned_target_sequences)
            perc_on_target = true_perc_cov / perc_cov
            conserved_igs_true_perc_coverage[guide_id] = true_perc_cov
            if not flexible_igs:
                all_indexes_of_target_data = [np.argwhere(np.array(aligned_target_sequences[target][2][2])[:, 1]
                                                          == ref_idx + 1)[0, 0] for target in
                                              orgs_with_igs_on_target]  # convert to base 1 indexing
                guides_to_optimize = [str(aligned_target_sequences[target][2][1][index_of_target_data]) for
                                      target, index_of_target_data in
                                      zip(orgs_with_igs_on_target, all_indexes_of_target_data)]
                # names_and_stuff = [(aligned_target_sequences[target][0], ig_idx, occurs_in_target)
                #                    for target, ig_idx, occurs_in_target in
                #                    zip(orgs_with_igs_on_target, all_indexes_of_target_data, counts)]
            else:
                orgs_with_u_on_target = conserved_u_sites[np.where(conserved_u_sites[:, 1] == ref_idx)][:, 0]
                all_indexes_of_target_data = [np.argwhere(np.array(aligned_target_sequences[target][2][2])[:, 1]
                                                          == ref_idx + 1)[0, 0] for target in
                                              orgs_with_u_on_target]  # convert to base 1 indexing
                guides_to_optimize = [str(aligned_target_sequences[target][2][1][index_of_target_data]) for
                                      target, index_of_target_data in
                                      zip(orgs_with_u_on_target, all_indexes_of_target_data)]
                # names_and_stuff = [(aligned_target_sequences[target][0], ig_idx, occurs_in_target)
                #                    for target, ig_idx, occurs_in_target in
                #                    zip(orgs_with_u_on_target, all_indexes_of_target_data, counts)]
            names_and_stuff = len(guides_to_optimize)
            if len(guides_to_optimize) > 1:
                to_optimize[guide_id] = [igs, ref_idx + 1, perc_cov, perc_on_target, true_perc_cov, names_and_stuff,
                                         guides_to_optimize, designed_guide]
            else:
                opti_seqs.append([igs, ref_idx + 1, 1, perc_cov, perc_on_target, true_perc_cov, 1 * true_perc_cov,
                                  names_and_stuff, Seq(guides_to_optimize[0]), designed_guide])
        else:
            conserved_igs_true_perc_coverage[guide_id] = 0
            # remove designs with too many Ns
            if designed_guide.count('N')/guide_length <= n_limit:
                # Follow the same naming as the outputs of replace_ambiguity. Pairwise score is set at 0 based on the
                # assumption that if there is no catalytic U there is no splicing activity
                initial_score = igs_and_guides_initial_comp_scores[guide_id]
                no_splices_in_background.append((designed_guide, initial_score, 0, 0, guide_id, 0))

    # Do an MSA
    opti_seqs.extend(optimize_sequences(to_optimize, identity_thresh, guide_length, '', [], gaps_allowed=gaps_allowed,
                                        score_type=score_type, msa_fast=msa_fast, for_comparison=True))

    # Now reduce ambiguity and score
    print('\nNow scoring ribozyme designs against background sequences...')
    start = time.perf_counter()
    in_data = [(seqs[-1], seqs[-2], do_not_target_background, score_type, key, n_limit, seqs[5]) for key, seqs in
               zip(to_optimize, opti_seqs)]

    replace_ambiguity(in_data[0][0], in_data[0][1], in_data[0][2], in_data[0][3], in_data[0][4], in_data[0][5],
                      in_data[0][6])


    with Pool() as pool:
        less_ambiguous_designs = pool.starmap(replace_ambiguity, in_data)

    less_ambiguous_designs.extend(no_splices_in_background)
    pairwise_composite_scores = {}
    background_guide_scores = {}
    igs_and_guides_comp_scores = {}

    # for new_seq, new_seq_score, background_seq_score, pairwise_score, key, true_perc_cov in less_ambiguous_designs:
    for new_seq, new_seq_score, background_seq_score, pairwise_score, key, true_perc_cov in filter(None, less_ambiguous_designs):

        # Update our scores with the new less ambiguous designs
        igs_and_guides_comp_scores[key] = [new_seq_score * true_perc_cov, new_seq]
        pairwise_composite_scores[key] = pairwise_score * true_perc_cov
        background_guide_scores[key] = background_seq_score

    delta_composite_scores = {}

    for key, good_data in igs_and_guides_comp_scores.items():
        bad_score = pairwise_composite_scores[key]
        delta_composite_scores[key] = good_data[0] - bad_score
    round_convert_time(start=start, end=time.perf_counter(), round_to=4,
                       task_timed='scoring against background sequences')

    return average_conserved, conserved_igs_true_perc_coverage, delta_composite_scores, background_guide_scores, \
           igs_and_guides_comp_scores