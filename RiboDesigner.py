"""
RiboDesigner is a Python program that takes in processed FASTA files containing sequences that we want to target with
RAM, and optionally sequences we want to avoid targeting, and generates Ribozyme designs to target a maximum amount
of targets with as few designs as possible.

Planned features:
- Return design pool: returns several designs to hit everything in a community
- Refine avoid function: make RiboDesigner avoid specific sequences in a pool
- Multiprocessing: cleverly utilize computer memory and cores to process more data, faster
"""

import glob
import os
import random
import warnings
from Bio.Seq import Seq
import re
import pandas as pd
from collections import defaultdict
import subprocess
from Bio.Align import AlignInfo
from math import exp, log
from multiprocessing import Pool
import time
with warnings.catch_warnings():
    # I know pairwise2 is being deprecated but I need to get this to work properly before attempting an update that
    # may break the program entirely.
    warnings.simplefilter('ignore')
    from Bio import SeqIO, pairwise2, AlignIO


def RiboDesigner(target_sequences_folder: str, barcode_seq_file: str, ribobody_file: str, igs_length: int = 5,
                 guide_length: int = 50, min_length: int = 35, ref_sequence_file=None, targeted: bool = False,
                 background_sequences_folder: str = '', optimize_seq: bool = True, min_true_cov: float = 0.7,
                 identity_thresh: float = 0.7, fileout: bool = False, folder_to_save: str = '',
                 score_type: str = 'quantitative', msa_fast: bool = False, keep_single_targets: bool = False):
    """Generates ribozyme designs to target a set of sequences.

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
    :param identity_thresh:
    :param fileout: whether we want a csv file output or not. Default is False.
    :param folder_to_save: the path where the folder we will save our outputs in if fileout = True
    :param score_type:
    :param msa_fast: whether to use super5 MUSCLE MSA or just regular MUSCLE MSA. Recommended for large datasets (over
    300 sequences) for faster data processing.
    :param keep_single_targets:
    """

    start = time.perf_counter()
    # Make the Ribozyme sequence by combining the main body and the barcode
    barcode_seq = transcribe_seq_file(barcode_seq_file)
    ribobody = transcribe_seq_file(ribobody_file)
    ribo_seq = ribobody + barcode_seq

    target_names_and_seqs = read_fasta(target_sequences_folder)

    try:
        target_names_and_seqs[0]
    # if we hit an index error we've run out of sequence and
    # should not add new residues
    except IndexError:
        print(f'No sequences found in {target_sequences_folder}. Please make sure your files are not empty!\n')
        return None

    if not msa_fast and len(target_names_and_seqs) > 300:
        print(f'Consider setting msa_fast=True for large datasets over 300 sequences for faster data processing!!\n')

    if not ref_sequence_file:
        # If we do not have a reference sequence, just choose one randomly
        print('No reference sequence provided. Picking a random sequence as reference...')
        ref_name_and_seq = random.choice(target_names_and_seqs)
    else:
        ref_name_and_seq = read_fasta(ref_sequence_file)[0]

    print(f'Found {len(target_names_and_seqs)} total sequences to analyze.')

    # find all catalytic U sites
    # Remember: data has one tuple per target sequence, the third entry is a tuple for each catalytic U site
    time1 = time.perf_counter()
    data = find_cat_sites(target_names_and_seqs, igs_length, guide_length, min_length)
    time2 = time.perf_counter()
    print(f'Time taken: {time2 - time1}s\n')

    # Now align sequences to reference sequences and get the conversion dictionaries for each
    time1 = time.perf_counter()
    new_data = align_to_ref(data, ref_name_and_seq)
    time2 = time.perf_counter()
    print(f'Time taken: {time2 - time1}s\n')
    time1 = time.perf_counter()

    if targeted:
        print('Now finding targeted designs...')
        time1 = time.perf_counter()
        # get all background sequences
        background_names_and_seqs = read_fasta(background_sequences_folder)

        # first find the IGSes and locations of the background sequences. We do not want to hit these.
        background_sequences_data = find_cat_sites(background_names_and_seqs, igs_length,
                                                   guide_length, min_length)

        # Then, use these data to only keep the designs that are least likely to work in background seqs

        out_data = filter_for_targets(new_data=new_data, background_sequences_data=background_sequences_data,
                                      guide_length=guide_length, ribo_seq=ribo_seq, ref_name_and_seq=ref_name_and_seq,
                                      optimize=optimize_seq, min_true_cov=min_true_cov, identity_thresh=identity_thresh,
                                      fileout=fileout, folder_to_save=folder_to_save, score_type=score_type,
                                      msa_fast=msa_fast, keep_single_targets=keep_single_targets)

        time2 = time.perf_counter()
        print(f'Time taken: {time2 - time1}s\n')
        return out_data

    # Now, we can optimize each sequence
    if optimize_seq:

        to_optimize, to_keep_single_targets = prep_for_optimizing(new_data, min_true_cov=min_true_cov,
                                                                  accept_single_targets=keep_single_targets)
        time2 = time.perf_counter()
        print(f'Time taken: {time2 - time1}s\n')


        time1 = time.perf_counter()
        opti_seqs = optimize_sequences(to_optimize, identity_thresh, guide_length, ribo_seq, to_keep_single_targets,
                                       fileout=fileout, file=folder_to_save, score_type=score_type, msa_fast=msa_fast)
        time2 = time.perf_counter()
        end = time.perf_counter()
        print(f'Time taken: {time2 - time1}s\n')
        print(f'Time taken overall: {end - start}s\n')
        print('########################################################\n')
        return opti_seqs

    else:
        ranked_sorted_IGS = find_repeat_targets(new_data, ribo_seq, fileout=fileout, file=folder_to_save)
        time2 = time.perf_counter()
        print(f'Time taken: {time2 - time1}s\n')
        print('All guide sequences generated.')
        end = time.perf_counter()
        print(f'Time taken overall: {end - start}s\n')
        print('########################################################\n')
        return ranked_sorted_IGS


def find_cat_sites(target_names_and_seqs: list[tuple], igs_length: int, guide_length: int, min_length: int):
    """
    Finds all instances of a U or T in a set of sequences and creates ribozyme designs for these sites.

    :param target_names_and_seqs: list of tuples, each formatted as (string, Seq object) for all sequences.
    :param igs_length: desired IGS sequence length.
    :param guide_length: desired guide binding sequence length.
    :param min_length: minimum guide sequence length from 3' end. Must be smaller than guide_length.
    :return:
    """
    #   (a.k.a. minimum guide sequence length from 3' end) ex: if you want your guide sequence to
    #   bind to at least 35 nt at the 3' end of the target sequence, set min_length = 35.

    # initialize final product - will be a list of tuples of the format:
    # [target_name, target_sequence, (IGS, guide_seq, IGS_idx]
    # where each big list is one target sequence, and each tuple is one individual catalytic site

    data = [None] * len(target_names_and_seqs)  # one entry per target sequence
    col = 0
    # run function to find the index of all U residues of each sequence in target_seqs
    for name, sequ in target_names_and_seqs:
        # find all possible splice sites
        idx = find(sequ, 'U')
        # remove indexes that are < guide_length or > len(sequ) - igs_length (must have enough residues to attach to)
        idx_new = [res for res in idx if igs_length <= res < (len(sequ) - min_length)]

        if not idx_new:
            print(f'No viable catalytic sites in {name}')
            col += 1
            continue
        IGSes = [None] * len(idx_new)
        guides = [None] * len(idx_new)
        indexes = [None] * len(idx_new)
        small_col = 0
        for i in idx_new:
            # generate complementary guide sequence guide_length residues *downstream* of U site
            guide = sequ[i + 1:i + guide_length + 1].reverse_complement_rna()
            # generate complementary IGS sequence igs_length bases long *upstream* of U site
            IGS = sequ[i - igs_length:i].reverse_complement_rna()

            # now join IGS + G (where U site should be for wobble pair) + guide sequence and
            # append to ribozyme template
            IGSes[small_col] = IGS
            guides[small_col] = guide
            indexes[small_col] = i + 1  # we're adding 1 to the idx because of 0 based indexing

            small_col += 1
        data[col] = [name, sequ, (IGSes, guides, indexes)]
        col += 1
    # now remove entries with no viable sites
    # from https://www.geeksforgeeks.org/python-remove-none-values-from-list/
    filtered_data = list(filter(lambda item: item is not None, data))

    return filtered_data


def align_to_ref(data, ref_name_and_seq, base_to_find: str = 'U'):
    """
    Aligns each target sequence to a reference sequence, finds all catalytic site indexes, and returns:

    a list of all analyzed data with new indices wherein each entry has the following formatting:
    [name, sequ, (ribozyme_seq, IGS_and_guide_seq, IGS, guide_seq, og_idx, ref_idx)]

    a nested dictionary containing the index of each catalytic site to match ite with its reference index site.
    Formatting: {name: {original_idx: reference_idx}}

    aligned sequences

    :param data:
    :param ref_name_and_seq:
    :param base_to_find:
    """

    # Ok doing this because biopython is deprecating pairwise2 and I do not have the bandwidth to deal with that rn

    in_data = [(name, sequ, cat_site_info, ref_name_and_seq, base_to_find) for name, sequ, cat_site_info in data]

    print(f'Now re-indexing sequences to reference {ref_name_and_seq[0].replace("_", " ")}...')
    with Pool() as pool:
        new_data = pool.starmap(align_to_ref_loop, in_data)

    return new_data


def align_to_ref_loop(name, sequ, cat_site_info, ref_name_and_seq, base_to_find):
    """

    :param name:
    :param sequ:
    :param cat_site_info:
    :param ref_name_and_seq:
    :param base_to_find:
    :return:
    """
    # prepare patterns to look for to extract individual sequences from pairwise alignment
    pattern_a = 'seqA=\'(.*?)\''
    pattern_b = 'seqB=\'(.*?)\''

    # will have to keep in mind the potential lengths of the sequences and add length igs_length to our final E.coli index
    alignments = pairwise2.align.globalxx(sequ, ref_name_and_seq[1])
    new_aligns = (''.join(str(alignments[0])))

    seq_a = re.search(pattern_a, new_aligns).group(1)
    seq_b = re.search(pattern_b, new_aligns).group(1)

    # obtain index of new U
    idx_seq_a = find(seq_a, base_to_find)

    # initialize lists
    temp_IGSes = [None] * len(cat_site_info[-1])
    temp_guide_sequences = [None] * len(cat_site_info[-1])
    temp_og_and_ref_idexes = [None] * len(cat_site_info[-1])

    small_col = 0

    for idx in idx_seq_a:
        if small_col >= len(cat_site_info[-1]):
            break
        og_idx = cat_site_info[2][small_col]
        seq_a_idx = len(seq_a[:idx].replace('-', '')) + 1
        if seq_a_idx != og_idx:
            continue
        # find what index that is based on the reference sequence
        ref_string = seq_b[:idx]
        ref_idx = len(ref_string.replace('-', '')) + 1  # turns zero based indexing to ref_seq numbering

        # Use location of index to pair with old U idx. remember cat_site_info is a tuple of the form
        # (ribozyme_seq, IGS_and_guide_seq, IGS, guide_seq, IGS_idx)
        # where each entry is a list

        # fill lists with data
        temp_IGSes[small_col] = cat_site_info[0][small_col]
        temp_guide_sequences[small_col] = cat_site_info[1][small_col]
        temp_og_and_ref_idexes[small_col] = (og_idx, ref_idx)
        small_col += 1

    new_data = [name, sequ, (temp_IGSes, temp_guide_sequences, temp_og_and_ref_idexes)]
    return new_data


def read_fasta(in_file: str, file_type: str = 'fasta') -> list[tuple]:
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
        fasta_iter = SeqIO.parse(in_file, file_type)
        for record in fasta_iter:
            target_seqs_and_names.append((record.id, record.seq.upper().transcribe()))
        return target_seqs_and_names
    else:
        for filename in glob.glob(os.path.join(in_file, '*.' + file_type)):
            fasta_iter = SeqIO.parse(filename, file_type)
            for record in fasta_iter:
                target_seqs_and_names.append((record.id, record.seq.upper().transcribe()))
        return target_seqs_and_names


def find(string_to_analyze: str, char_to_find: str) -> list[int, str]:
    """
    Finds all instances of a character char_to_find in a string string_to_analyze.

    :param string_to_analyze:
    :param char_to_find:
    :return: A list of the indices and strings where the desired character is found.
    """
    return [i for i, ltr in enumerate(string_to_analyze) if ltr == char_to_find]


def transcribe_seq_file(seq_file: str) -> Seq:
    """
    Converts DNA sequences into RNA transcripts

    :param seq_file: .txt file containing one or more DNA sequences
    :return: Seq object
    """
    # seq_file must be a .txt file
    with open(seq_file) as f:
        for i in f:
            out_seq = Seq(i).upper().transcribe()
    return out_seq


def filter_for_targets(new_data: list[list], background_sequences_data: list[list], guide_length: int, ribo_seq: str,
                       ref_name_and_seq, optimize=True, min_true_cov=0.7, identity_thresh=0.7, fileout=False,
                       folder_to_save='', score_type='quantitative', msa_fast=False, keep_single_targets=False,
                       max_off_target=0.1):
    # Takes a list full of aligned target sequences and another with the background sequences which we do not want.
    # new_data items are formatted as: [name, sequ, (IGSes, guides, (og_idx, ref_idx))]
    # background_sequences_data items are formatted as: [name, sequ, (IGSes, guides, og_indexes)]

    # Find IGSes that only occur in new_data regardless of position
    # Extract IGSes
    print('First attempting to find IGSes that only occur in target sequences...')
    target_igs_set = {y for (_, _, cat_site_data) in new_data for y in cat_site_data[0]}
    background_igs_set = {y for (_, _, cat_site_data) in background_sequences_data for y in cat_site_data[0]}
    ideal_target_igs_set = target_igs_set.difference(background_igs_set)

    # If there are IGSes that occur only in targeted sequences, then awesome those are highest priority! check if they
    # meet our min_true_cov.
    # Keep track of how many targets don't have an IGS in the list. It's possible that we can skip to the next check
    # without having to call prep_for_optimizing yet.
    if ideal_target_igs_set:
        out_data = filter_for_targets_loop(new_data, ideal_target_igs_set, guide_length, ribo_seq, True, optimize,
                                           min_true_cov, identity_thresh, fileout, folder_to_save, score_type,
                                           msa_fast, keep_single_targets)
        if out_data:
            return out_data

    # If not, find IGSes that occur only in targeted sequences at the same position.
    # find the position of background sequence IGSes:
    print('\nChecking for IGSes that occur only in targeted sequences at the same position... ')
    aligned_background_sequences = align_to_ref(background_sequences_data, ref_name_and_seq)
    target_igs_set = {(y, x) for (_, _, cat_site_data) in new_data for y, x in zip(cat_site_data[0], cat_site_data[2][1])}
    background_igs_set = {(y, x) for (_, _, cat_site_data) in aligned_background_sequences for y, x
                          in zip(cat_site_data[0], cat_site_data[2][1])}

    ideal_target_igs_set = target_igs_set.difference(background_igs_set)

    if ideal_target_igs_set:
        out_data = filter_for_targets_loop(new_data, ideal_target_igs_set, guide_length, ribo_seq, False, optimize,
                                           min_true_cov, identity_thresh, fileout, folder_to_save, score_type,
                                           msa_fast, keep_single_targets)
        if out_data:
            return out_data

    # If still not, find IGSes with the most occurrences in the targeted sequences and the least occurrences in the
    # background sequences.

    if optimize:
        print('Could not find targets that do not hit background sequences. Finding targets with the least occurrences '
              'in background sequences...')
        new_data.extend(aligned_background_sequences)
        to_optimize, to_keep_single_targets = prep_for_optimizing(new_data, min_true_cov=min_true_cov,
                                                                  accept_single_targets=keep_single_targets)
        # Here add a step that goes through designs and finds those that occur mostly in target sequences


        opti_seqs = optimize_sequences(to_optimize, identity_thresh, guide_length, ribo_seq, to_keep_single_targets,
                                       fileout=fileout, file=folder_to_save, score_type=score_type,
                                       msa_fast=msa_fast)
        return opti_seqs

    else:
        print('Could not find targets that do not hit background sequences. Consider optimizing to find targets with '
              'the least occurrences in background sequences.\nNow finding single targets...')
        out_data = find_repeat_targets(new_data, ribo_seq, fileout=fileout, file=folder_to_save)
        return out_data


def filter_for_targets_loop(new_data: list, ideal_target_igs_set: set, guide_length: int, ribo_seq: str, testing_round_1: bool,
                            optimize=True, min_true_cov=0.7, identity_thresh=0.7, fileout=False, folder_to_save='',
                            score_type='quantitative', msa_fast=False, keep_single_targets=False):
    print('Found candidate IGSes! Now checking for minimum coverage requirements...\n')
    no_matches = len(new_data)
    filtered_data = []
    for name, sequ, target_data in new_data:
        igs_list_temp = []
        guide_list_temp = []
        idx_list_temp = []
        igs_candidates, guide_candidates, indexes_candidates = target_data
        for igs, guide, indexes in zip(igs_candidates, guide_candidates, indexes_candidates):
            if testing_round_1:
                if igs in ideal_target_igs_set:
                    igs_list_temp.append(igs)
                    guide_list_temp.append(guide)
                    idx_list_temp.append(indexes)
            else:
                if (igs, indexes[1]) in ideal_target_igs_set:
                    igs_list_temp.append(igs)
                    guide_list_temp.append(guide)
                    idx_list_temp.append(indexes)

        if igs_list_temp:
            filtered_data.append([name, sequ, (igs_list_temp, guide_list_temp, idx_list_temp)])
        else:
            # If the target did not contain any igs, then track that
            no_matches -= 1
            if no_matches / len(new_data) < min_true_cov:
                print('No candidates within the minimum true coverage were found.')
                return None

    if optimize:
        print('Candidates within minimum true coverage found! Now optimizing...')
        to_optimize, to_keep_single_targets = prep_for_optimizing(filtered_data, min_true_cov=min_true_cov,
                                                                  accept_single_targets=keep_single_targets)
        if to_optimize:
            opti_seqs = optimize_sequences(to_optimize, identity_thresh, guide_length, ribo_seq, to_keep_single_targets,
                                           fileout=fileout, file=folder_to_save, score_type=score_type,
                                           msa_fast=msa_fast)
            return opti_seqs
        elif to_keep_single_targets:
            print('Only found single targets.')
            return to_keep_single_targets
        else:
            print('No candidate IGSes meet minimum coverage requirements. Continuing to search...')
            return None
    print('Candidates within minimum true coverage found! Now ranking...')
    ranked_sorted_IGS = find_repeat_targets(new_data, ribo_seq, fileout=fileout, file=folder_to_save)
    if ranked_sorted_IGS:
        return ranked_sorted_IGS
    else:
        print('No candidate IGSes meet minimum coverage requirements. Continuing to search...')
        return None

def ribochecker(opti_seqs, bad_seqs):


    return


def prep_for_optimizing(new_data: list[list], min_true_cov: float = 0, accept_single_targets: bool = True):
    """
    :param new_data:
    :param min_true_cov:
    :param accept_single_targets:
    :return:
    """
    # recall that new_data is a list of lists with one entry per target sequence where each entry is of the form:
    # [name, sequ, (IGSes, guide_sequences, (og_idx, ref_idx))]
    big_repeats = []
    start = time.perf_counter()
    print('Finding repeat subsets...')

    igs_subsets = [set(cat_site_data[0]) for i, (_, _, cat_site_data) in enumerate(new_data)]
    # igs_subsets_pairs = [(igs_subsets[i], igs_subsets[i+1:]) for i in range(len(igs_subsets) -1)]

    for i, igs_data_a in enumerate(igs_subsets):
        if i == len(igs_subsets):
            break
        for igs_data_b in igs_subsets[i + 1:]:
            # make a subset of second column of unique IGS values
            # and find the shared values between sets with no duplicates
            no_dupes = igs_data_a & igs_data_b
            # remove any nans
            repeats = [item for item in no_dupes if str(item) != 'nan']
            if repeats:
                big_repeats.extend(repeats)

    print(f'Time taken: {time.perf_counter() - start}s\n')
    print('Found repeat subsets. Now analyzing sequences...')

    # remove duplicates of all IGSes found
    filtered_list = list(set(big_repeats))

    to_optimize = defaultdict(list)  # will keep those IGSes that meet our minimum percent coverage later
    to_keep_single_targets = defaultdict(list)

    # Now find the IGS sequences in the original target sequences and extract matching data
    big_temp_list = []
    for IGS_sequ in filtered_list:  # for each IGS sequence
        temp_list = []
        coverage_count = 0  # this will track how many target sequences contain the IGS
        col = 0
        for org, sequ, cat_site_data in new_data:  # in each target sequence (target seq is org)
            pos = [i for i, e in enumerate(cat_site_data[0]) if e == IGS_sequ]  # extract positions of matching IGSes

            if not pos:  # if not in this particular seq, try next seq
                col += 1
                continue
            else:  # if the IGS is found in the target seq, get all the data we want
                # To calculate whether a particular IGS is in the same position in multiple targets or not I am making a
                # HUGE assumption: as the sequences are already aligned to a single reference sequence, any base
                # position will have EXACTLY the same position numbering. This allows my code to use dictionaries to
                # determine whether something is or is not on target as the position must match exactly. If we want some
                # tolerance in how many base pairs away a position can be to be considered on-target I will need to find
                # another way to calculate this which would probably use a LOT more computational power.
                if coverage_count == 0:
                    on_target_count = defaultdict(list)
                target_num = len(pos)

                for p in pos:
                    ref_pos = cat_site_data[2][p][1]
                    if on_target_count[ref_pos]:
                        on_target_count[ref_pos] += 1
                    else:
                        on_target_count[ref_pos] = 1

                    guide_id = str(IGS_sequ) + str(
                        ref_pos)  # save a guide ID that is the IGS sequence and the reference position
                    temp_list.append([IGS_sequ, None, None, None, org, target_num, cat_site_data[2][p][0],
                                      ref_pos, cat_site_data[1][p], sequ, guide_id])

                coverage_count += 1
            col += 1

        # if the IGS sequence wasn't found in enough target seqs, add to list of seqs to filter
        perc_coverage = coverage_count / len(new_data)

        for item in temp_list:
            item[1] = perc_coverage
            item[2] = on_target_count[item[7]] / coverage_count  # out of all organisms with this IGS at this position
            true_coverage = on_target_count[item[7]] / len(new_data)
            item[3] = true_coverage
            big_temp_list.append(list(item))
            # if the true % coverage is above min_true_cov, and more than one sequence, mark this IGS for optimization
            if true_coverage >= min_true_cov and on_target_count[item[7]] > 1:
                # save the IGS sequence AND the matching index in the reference sequence
                to_optimize[item[10]].append(item)
            # Else if the true % coverage is above min_true_cov but only one sequence still keep but will not optimize
            elif true_coverage >= min_true_cov and accept_single_targets:
                to_keep_single_targets[item[10]].append(item)

    return to_optimize, to_keep_single_targets


def find_repeat_targets(new_data: list[list], ribo_seq: str, fileout: bool = False, file: str = ''):
    """

    :param new_data:
    :param ribo_seq:
    :param fileout:
    :param file:
    :return:
    """

    # recall that new_data is a list of lists with one entry per target sequence where each entry is of the form:
    # [name, sequ, (IGSes, guide_sequences, (og_idx, ref_idx))]
    big_repeats = []
    start = time.perf_counter()
    print('Finding repeat subsets...')

    igs_subsets = [set(cat_site_data[0]) for i, (_, _, cat_site_data) in enumerate(new_data)]
    # igs_subsets_pairs = [(igs_subsets[i], igs_subsets[i+1:]) for i in range(len(igs_subsets) -1)]

    for i, igs_data_a in enumerate(igs_subsets):
        if i == len(igs_subsets):
            break
        for igs_data_b in igs_subsets[i+1:]:
            # make a subset of second column of unique IGS values
            # and find the shared values between sets with no duplicates
            no_dupes = igs_data_a & igs_data_b
            # remove any nans
            repeats = [item for item in no_dupes if str(item) != 'nan']
            if repeats:
                big_repeats.extend(repeats)

    print(f'Time taken: {time.perf_counter() - start}s\n')
    print('Found repeat subsets. Now analyzing sequences...')

    # remove duplicates of all IGSes found
    filtered_list = list(set(big_repeats))

    # Now find the IGS sequences in the original target sequences and extract matching data
    big_temp_list = []
    for IGS_sequ in filtered_list:  # for each IGS sequence
        temp_list = []
        coverage_count = 0  # this will track how many target sequences contain the IGS
        col = 0
        for org, sequ, cat_site_data in new_data:  # in each target sequence (target seq is org)
            pos = [i for i, e in enumerate(cat_site_data[0]) if e == IGS_sequ]  # extract positions of matching IGSes

            if not pos:  # if not in this particular seq, try next seq
                col += 1
                continue
            else:  # if the IGS is found in the target seq, get all the data we want
                # To calculate whether a particular IGS is in the same position in multiple targets or not I am making a
                # HUGE assumption: as the sequences are already aligned to a single reference sequence, any base
                # position will have EXACTLY the same position numbering. This allows my code to use dictionaries to
                # determine whether something is or is not on target as the position must match exactly. If we want some
                # tolerance in how many base pairs away a position can be to be considered on-target I will need to find
                # another way to calculate this which would probably use a LOT more computational power.
                if coverage_count == 0:
                    on_target_count = defaultdict(list)
                target_num = len(pos)

                for p in pos:
                    ref_pos = cat_site_data[2][p][1]
                    if on_target_count[ref_pos]:
                        on_target_count[ref_pos] += 1
                    else:
                        on_target_count[ref_pos] = 1

                    guide_id = str(IGS_sequ) + str(
                        ref_pos)  # save a guide ID that is the IGS sequence and the reference position

                    guide_and_igs = cat_site_data[1][p] + 'G' + IGS_sequ
                    full_design = guide_and_igs + ribo_seq
                    temp_list.append([IGS_sequ, None, None, None, org, target_num, cat_site_data[2][p][0],
                                      ref_pos, cat_site_data[1][p], guide_and_igs,
                                      full_design, sequ, guide_id])
                coverage_count += 1
            col += 1

        # if the IGS sequence wasn't found in enough target seqs, add to list of seqs to filter
        perc_coverage = coverage_count / len(new_data)

        for item in temp_list:
            item[1] = perc_coverage
            item[2] = on_target_count[item[7]] / coverage_count  # out of all organisms with this IGS at this position
            true_coverage = on_target_count[item[7]] / len(new_data)
            item[3] = true_coverage
            big_temp_list.append(list(item))

    # Make dataframes for excel ig
    ranked_IGS = pd.DataFrame(data=big_temp_list, index=None,
                              columns=['IGS sequence', '%' + ' coverage', '% on target in targets covered',
                                       'True % coverage', 'Target name', 'Occurrences in Target Sequence',
                                       'Index of Splice Site', 'Equivalent Reference Index of Splice Site',
                                       'Just Guide',
                                       'Guide + G + IGS', 'Ribozyme Design', 'Original sequence', 'ID'])

    ranked_sorted_IGS = ranked_IGS.sort_values(by=['%' + ' coverage', '% on target in targets covered'],
                                               ascending=[False, False])

    if fileout:
        # make a new file
        if os.path.exists(f'{file}/Ranked Ribozyme Designs with Raw Guide Sequence Designs.csv'):
            os.remove(f'{file}/Ranked Ribozyme Designs with Raw Guide Sequence Designs.csv')

        ranked_sorted_IGS.to_csv(f'{file}/Ranked Ribozyme Designs with Raw Guide Sequence Designs.csv',
                                 index=False)


    return ranked_sorted_IGS


def optimize_sequences(to_optimize, thresh: float, guide_length: int, ribo_seq: str, single_targets,
                       fileout: bool = False, file: str = '', score_type: str = 'quantitative',
                       gaps_allowed: bool = True, msa_fast: bool = False):
    """

    :param to_optimize:
    :param thresh:
    :param guide_length:
    :param ribo_seq:
    :param single_targets:
    :param fileout:
    :param file:
    :param score_type:
    :param gaps_allowed:
    :param msa_fast:
    :return:
    """
    print(f'Optimizing {len(to_optimize)} guide sequences...')

    in_data = [(key, to_optimize[key], thresh, score_type, gaps_allowed, guide_length, ribo_seq, msa_fast) for key in to_optimize]

    with Pool() as pool:
        opti_seqs = pool.starmap(optimize_sequences_loop, in_data)

    # If there are any single targets to keep (depends on min true coverage setting) add them here
    if single_targets:
        print(f'Storing {len(single_targets)} single target sequences. If you do not want single target guides, '
              f'please increase your min true coverage parameter.')
        for key in single_targets.keys():
            guide = single_targets[key][0][8]
            design_sequence = guide + 'G' + re.sub(r'\d+', '', key)
            ribo_design = design_sequence + ribo_seq

            # set score to nan
            opti_seqs.append([single_targets[key][0][0], single_targets[key][0][7], 1, single_targets[key][0][1],
                              single_targets[key][0][2], single_targets[key][0][3],
                              [(target[4], target[6], target[5] - 1) for target in single_targets[key]], guide,
                              design_sequence, ribo_design])

    if fileout:
        if os.path.exists(f'{file}/Ranked Ribozyme Designs with Optimized Guide Sequence Designs {score_type}.csv'):
            os.remove(f'{file}/Ranked Ribozyme Designs with Optimized Guide Sequence Designs {score_type}.csv')

        with open(f'{file}/Ranked Ribozyme Designs with Optimized Guide Sequence Designs {score_type}.csv', 'w') as f:
            f.write('IGS,Reference index,Score,% cov,% on target,True % cov,(Target name| Target idx| Other occurrences '
                    'of IGS in target sequence),Optimized guide,Optimized guide + G + IGS,Full Ribozyme design\n')
            for item in opti_seqs:
                list_for_csv = str(item[6]).replace(',', '|')
                f.write(f'{item[0]},{item[1]},{item[2]},{item[3]},{item[4]},{item[5]},{list_for_csv},{item[7]},{item[8]},'
                        f'{item[9]}\n')

    print('All guide sequences optimized.')
    return opti_seqs


def optimize_sequences_loop(name, items, thresh, score_type, gaps_allowed, guide_length, ribo_seq, msa_fast):
    guides_to_optimize = [target[8] for target in items]

    # do a MSA of the sequences and optimize the guide sequence with this MSA
    opti_seq, score = msa_and_optimize(name, guides_to_optimize, thresh, score_type, gaps_allowed, msa_fast)

    # truncate optimized sequence to the desired guide sequence length
    truncated_guide = opti_seq[-guide_length:]

    # Now make an actual full design with the length needed
    design_sequence = truncated_guide + 'G' + re.sub(r'\d+', '', name)  # our ID is the IGS and the ref index so remove the index
    ribo_design = design_sequence + ribo_seq

    # store these designs and score
    # store the data as follows: IGS, reference idx, score, % cov, % on target, true % cov, [(target name, target idx,
    # occurrences of IGS in target)], truncated_guide, design_sequence, ribo_design
    opti_seqs = [items[0][0], items[0][7], score, items[0][1], items[0][2], items[0][3],
                 [(target[4], target[6], target[5] - 1) for target in items], truncated_guide, design_sequence,
                 ribo_design]

    return opti_seqs


def msa_and_optimize(name, seqs_to_align, thresh: float = 0.7, score_type: str = 'quantitative',
                     gaps_allowed: bool = True, msa_fast: bool = False):
    """

    :param name:
    :param seqs_to_align:
    :param thresh:
    :param score_type:
    :param gaps_allowed:
    :param msa_fast:
    :return:
    """
    # seqs_to_align is a list containing a sequence from an individual organism in each position.

    with open(f'to_align_{name}.fasta', 'w') as f:
        for line in range(len(seqs_to_align)):
            f.write('>seq' + str(line) + '\n' + str(seqs_to_align[line]) + '\n')

    muscle_exe = 'muscle5'

    seqs_str = ''
    for i, seq in enumerate(seqs_to_align):
        seqs_str += f'>seq{i}\n{str(seqs_to_align[i])}\n'

    # seqs_str_write = (SeqRecord(seq, id=f'seq{i}', name='_', description='_') for i, seq in enumerate(seqs_to_align))
    #
    # # now align the sequences.
    # with subprocess.Popen([muscle_exe], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    #                       universal_newlines=True) as child:
    #
    #     SeqIO.write(seqs_str_write, child.stdin, 'fasta')
    #     msa = AlignIO.read(child.stdout, 'fasta')
    #     print(msa)
    #
    #     child.stdin.write(seqs_str)
    #     child_out = StringIO(child.communicate()[0])
    #     print(child_out)
    #     seqs_aligned = list(SeqIO.parse(child_out, format='fasta'))
    #     print(seqs_aligned)
    #     msa = AlignIO.read(child_out, 'fasta')
    #     print(msa)
    # # msa = AlignIO.read(seqs_aligned, format='fasta')
    # # msa = AlignIO.read(child_out, 'fasta')
    # # print(msa)
    # summary_align = AlignInfo.SummaryInfo(msa)
    # print(summary_align.gap_consensus(threshold=thresh, ambiguous='N'))
    # child = subprocess.check_output([muscle_exe, "-align"],
    #                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    # child.stdin.write(seqs_str.encode())
    # child_out = child.communicate()[0].decode('utf8')

    if msa_fast:
        subprocess.check_output([muscle_exe, "-super5", f'to_align_{name}.fasta', "-output", f'aln_{name}.afa'],
                                stderr=subprocess.DEVNULL)
    else:
        subprocess.check_output([muscle_exe, "-align", f'to_align_{name}.fasta', "-output", f'aln_{name}.afa'],
                                stderr=subprocess.DEVNULL)

    msa = AlignIO.read(open(f'aln_{name}.afa'), 'fasta')

    # delete files:
    if os.path.exists(f'aln_{name}.afa'):
        os.remove(f'aln_{name}.afa')
    if os.path.exists(f'to_align_{name}.fasta'):
        os.remove(f'to_align_{name}.fasta')

    # if needed at this point, get alignment
    if score_type != 'quantitative':
        summary_align = AlignInfo.SummaryInfo(msa)
        if gaps_allowed:
            opti_seq = summary_align.gap_consensus(threshold=thresh, ambiguous='N')
        else:
            opti_seq = summary_align.dumb_consensus(threshold=thresh, ambiguous='N')
        n = len(opti_seq)

    if score_type == 'naive':
        # Naively tell me what the percent identity is: 1- ambiguity codons/ length
        score = 1 - str(opti_seq).count('N') / n

    elif score_type == 'weighted':
        # prepare scoring matrix a(x) (currently only uses N and gap)
        a = {
            'A': 1, 'T': 1, 'G': 1, 'C': 1, 'U': 1,
            'R': 0.5, 'Y': 0.5, 'M': 0.5, 'K': 0.5, 'S': 0.5, 'W': 0.5,
            'H': 0.3, 'B': 0.3, 'D': 0.3, 'V': 0.3,
            'N': 0, '-': 0
        }

        # score based on extended identity and position of the bases
        weight_sum = 0
        weight_score_sum = 0
        i = 0  # first position

        for x in str(opti_seq):
            w = exp(n - i)  # Calculate weight based on proximity to IGS
            weight_sum += w
            weight_score_sum += w * a[x]
            i += 1  # move downstream
        score = weight_score_sum / weight_sum

    elif score_type == 'quantitative':
        score, opti_seq = get_quantitative_score(msa, count_gaps=gaps_allowed)

    # return optimized_guide, score
    return opti_seq, score


def get_quantitative_score(msa, thresh: float = 0.7, chars_to_ignore: list[str] = None,
                           count_gaps: bool = True, penalize_trailing_gaps: bool = False):
    """

    :param msa:
    :param thresh:
    :param chars_to_ignore:
    :param count_gaps:
    :param penalize_trailing_gaps:
    :return:
    """
    # a modified version of BioPython'string_to_analyze PSSM funtction that counts gaps by default and also gives us
    # the quantitative score. This quantitative score is as follows:

    if chars_to_ignore is None:
        chars_to_ignore = []
    if not isinstance(chars_to_ignore, list):
        raise TypeError('chars_to_ignore should be a list.')

    gap_char = '-'
    if not count_gaps:
        chars_to_ignore.append(gap_char)

    summary_align = AlignInfo.SummaryInfo(msa)
    if count_gaps:
        left_seq = summary_align.gap_consensus(threshold=thresh, ambiguous='N')
    else:
        left_seq = summary_align.dumb_consensus(threshold=thresh, ambiguous='N')

    sum_probabilities = 0
    seqs_list = [record.seq for record in msa]  # had to convert to list because biopython refuses to work
    # now start looping through all of the sequences and getting info
    seq_num = len(seqs_list)
    for residue_num in range(len(left_seq)):
        pos_prob = 0
        for record in seqs_list:
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

    opti_seq = left_seq.strip('-')
    score = sum_probabilities / len(left_seq)

    return score, opti_seq


def calc_shannon_entropy(target_names_and_seqs, ref_name_and_seq, base: float = None, count_gap: bool = False,
                         fileout: bool = False, file: str = ''):
    """

    :param target_names_and_seqs:
    :param ref_name_and_seq:
    :param base:
    :param count_gap:
    :param fileout:
    :param file:
    :return:
    """
    # prepare patterns to look for to extract individual sequences from pairwise alignment
    pattern_a = 'seqA=\'(.*?)\''
    pattern_b = 'seqB=\'(.*?)\''

    if count_gap:
        states = ['A', 'U', 'C', 'G', '-']
    else:
        states = ['A', 'U', 'C', 'G']

    probs = {pos: {'A': 0, 'U': 0, 'C': 0, 'G': 0, '-': 0} for pos in range(len(ref_name_and_seq[1]))}

    col = 0
    new_aligns = [None] * len(target_names_and_seqs)

    for name, sequ in target_names_and_seqs:
        # will have to keep in mind the potential lengths of the sequences and add length igs_length to our final E.coli index

        alignments = pairwise2.align.globalxx(sequ, ref_name_and_seq[1])
        new_aligns[col] = (''.join(str(alignments[0])))

        # will have to keep in mind the potential lengths of the sequences and add length igs_length to our final E.coli index
        seq_a = re.search(pattern_a, new_aligns[col]).group(1)
        seq_b = re.search(pattern_b, new_aligns[col]).group(1)

        # set up range to check over
        idx_seq_a = range(len(seq_a))

        for idx in idx_seq_a:
            # find what index that is based on the reference sequence
            ref_string = seq_b[:idx + 1]
            nucleotide = seq_a[idx]
            ref_idx = len(
                ref_string.replace('-', '')) - 1  # remember that the probs dictionary is in zero-based indexing
            if ref_idx < 0:  # if the sequence is bigger than the original sequence
                continue
            if nucleotide not in ['A', 'U', 'C', 'G', '-']:
                if nucleotide == 'T':
                    probs[ref_idx]['U'] += 1
                elif nucleotide == 'M':
                    probs[ref_idx]['A'] += 0.5
                    probs[ref_idx]['C'] += 0.5
                elif nucleotide == 'R':
                    probs[ref_idx]['A'] += 0.5
                    probs[ref_idx]['G'] += 0.5
                elif nucleotide == 'W':
                    probs[ref_idx]['A'] += 0.5
                    probs[ref_idx]['U'] += 0.5
                elif nucleotide == 'S':
                    probs[ref_idx]['G'] += 0.5
                    probs[ref_idx]['C'] += 0.5
                elif nucleotide == 'Y':
                    probs[ref_idx]['U'] += 0.5
                    probs[ref_idx]['C'] += 0.5
                elif nucleotide == 'K':
                    probs[ref_idx]['G'] += 0.5
                    probs[ref_idx]['U'] += 0.5
                elif nucleotide == 'V':
                    probs[ref_idx]['A'] += 1 / 3
                    probs[ref_idx]['C'] += 1 / 3
                    probs[ref_idx]['G'] += 1 / 3
                elif nucleotide == 'H':
                    probs[ref_idx]['A'] += 1 / 3
                    probs[ref_idx]['C'] += 1 / 3
                    probs[ref_idx]['U'] += 1 / 3
                elif nucleotide == 'D':
                    probs[ref_idx]['A'] += 1 / 3
                    probs[ref_idx]['G'] += 1 / 3
                    probs[ref_idx]['U'] += 1 / 3
                elif nucleotide == 'B':
                    probs[ref_idx]['C'] += 1 / 3
                    probs[ref_idx]['G'] += 1 / 3
                    probs[ref_idx]['U'] += 1 / 3
                elif nucleotide == 'N':
                    probs[ref_idx]['A'] += 0.25
                    probs[ref_idx]['C'] += 0.25
                    probs[ref_idx]['U'] += 0.25
                    probs[ref_idx]['G'] += 0.25
            else:
                probs[ref_idx][nucleotide] += 1

        col += 1

    shannon_entropy_list = {pos: 0 for pos in range(len(ref_name_and_seq[1]))}

    # The original Shannon entropy assumes a log base 2 as it deals with bits, but we can use any tbh. We'll default to
    # the natural log
    for pos in range(len(shannon_entropy_list)):
        for nucleotide in states:
            # the probability or base is the number of occurrences of said base at position pos over the
            # total possible bases
            p = probs[pos][nucleotide] / len(states)
            if p != 0:
                if base:
                    shannon_entropy_list[pos] += -p * log(p, base)
                else:
                    shannon_entropy_list[pos] += -p * log(p)
    xaxis_list = [pos + 1 for pos in shannon_entropy_list.keys()]
    yaxis_list = [shannon_entropy_list[pos] for pos in shannon_entropy_list.keys()]

    if fileout:
        sorted_opti_seqs = pd.DataFrame(data={'Reference index': xaxis_list, 'Shannon entropy': yaxis_list}, index=None)
        sorted_opti_seqs.to_csv(f'{file}/Shannon entropy of aligned sequences.csv', index=False)

    return shannon_entropy_list, xaxis_list, yaxis_list


    def ribo_checker(designs, target_sequences):
        """Checks generated designs against a set of sequences to either find how well they align to a given dataset.
        Will return a set of sequences that match for each design as well as a score showing how well they matched.
        :param target_sequences:
        :param designs:
        """
    #     find IGSes

    return

# def ambiguity_consensus(summary_align, threshold=0.7, gaps_allowed=False):
#     # I've made a modified version of BioPython'string_to_analyze dumb consensus that will allow IUPAC ambiguity codes
#
#     consensus = ''
#
#     # find the length of the consensus we are creating
#     con_len = summary_align.alignment.get_alignment_length()
#
#     # go through each seq item
#     for guide_length in range(con_len):
#         # keep track of the counts of the different atoms we get
#         atom_dict = Counter()
#         num_atoms = 0
#
#         for record in summary_align.alignment:
#             # make sure we haven't run past the end of any sequences
#             # if they are of different lengths
#             try:
#                 c = record[guide_length]
#             except IndexError:
#                 continue
#             if c != '-' and c != '.':
#                 atom_dict[c] += 1
#
#                 num_atoms += 1
#
#         max_atoms = []
#         max_size = 0
#
#         for atom in atom_dict:
#             if atom_dict[atom] > max_size:
#                 max_atoms = [atom]
#                 max_size = atom_dict[atom]
#             elif atom_dict[atom] == max_size:
#                 max_atoms.append(atom)
#
#         # Here is where I will change the code. Instead of looking at a threshold, we'll try to find the ambiguity code
#         # that matches the sequence
#         if (len(max_atoms) == 1) and (
#                 (float(max_size) / float(num_atoms)) >= threshold):
#             consensus += max_atoms[0]
#         else:
#             # Check whether it matches ambiguous codons
#             consensus += ambiguous
#
#     return Seq(consensus)
