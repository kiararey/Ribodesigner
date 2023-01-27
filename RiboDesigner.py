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

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SeqIO, pairwise2, AlignIO
from Bio.Seq import Seq
import re
import pandas as pd
from collections import defaultdict
import subprocess
from Bio.Align import AlignInfo
from math import exp, log
from multiprocessing import Pool
import time



def RiboDesigner(igs_length: int, guide_length: int, min_length: int, barcode_seq_file: str, ribobody_file: str,
                 target_sequences_folder: str, ref_sequence_file=None,
                 optimize_seq=True, min_true_cov=0.7, identity_thresh=0.7, fileout=False, folder_to_save='',
                 score_type='quantitative'):
    # RiboDesigner is a function that generates Ribozyme designs to target a set of sequences (target_seqs)

    # target_seqs: list of seqs, containing the target RNA sequences to target (5' -> 3'). Generate with read_fasta or
    # read_fasta_folder
    # target_names: list of strings, containing the corresponding names to target_seqs. Generate with read_fasta or
    # read_fasta_folder
    # ribobody: string, containing the Ribozyme sequence to use as the body template (5' -> 3')
    # bar_seq: string, with the desired insertion barcode sequence (5' -> 3')
    # min_length: int, must be smaller than guide_length. nucleotide tolerance at 3' end of target sequence
    #   (a.k.a. minimum guide sequence length from 3' end) ex: if you want your guide sequence to
    #   bind to at least 35 nt at the 3' end of the target sequence, set min_length = 35.
    # fileout: boolean, whether we want a csv file output or not (default True)
    # file_name: string, the path where the folder we will save our outputs in if fileout = True

    start = time.time()
    barcode_seq = transcribe_seq_file(barcode_seq_file)
    ribobody = transcribe_seq_file(ribobody_file)

    if target_sequences_folder[-6:] == '.fasta':
        target_names_and_seqs = read_fasta(target_sequences_folder)
    else:
        target_names_and_seqs = read_fasta_folder(target_sequences_folder)

    try:
        target_names_and_seqs[0]
    # if we hit an index error we've run out of sequence and
    # should not add new residues
    except IndexError:
        print(f'No sequences found in {target_sequences_folder}. Please make sure your files are not empty!\n')
        return None

    if not ref_sequence_file:
        # If we do not have a reference sequence, just choose one randomly
        print('No reference sequence provided. Picking a random sequence as reference...')
        ref_name_and_seq = random.choice(target_names_and_seqs)
    else:
        ref_name_and_seq = read_fasta(ref_sequence_file)[0]

    print(f'Found {len(target_names_and_seqs)} total sequences to analyze.')

    # Make the Ribozyme sequence by combining the main body and the barcode
    ribo_seq = ribobody + barcode_seq

    # find all catalytic U sites
    # Remember: data has one tuple per target sequence, the third entry is a tuple for each catalytic U site
    time1 = time.time()
    data = find_cat_sites(target_names_and_seqs, ribo_seq, igs_length, guide_length, min_length)
    time2 = time.time()
    print(f'Time taken: {time2 - time1}s\n')

    # Now align sequences to reference sequences and get the conversion dictionaries for each
    time1 = time.time()
    new_data = align_to_ref(data, ref_name_and_seq)
    time2 = time.time()
    print(f'Time taken: {time2 - time1}s\n')

    time1 = time.time()
    big_temp_list, to_optimize, filtered_list, \
    ranked_IGS, ranked_sorted_IGS, to_keep_single_targets = find_repeat_targets(new_data, min_true_cov=min_true_cov,
                                                                                fileout=fileout, file=folder_to_save)
    time2 = time.time()
    print(f'Time taken: {time2 - time1}s\n')

    # Now, we can optimize each sequence
    if optimize_seq:
        time1 = time.time()
        opti_seqs = optimize_sequences(to_optimize, identity_thresh, guide_length, ribo_seq, to_keep_single_targets,
                                       fileout=fileout, file=folder_to_save, score_type=score_type)
        time2 = time.time()
        end = time.time()
        print(f'Time taken: {time2 - time1}s\n')
        print(f'Time taken overall: {end - start}s\n')
        return opti_seqs

    else:
        print('All guide sequences generated.')
        end = time.time()
        print(f'Time taken overall: {end - start}s\n')
        return ranked_sorted_IGS


def find_cat_sites(target_names_and_seqs, ribo_seq, igs_length: int, guide_length, min_length: int):
    # find_cat_sites finds all instances of a U or T in a set of sequences and creates ribozyme designs for these
    # target_seqs_and_names (list of tuples: each with format (string, Seq object) are the sequences we want to analyze
    # ribo_seq (Seq object) is the sequence of the ribozyme body and barcode we are using
    # igs_length (int), desired IGS sequence length
    # guide_length (int), desired guide binding sequence length
    # min_length (int), must be smaller than guide_length. nucleotide tolerance at 3' end of target sequence
    #   (a.k.a. minimum guide sequence length from 3' end) ex: if you want your guide sequence to
    #   bind to at least 35 nt at the 3' end of the target sequence, set min_length = 35.

    # initialize final product - will be a list of tuples of the format:
    # [target_name, target_sequence, (IGS_and_guide_seq + ribo_seq, IGS_and_guide_seq, IGS, guide_seq, IGS_idx]
    # where each big list is one target sequence, and each tuple is one individual catalytic site

    data = [None] * len(target_names_and_seqs)  # one entry per target sequence
    col = 0
    # run function to find the index of all U residues of each sequence in target_seqs
    for name, sequ in target_names_and_seqs:
        # find all possible splice sites
        idx = find(sequ, 'U')

        if not idx:  # in case it'string_to_analyze a DNA sequence (so no U sites in sequence)
            idx = find(sequ, 'T')  # because U is T duh
        # remove indexes that are < guide_length or > len(sequ) - igs_length (must have enough residues to attach to)
        idx_new = [res for res in idx if igs_length <= res < (len(sequ) - min_length)]

        if not idx_new:
            print(f'No viable catalytic sites in {name}')
            col += 1
            continue

        IGS_and_ribo_seqs = [None] * len(idx_new)
        IGS_and_guide_seqs = [None] * len(idx_new)
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
            IGS_and_guide_seq = guide + 'G' + IGS  # Make 5'-> 3' binding seq: attaches to 5' end of ribo_seq correctly

            IGS_and_ribo_seqs[small_col] = IGS_and_guide_seq + ribo_seq
            IGS_and_guide_seqs[small_col] = IGS_and_guide_seq
            IGSes[small_col] = IGS
            guides[small_col] = guide
            indexes[small_col] = i + 1  # we're adding 1 to the idx because of 0 based indexing

            small_col += 1
        data[col] = [name, sequ, (IGS_and_ribo_seqs, IGS_and_guide_seqs, IGSes, guides, indexes)]
        col += 1
    # now remove entries with no viable sites
    # from https://www.geeksforgeeks.org/python-remove-none-values-from-list/
    filtered_data = list(filter(lambda item: item is not None, data))

    return filtered_data


def align_to_ref(data, ref_name_and_seq, base_to_find='U'):
    # align_to_ref will first align each target sequence to a reference sequence, find all catalytic site indexes, and
    # return a list of all the data analyzed with their new indexes, with one entry per target sequence as below:
    # [name, sequ, (ribozyme_seq, IGS_and_guide_seq, IGS, guide_seq, og_idx, ref_idx)]
    # will also return a dictionary of dictionaries (one inner dictionary for each target sequence) containing the
    # index of each catalytic site to match it with its reference index site.
    # the dictionary will have the form: conversion_dicts = {name: {original_idx: reference_idx}}
    # finally it will also return the alignments in case we need them downstream.

    # Ok doing this because biopython is deprecating pairwise2 and I do not have the bandwidth to deal with that rn

    in_data = [(name, sequ, cat_site_info, ref_name_and_seq, base_to_find) for name, sequ, cat_site_info in data]

    print(f'Now re-indexing target sequences to reference {ref_name_and_seq[0].replace("_", " ")}...')
    with Pool() as pool:
        new_data = pool.starmap(align_to_ref_loop, in_data)

    return new_data


def align_to_ref_loop(name, sequ, cat_site_info, ref_name_and_seq, base_to_find):
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
    temp_ribozyme_sequences = [None] * len(cat_site_info[4])
    temp_IGS_and_guide_sequences = [None] * len(cat_site_info[4])
    temp_IGSes = [None] * len(cat_site_info[4])
    temp_guide_sequences = [None] * len(cat_site_info[4])
    temp_og_and_ref_idexes = [None] * len(cat_site_info[4])

    small_col = 0

    for idx in idx_seq_a:
        if small_col >= len(cat_site_info[4]):
            break
        og_idx = cat_site_info[4][small_col]
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
        temp_ribozyme_sequences[small_col] = cat_site_info[0][small_col]
        temp_IGS_and_guide_sequences[small_col] = cat_site_info[1][small_col]
        temp_IGSes[small_col] = cat_site_info[2][small_col]
        temp_guide_sequences[small_col] = cat_site_info[3][small_col]
        temp_og_and_ref_idexes[small_col] = (og_idx, ref_idx)
        small_col += 1

    new_data = [name, sequ, (temp_ribozyme_sequences, temp_IGS_and_guide_sequences, temp_IGSes,
                             temp_guide_sequences, temp_og_and_ref_idexes)]
    return new_data


def read_fasta_folder(path, file_type='fasta'):
    # for now, filepath is where we are keeping all of the .FASTA files
    # will return a list of tuples as (ID, sequence)

    target_seqs_and_names = []
    for filename in glob.glob(os.path.join(path, '*.' + file_type)):
        fasta_iter = SeqIO.parse(filename, file_type)
        for record in fasta_iter:
            target_seqs_and_names.append((record.id, record.seq.upper().transcribe()))
    return target_seqs_and_names


def read_fasta(in_file: str, file_type='fasta'):
    # filepath here is a fasta file
    # will return a list of tuples as (ID, sequence)

    target_seqs_and_names = []
    fasta_iter = SeqIO.parse(in_file, file_type)
    for record in fasta_iter:
        target_seqs_and_names.append((record.id, record.seq.upper().transcribe()))
    return target_seqs_and_names


def find(string_to_analyze: str, char_to_find: str):
    """
    Finds all instances of a character char_to_find in a string string_to_analyze.
    """
    return [i for i, ltr in enumerate(string_to_analyze) if ltr == char_to_find]


def transcribe_seq_file(seq_file: str):
    # seq_file must be a .txt file
    with open(seq_file) as f:
        for i in f:
            out_seq = Seq(i).upper().transcribe()
    return out_seq


def find_repeat_targets(new_data, min_true_cov=0, fileout=False, file=''):
    # recall that new_data is a list of lists with one entry per target sequence where each entry is of the form:
    # [name, sequ, (ribozyme_sequences, IGS_and_guide_sequences, IGSes, guide_sequences, (og_idx, ref_idx))]
    big_repeats = []

    col = 0  # initialize the current column counter
    for org, sequ, cat_site_data in new_data:  # for each target sequence
        IGS_subset_a = list(set(cat_site_data[2]))  # extract each IGS sequence in this organism
        for j in range(col + 1, len(new_data)):  # this is to avoid repeat combinations
            # make a subset of second column of unique IGS values
            # and find the shared values between sets with no duplicates
            org_b, sequ_b, cat_site_data_b = new_data[j]
            no_dupes = (set(IGS_subset_a) & set(list(set(cat_site_data_b[2]))))

            # remove any nans
            repeats = [item for item in no_dupes if str(item) != 'nan']
            # append to our big list
            if repeats:
                big_repeats.append(repeats)
        col += 1  # move to the next column

    print('Found repeat subsets. Now analyzing sequences...')

    # now flatten the list to compare all IGS sequences against one another
    flat_repeats = [item for sublist in big_repeats for item in sublist]
    # remove duplicates of all IGSes found
    filtered_list = list(set(flat_repeats))

    to_optimize = defaultdict(list)  # will keep those IGSes that meet our minimum percent coverage later
    to_keep_single_targets = defaultdict(list)

    # Now find the IGS sequences in the original target sequences and extract matching data
    big_temp_list = []
    for IGS_sequ in filtered_list:  # for each IGS sequence
        temp_list = []
        coverage_count = 0  # this will track how many target sequences contain the IGS
        col = 0
        for org, sequ, cat_site_data in new_data:  # in each target sequence (target seq is org)
            pos = [i for i, e in enumerate(cat_site_data[2]) if e == IGS_sequ]  # extract positions of matching IGSes

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
                    ref_pos = cat_site_data[4][p][1]
                    if on_target_count[ref_pos]:
                        on_target_count[ref_pos] += 1
                    else:
                        on_target_count[ref_pos] = 1

                    guide_id = str(IGS_sequ) + str(
                        ref_pos)  # save a guide ID that is the IGS sequence and the reference position
                    temp_list.append([IGS_sequ, None, None, None, org, target_num, cat_site_data[4][p][0],
                                 ref_pos, cat_site_data[3][p], cat_site_data[1][p],
                                 cat_site_data[0][p], sequ, guide_id])
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
                to_optimize[item[12]].append(item)
            # Else if the true % coverage is above min_true_cov but only one sequence still keep but will not optimize
            elif true_coverage >= min_true_cov:
                to_keep_single_targets[item[12]].append(item)

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
        ranked_sorted_IGS.to_csv(f'{file}/Ranked Ribozyme Designs with Raw Guide Sequence Designs.csv',
                                 index=False)

        new_data_df = pd.DataFrame.from_records(new_data).T
        new_data_df.to_csv(f'{file}/All catalytic U data unsorted.csv', index=True)

    return big_temp_list, to_optimize, filtered_list, ranked_IGS, ranked_sorted_IGS, to_keep_single_targets


def find_repeat_targets_initialize_loop(igs_data_a, to_compare):
    for igs_data_b in to_compare:
        # make a subset of second column of unique IGS values
        # and find the shared values between sets with no duplicates
        no_dupes = igs_data_a & igs_data_b
        # remove any nans
        repeats = [item for item in no_dupes if str(item) != 'nan']
        if repeats:
            return repeats


def find_repeat_targets_middle_loop(igs_sequ, org, sequ, cat_site_data):
    # in each target sequence (target seq is org)
    pos = [i for i, e in enumerate(cat_site_data[2]) if e == igs_sequ]  # extract positions of matching IGSes
    if pos:  # if not in this particular seq, try next seq
        # if the IGS is found in the target seq, get all the data we want
        # To calculate whether a particular IGS is in the same position in multiple targets or not I am making a
        # HUGE assumption: as the sequences are already aligned to a single reference sequence, any base
        # position will have EXACTLY the same position numbering. This allows my code to use dictionaries to
        # determine whether something is or is not on target as the position must match exactly. If we want some
        # tolerance in how many base pairs away a position can be to be considered on-target I will need to find
        # another way to calculate this which would probably use a LOT more computational power.
        on_target_count = defaultdict(list)
        target_num = len(pos)

        for p in pos:
            ref_pos = cat_site_data[4][p][1]
            if on_target_count[ref_pos]:
                on_target_count[ref_pos] += 1
            else:
                on_target_count[ref_pos] = 1

            guide_id = str(igs_sequ) + str(
                ref_pos)  # save a guide ID that is the IGS sequence and the reference position
            temp_list = [igs_sequ, None, None, None, org, target_num, cat_site_data[4][p][0],
                         ref_pos, cat_site_data[3][p], cat_site_data[1][p],
                         cat_site_data[0][p], sequ, guide_id]
        return (temp_list, on_target_count)


def find_repeat_targets_last_loop(item, new_data, min_true_cov, perc_coverage, on_target_count, coverage_count):
    to_optimize = defaultdict(list)  # will keep those IGSes that meet our minimum percent coverage later
    to_keep_single_targets = defaultdict(list)

    item[1] = perc_coverage
    item[2] = on_target_count[item[7]] / coverage_count  # out of all organisms with this IGS at this position
    true_coverage = on_target_count[item[7]] / len(new_data)
    item[3] = true_coverage
    big_temp_list = list(item)
    # if the true % coverage is above min_true_cov, and more than one sequence, mark this IGS for optimization
    if true_coverage >= min_true_cov and on_target_count[item[7]] > 1:
        # save the IGS sequence AND the matching index in the reference sequence
        to_optimize[item[12]].append(item)
    # Else if the true % coverage is above min_true_cov but only one sequence still keep but will not optimize
    elif true_coverage >= min_true_cov:
        to_keep_single_targets[item[12]].append(item)

    return (big_temp_list, to_optimize, to_keep_single_targets)


def find_repeat_targets_loop(igs_sequ, new_data, min_true_cov):
    to_optimize = defaultdict(list)  # will keep those IGSes that meet our minimum percent coverage later
    to_keep_single_targets = defaultdict(list)
    temp_list = []
    big_temp_list = []
    coverage_count = 0  # this will track how many target sequences contain the IGS
    col = 0
    for org, sequ, cat_site_data in new_data:  # in each target sequence (target seq is org)
        pos = [i for i, e in enumerate(cat_site_data[2]) if e == igs_sequ]  # extract positions of matching IGSes

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
                ref_pos = cat_site_data[4][p][1]
                if on_target_count[ref_pos]:
                    on_target_count[ref_pos] += 1
                else:
                    on_target_count[ref_pos] = 1

                guide_id = str(igs_sequ) + str(
                    ref_pos)  # save a guide ID that is the IGS sequence and the reference position
                temp_list.append([igs_sequ, None, None, None, org, target_num, cat_site_data[4][p][0],
                                  ref_pos, cat_site_data[3][p], cat_site_data[1][p],
                                  cat_site_data[0][p], sequ, guide_id])
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
            to_optimize[item[12]].append(item)
        # Else if the true % coverage is above min_true_cov but only one sequence still keep but will not optimize
        elif true_coverage >= min_true_cov:
            to_keep_single_targets[item[12]].append(item)

    return (big_temp_list, to_optimize, to_keep_single_targets)


def optimize_sequences(to_optimize, thresh, guide_length: int, ribo_seq, single_targets, fileout=False, file='',
                       score_type='quantitative',
                       gaps_allowed=True):
    print(f'Optimizing {len(to_optimize)} guide sequences...')

    in_data = [(key, to_optimize[key], thresh, score_type, gaps_allowed, guide_length, ribo_seq) for key in to_optimize]

    with Pool() as pool:
        opti_seqs = pool.starmap(optimize_sequences_loop, in_data)

    # If there are any single targets to keep (depends on min true coverage setting) add them here
    if single_targets:
        print(f'Storing {len(single_targets)} single target sequences. If you do not want single target guides, '
              f'please increase your min true coverage parameter.')
        for key in single_targets.keys():
            guide = single_targets[key][0][8]
            design_sequence = guide + 'G' + re.sub(r'[0-9]+', '', key)
            ribo_design = design_sequence + ribo_seq

            # set score to nan
            opti_seqs.append([single_targets[key][0][0], single_targets[key][0][7], 1, single_targets[key][0][1],
                              single_targets[key][0][2], single_targets[key][0][3],
                              [(target[4], target[6], target[5] - 1) for target in single_targets[key]], guide,
                              design_sequence, ribo_design])

    if fileout:
        sorted_opti_seqs = pd.DataFrame(data=opti_seqs, index=None, columns=['IGS', 'Reference index', 'Score', '% cov',
                                                                             '% on target', 'True % cov',
                                                                             '(Target name, Target idx, Other '
                                                                             'occurrences of IGS in target sequence)',
                                                                             'Optimized guide',
                                                                             'Optimized guide + G + IGS',
                                                                             'Full Ribozyme design'],
                                        dtype=object).sort_values(by=['True % cov', 'Score'], ascending=[False, False])
        sorted_opti_seqs.to_csv(
            f'{file}/Ranked Ribozyme Designs with Optimized Guide Sequence Designs {score_type}.csv', index=False)

    print('All guide sequences optimized.\n')
    return opti_seqs


def optimize_sequences_loop(name, items, thresh, score_type, gaps_allowed, guide_length, ribo_seq):
    guides_to_optimize = [target[8] for target in items]

    # do a MSA of the sequences and optimize the guide sequence with this MSA
    opti_seq, score = msa_and_optimize(name, guides_to_optimize, thresh, score_type, gaps_allowed)

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


def msa_and_optimize(name, seqs_to_align, thresh=0.7, score_type='quantitative', gaps_allowed='True'):
    # seqs_to_align is a list containing a sequence from an individual organism in each position.

    with open(f'to_align_{name}.fasta', 'w') as f:
        for line in range(len(seqs_to_align)):
            f.write('>seq' + str(line) + '\n' + str(seqs_to_align[line]) + '\n')

    # now align the sequences.
    muscle_exe = 'muscle5'

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


def get_quantitative_score(msa, thresh=0.7, chars_to_ignore=None, count_gaps=True, penalize_trailing_gaps=False):
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


def calc_shannon_entropy(target_names_and_seqs, ref_name_and_seq, base=None, count_gap=False, fileout=False, file=''):
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
                elif nucleotide == 'R':
                    probs[ref_idx]['A'] += 0.5
                    probs[ref_idx]['G'] += 0.5
                elif nucleotide == 'Y':
                    probs[ref_idx]['U'] += 0.5
                    probs[ref_idx]['C'] += 0.5
                elif nucleotide == 'M':
                    probs[ref_idx]['A'] += 0.5
                    probs[ref_idx]['C'] += 0.5
                elif nucleotide == 'K':
                    probs[ref_idx]['G'] += 0.5
                    probs[ref_idx]['U'] += 0.5
                elif nucleotide == 'S':
                    probs[ref_idx]['G'] += 0.5
                    probs[ref_idx]['C'] += 0.5
                elif nucleotide == 'W':
                    probs[ref_idx]['A'] += 0.5
                    probs[ref_idx]['U'] += 0.5
                elif nucleotide == 'H':
                    probs[ref_idx]['A'] += 1 / 3
                    probs[ref_idx]['C'] += 1 / 3
                    probs[ref_idx]['U'] += 1 / 3
                elif nucleotide == 'B':
                    probs[ref_idx]['C'] += 1 / 3
                    probs[ref_idx]['G'] += 1 / 3
                    probs[ref_idx]['U'] += 1 / 3
                elif nucleotide == 'D':
                    probs[ref_idx]['A'] += 1 / 3
                    probs[ref_idx]['G'] += 1 / 3
                    probs[ref_idx]['U'] += 1 / 3
                elif nucleotide == 'V':
                    probs[ref_idx]['A'] += 1 / 3
                    probs[ref_idx]['C'] += 1 / 3
                    probs[ref_idx]['G'] += 1 / 3
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
