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


# Expanded version of SilvaSequence
class TargetSeq:
    # Important: TargetSeq id must be formatted like a Silva sequence if we want the give_taxonomy function to work!

    def __init__(self, id: str, seq: Seq):
        self.id = id
        self.seq = seq
        self.cat_sites = []

    def __str__(self):
        return f'{self.id}({self.seq})'

    def __repr__(self):
        return f'{self.id}({self.seq})'

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

    def find_cat_sites(self, igs_length: int = 5, guide_length: int = 50,
                       min_length: int = 35):
        """
        Finds all instances of a U or T in a set of sequences and creates ribozyme designs for these sites.

        :param igs_length: desired IGS sequence length.
        :param guide_length: desired guide binding sequence length.
        :param min_length: minimum guide sequence length from 3' end. Must be smaller than guide_length.
        :return:
        """
        if self.cat_sites:
            if len(self.cat_sites[0]) < 4:
                print('TargetSeq object already has catalytic sites.')
            else:
                print('TargetSeq object was already aligned and has catalytic sites.')
            return

        #   (a.k.a. minimum guide sequence length from 3' end) ex: if you want your guide sequence to
        #   bind to at least 35 nt at the 3' end of the target sequence, set min_length = 35.

        # run function to find the index of all U residues of each sequence in target_seqs
        # find all possible splice sites
        idx = [i for i, ltr in enumerate(self.seq) if ltr == 'T']
        # remove indexes that are < guide_length or > len(sequ) - igs_length (must have enough residues to attach to)
        idx_new = [res for res in idx if igs_length <= res < (len(self.seq) - min_length)]

        if not idx_new:
            print(f'No viable catalytic sites in {self.id}')
            self.cat_sites = None
        else:
            data = []
            for i in idx_new:
                # generate complementary guide sequence guide_length residues *downstream* of U site
                guide = self.seq[i + 1:i + guide_length + 1].reverse_complement()
                # generate complementary IGS sequence igs_length bases long *upstream* of U site
                igs = self.seq[i - igs_length:i].reverse_complement()
                # Store as (idx, igs, guide)
                data.append((i + 1, igs, guide))  # we're adding 1 to the idx because of 0 based indexing
            self.cat_sites = data

    def align_to_ref(self, ref_name_and_seq, igs_length: int = 5, guide_length: int = 50, min_length: int = 35):
        """
        :param ref_name_and_seq: TargetSequence object
        :return:
        """
        if not self.cat_sites:
            self.find_cat_sites(igs_length, guide_length, min_length)
        elif len(self.cat_sites[0]) == 4:
            print('TargetSeq object has already been aligned to a reference sequence.')
            return

        # will have to keep in mind the potential lengths of the sequences and add length igs_length to our final
        # E.coli index
        aligner = PairwiseAligner(mode='global')
        alignments = aligner.align(self.seq, ref_name_and_seq[1])[0]

        # seq_a is the test sequence, seq_b is the reference sequence
        seq_a, seq_b = alignments
        # seq_a is the test sequence, seq_b is the reference sequence

        # obtain index of new U
        idx_seq_a = [i for i, ltr in enumerate(seq_a) if ltr == 'T']

        data = []
        current_spot = 0
        for og_idx, igs, guide in self.cat_sites:
            seq_a_idx = len(seq_a[:idx_seq_a[current_spot]].replace('-', '')) + 1
            while seq_a_idx != og_idx:
                current_spot += 1
                seq_a_idx = len(seq_a[:idx_seq_a[current_spot]].replace('-', '')) + 1
            # find what index that is based on the reference sequence
            ref_string = seq_b[:idx_seq_a[current_spot]]
            ref_idx = len(ref_string.replace('-', '')) + 1  # turns zero based indexing to ref_seq numbering

            data.append((og_idx, ref_idx, igs, guide))

        self.cat_sites = data