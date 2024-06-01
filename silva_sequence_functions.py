import os
import numpy as np
from numpy.random import default_rng
from alive_progress import alive_bar
from collections import defaultdict
import math


class SilvaSequence:
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

    def __str__(self):
        return f'{self.id}({self.seq})'

    def __repr__(self):
        return f'{self.id}({self.seq})'

    def give_taxonomy(self, level):
        level_names = {'Domain': 0, 'Phylum': 1, 'Class': 2, 'Order': 3,
                       'Family': 4, 'Genus': 5, 'Species': 6, 'Taxon': 7}
        # Eukaryotic levels:
        level_names_euk = {'Domain': 0, 'Kingdom': 1, 'Subkingdom': 2, 'Phylum': 3, 'Subphylum': 4, 'Class': 5,
                           'Order': 6, 'Family': 7, 'Genus': 8, 'Species': -1}
        if level == 'any':
            return self.id[1:]
        if level not in level_names or level not in level_names_euk:
            print('Taxonomic level not found')
        # gives us the name at a certain taxonomic level
        all_levels = self.id.split(';')
        all_levels[0] = all_levels[0].split(' ')[-1]
        if 'Eukaryota' not in all_levels[0]:
            try:
                here_is_taxonomy = all_levels[level_names[level]]
            except IndexError:
                here_is_taxonomy = ''
        else:
            try:
                here_is_taxonomy = all_levels[level_names_euk[level]]
            except IndexError:
                here_is_taxonomy = ''
        return here_is_taxonomy


def selection_routine(target_seqs, putative_sequence, exclude, include, exclude_include_taxonomy_level):
    if exclude_include_taxonomy_level == 'any':
        the_level_to_check_name = putative_sequence.id
    else:
        the_level_to_check_name = putative_sequence.give_taxonomy(level=exclude_include_taxonomy_level)
    # First check if we should exclude
    for check in exclude:
        if check in the_level_to_check_name:
            continue
    # If not excluded, check if we should include
    if include:
        for check in include:
            if check in the_level_to_check_name:
                target_seqs.append(putative_sequence)
                continue
    else:
        target_seqs.append(putative_sequence)
def read_fasta_file_full_name(filename, exclude=None, include=None, exclude_include_taxonomy_level='any'):
    # Include and exclude are mutually exclusive, and will take into account exclude first if both are given.
    if include is None:
        include = []
    if exclude is None:
        exclude = []
    target_seqs = []
    started_seq = False
    with open(filename) as f:
        for line in f:  # read each line
            if line[0] == '>':  # when we reach the ID of a line
                if started_seq:
                    putative_sequence = SilvaSequence(name, sequ)
                    selection_routine(target_seqs, putative_sequence, exclude, include, exclude_include_taxonomy_level)
                name = line
                sequ = ''
                started_seq = True
            else:  # if we are not at an ID line we must be in the next sequ
                sequ += line.strip('\n')
        putative_sequence = SilvaSequence(name, sequ)
        selection_routine(target_seqs, putative_sequence, exclude, include, exclude_include_taxonomy_level)
    return target_seqs


def get_unique_members(list_of_sequences, set_of_species_represented, level):
    for current_spot, sequence in enumerate(list_of_sequences):
        sequences_found = sequence.give_taxonomy(level)
        if len(sequences_found) != len(set(sequences_found)):
            list_of_sequences = [seq for spot, seq in enumerate(list_of_sequences) if spot != current_spot]
            if len(list_of_sequences) == len(set_of_species_represented):
                break
    return list_of_sequences
def generate_silva_datasets(silva_by_taxonomy_path: str, output_path: str, num_of_sequences: int = 5, divide_by='Order',
                            unique_at='Genus', exclude_only: list = None, include_only: list = None,
                            exclude_taxonomy_level: str = 'any', seed: int = 1, pick_from_file: bool = False,
                            repeat_species=False, only_test = False):
    """
    Will generate a squished dataset given an input file that has already been made through SequencePrepper.py. Will
    return test and target datasets containing ideally num_of_sequences in each without any repeated species.
    :param silva_by_taxonomy_path:  Here is the file where the dataset to squish is. To make this, run
    SequencePrepper.py on the Silva database file or similarly structured data to separate it by levels of taxonomy.
    :param output_path: Here is the file where we want to save the dataset
    :param num_of_sequences: This is the maximum of how many sequences per taxonomy file we want to keep
    :param exclude_only: this is a list for sequences we want to exclude at the taxonomic level given by
    exclude_taxonomy_level. Mutually exclusive with include_only. Note that if both exclude_only and include_only are
    given exclude_only will take preference.
    :param include_only: this is a list for sequences we want to include at the taxonomic level given by
    exclude_taxonomy_level. Mutually exclusive with exclude_only. Note that if both exclude_only and include_only are
    given exclude_only will take preference.
    :param exclude_taxonomy_level: the taxonomic level that we want to exclude or include things.
    Must be one of the following: ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Taxon']
    :param seed: random seed, for consistency across runs to debug.
    :return:
    """

    if type(include_only) == str:
        include_only = [include_only]
    if type(exclude_only) == str:
        exclude_only = [exclude_only]

    # From the given path extract only the .fasta files
    if silva_by_taxonomy_path.endswith('.fasta'):
        fasta_file_names = [silva_by_taxonomy_path]
    elif not pick_from_file:
        fasta_file_names = np.array(
            [f'{silva_by_taxonomy_path}/{file_name}' for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta'
             in file_name])
    else:
        fasta_file_names = [silva_by_taxonomy_path]

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Prepare list to store data in:
    target_seqs_and_names = []
    print('Importing sequences...')
    with alive_bar(len(fasta_file_names), spinner='fishes') as bar:
        for i, seq_file in enumerate(fasta_file_names):
            # Pick five unique organisms per genus that are not the same species
            target_seqs_and_names.extend(read_fasta_file_full_name(seq_file, exclude=exclude_only, include=include_only,
                                                                   exclude_include_taxonomy_level=exclude_taxonomy_level))
            bar()

    print('Making datasets...')
    with alive_bar(unknown='fish', spinner='fishes') as bar:
        if repeat_species:
            array = np.array(target_seqs_and_names, dtype=tuple)
            rng = default_rng(seed=seed)
            if num_of_sequences > len(target_seqs_and_names):
                print(f'There are not enough sequences available. Defaulting to {len(target_seqs_and_names)} sequences')
                num_of_sequences = len(target_seqs_and_names)
            test_seqs = rng.choice(array, num_of_sequences, replace=False)
            if not only_test:
                target_seqs = rng.choice(array, num_of_sequences, replace=False)
        else:
            # Get the species of everyone so that we can make sure we can have all species represented
            species_represented_dict = defaultdict(lambda: defaultdict(lambda: []))
            species_per_order = defaultdict(lambda: set())
            for seq in target_seqs_and_names:
                species_represented_dict[seq.give_taxonomy(divide_by)][seq.give_taxonomy(unique_at)].append(seq)
                species_per_order[seq.give_taxonomy(divide_by)].add(seq.give_taxonomy(unique_at))

            # Sort by species, see if there is enough for both test and target. If there is not enough then add all but
            # one into target. If there is only one sequence add to test
            test_seqs = []
            target_seqs = []
            rng = default_rng(seed=seed)
            targeting_sets_meeting_criteria = 0
            testing_sets_meeting_criteria = 0
            orders_in_targeting_set = 0
            orders_in_testing_set = 0
            for key, items in species_represented_dict.items():
                to_pick = list(species_per_order[key])
                if len(items) >= num_of_sequences:
                    # Sequences between target and test don't have to be different for most user cases.
                    sequences_to_take_target = rng.choice(to_pick, num_of_sequences, replace=False)
                    sequences_to_take_test = rng.choice(to_pick, num_of_sequences, replace=False)
                    targeting_sets_meeting_criteria += 1
                    testing_sets_meeting_criteria += 1
                    orders_in_targeting_set += 1
                    orders_in_testing_set += 1
                else:
                    sequences_to_take_target = rng.choice(to_pick, len(items), replace=False)
                    sequences_to_take_test = rng.choice(to_pick, len(items), replace=False)
                    orders_in_targeting_set += 1
                    orders_in_testing_set += 1
                for col in range(len(sequences_to_take_target)):
                    if not only_test:
                        target_seqs.extend(rng.choice(items[sequences_to_take_target[col]], 1, replace=False))
                    test_seqs.extend(rng.choice(items[sequences_to_take_test[col]], 1, replace=False))
    # get plural:
    plural_dict = {'Domain': 'domains', 'Phylum': 'phyla', 'Class': 'classes', 'Order': 'orders', 'Family': 'families',
                   'Genus': 'genera', 'Species': 'species', 'Taxon': 'taxa', 'Kingdom': 'kingdoms',
                   'Subkingdom': 'subkingdoms', 'Subphylum': 'subphyla'}
    if not repeat_species:
        if not only_test:
            print(f'\n{targeting_sets_meeting_criteria} out of {len(species_represented_dict)} '
                  f'{plural_dict[divide_by]} had at least {num_of_sequences} unique {plural_dict[unique_at]} available.'
                  f'\nThere are {orders_in_targeting_set} {plural_dict[divide_by]} in targeting set and '
                  f'{orders_in_testing_set} {plural_dict[divide_by]} in testing set represented.'
                  f'\nTarget dataset has {len(target_seqs)} sequences and test dataset has {len(test_seqs)} sequences '
                  f'of unique {plural_dict[unique_at]}.\nNow saving in {output_path}...')
        else:
            print(f'\n{testing_sets_meeting_criteria} out of {len(species_represented_dict)} '
                  f'{plural_dict[divide_by]} had at least {num_of_sequences} unique {plural_dict[unique_at]} available.'
                  f'\nThere are {orders_in_testing_set} {plural_dict[divide_by]} in test set represented.'
                  f'\nTest dataset has {len(test_seqs)} sequences of unique {plural_dict[unique_at]}.'
                  f'\nNow saving in {output_path}...')
    else:
        if not only_test:
            print(f'Target dataset has {len(target_seqs)} sequences and test dataset has {len(test_seqs)} sequences. '
                  f'\nNow saving in {output_path}...')
        else:
            print(f'Test dataset has {len(test_seqs)} sequences.\nNow saving in {output_path}...')

    # Add to FASTA file and save
    print('Generating dataset files...')
    if include_only:
        kingdom_level = '_'.join(include_only) + '_Only'
    elif exclude_only:
        kingdom_level = 'All_but_' + '_'.join(exclude_only)
    else:
        kingdom_level = 'All'
    if only_test:
        save_this = [('test', test_seqs)]
    else:
        save_this = [('test', test_seqs), ('target', target_seqs)]
    for name, dataset in save_this:
        if repeat_species:
            file_name = f'{output_path}/{kingdom_level}_by_{divide_by}_{name}.fasta'
        else:
            file_name = f'{output_path}/{kingdom_level}_by_{divide_by}_unique_{unique_at}_{name}.fasta'
        with open(file_name, 'w') as f:
            for sequence in set(dataset):
                f.writelines([sequence.id, sequence.seq, '\n'])

    print('Files saved!\n#######################################################\n')
    return

