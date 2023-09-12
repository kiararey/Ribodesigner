import os
import numpy as np
from numpy.random import default_rng


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
        if level not in level_names:
            print('Taxonomic level not found')
        # gives us the name at a certain taxonomic level
        all_levels = self.id.split(';')
        try:
            here_is_taxonomy = all_levels[level_names[level]]
        except IndexError:
            here_is_taxonomy = ''
        return here_is_taxonomy


def read_fasta_file_full_name(filename, exclude=[], include=[], exclude_include_taxonomy_level=''):
    # Include and exclude are mutually exclusive, and will take into account exclude first if both are given.
    target_seqs = []
    with open(filename) as f:
        for line in f:  # read each line
            if line[0] == '>':  # when we reach the ID of a line
                name = line
            else:  # if we are not at an ID line we must be in the next sequ
                if exclude and exclude_include_taxonomy_level:
                    putative_sequence = SilvaSequence(name, line)
                    if putative_sequence.give_taxonomy(level=exclude_include_taxonomy_level) not in exclude:
                        target_seqs.append(putative_sequence)
                elif include and exclude_include_taxonomy_level:
                    putative_sequence = SilvaSequence(name, line)
                    if putative_sequence.give_taxonomy(level=exclude_include_taxonomy_level) in include:
                        target_seqs.append(putative_sequence)
                else:
                    target_seqs.append(SilvaSequence(name, line))
    return target_seqs


def generate_silva_datasets(silva_by_taxonomy_path: str, output_path: str, num_of_datasets: int = 5,
                            num_of_sequences: int = 5, exclude_only: list = [], include_only: list = [],
                            exclude_taxonomy_level: str = '', seed: int = 1):
    """
    Will generate a squished dataset given an input file that has already been made through SequencePrepper.py. Will
    return num_of_datasets datasets containing num_of_sequences in each without any repeated species.
    :param silva_by_taxonomy_path:  Here is the file where the dataset to squish is. To make this, run
    SequencePrepper.py on the Silva database file or similarly structured data to separate it by levels of taxonomy.
    :param output_path: Here is the file where we want to save the dataset
    :param num_of_datasets: This is how many total squished datasets we want
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

    if (exclude_only or include_only) and exclude_taxonomy_level not in ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Taxon']:
        taxonomy_choose = {0: 'Domain', 1: 'Phylum', 2: 'Class', 3: 'Order',
                           4: 'Family', 5: 'Genus', 6: 'Species', 7: 'Taxon'}
        ini = -1
        while ini < 0 or ini > 7:
            ini = int(input(f'Taxonomy level not recognized. Please choose a number to assign:\n{taxonomy_choose}'))
        exclude_taxonomy_level = taxonomy_choose[ini]

    # From the given path extract only the .fasta files
    fasta_file_names = np.array(
        [file_name for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta' in file_name])

    # Prepare list to store data in:
    seqs_to_write = [np.empty(0)] * num_of_datasets

    old_perc = 0
    perc_under_seq_num = 0
    seqs_kept = 0
    for i, seq_file in enumerate(fasta_file_names):
        # Pick five unique organisms per genus that are not the same species
        target_seqs_and_names = read_fasta_file_full_name(f'{silva_by_taxonomy_path}/{seq_file}', exclude=exclude_only,
                                                          include=include_only,
                                                          exclude_include_taxonomy_level=exclude_taxonomy_level)
        if not target_seqs_and_names:
            continue
        seqs_kept += 1
        if len(target_seqs_and_names) <= num_of_sequences:
            # if there are only num_of_sequences or less sequences, all of the sequences will be in all datasets.
            for row in range(len(seqs_to_write)):
                seqs_to_write[row] = np.append(seqs_to_write[row], target_seqs_and_names)
            perc_under_seq_num += 1
            continue
        # Now pick num_of_sequences unique organisms per taxonomy level num_of_datasets times
        j = 0
        idx_used = []  # sequences we already used, ideally we won't reuse any of these but
        while j < num_of_datasets:
            # take num_of_sequences sequences randomly
            array = np.array(target_seqs_and_names, dtype=tuple)
            rng = default_rng(seed=seed)
            sequences_taken = rng.choice(array, num_of_sequences, replace=False)

            # check that all sequences are of unique species
            species_of_sequences_taken_set = set(sequence.give_taxonomy('Species') for sequence in sequences_taken)

            if len(species_of_sequences_taken_set) == num_of_sequences:
                seqs_to_write[j] = np.append(seqs_to_write[j], sequences_taken)
                j += 1
                continue
            else:
                total_species_set = set(sequence.give_taxonomy('Species') for sequence in target_seqs_and_names)
                if len(total_species_set) < num_of_sequences:
                    # if there is no way to get num_of_sequences unique species, just take the L:
                    seqs_to_write[j] = np.append(seqs_to_write[j], sequences_taken)
                    j += 1
                    continue
            # Make a list containing species names that are not represented in the current pool
            species_not_represented = [species_name for species_name in total_species_set
                                       if species_name not in species_of_sequences_taken_set]

            sequences_of_species_not_represented = [sequence for sequence in target_seqs_and_names
                                                    if sequence.give_taxonomy('Species') in species_not_represented]

            # Find repeated species, remove one until we have enough unique members:
            def get_unique_members(list_of_sequences, set_of_species_represented, level):
                for current_spot, sequence in enumerate(list_of_sequences):
                    sequences_found = sequence.give_taxonomy(level)
                    if len(sequences_found) != len(set(sequences_found)):
                        list_of_sequences = [seq for spot, seq in enumerate(list_of_sequences) if spot != current_spot]
                        if len(list_of_sequences) == len(set_of_species_represented):
                            break
                return list_of_sequences

            sequences_taken = get_unique_members(sequences_taken.copy(), species_of_sequences_taken_set, 'Species')

            # Finally, find enough sequences that are not the same species to fill out our current dataset
            num_needed = num_of_sequences - len(sequences_taken)
            new_lines_to_take = np.random.randint(0, len(sequences_of_species_not_represented), num_needed)
            candidate_sequences = [sequences_of_species_not_represented[num] for num in new_lines_to_take]
            remaining_candidates = [sequences_of_species_not_represented[num] for num in range(len(candidate_sequences))
                                    if num not in new_lines_to_take]
            species_of_candidate_sequences = [sequence.give_taxonomy('Species') for sequence in candidate_sequences]
            # If we only need one species or there is no way to get all unique species just take as many needed
            if not num_needed == 1 | len(set(species_not_represented)) < num_needed:
                # first check if the sequences we got randomly are not the same species
                tries = 0
                while len(species_of_candidate_sequences) > len(set(species_of_candidate_sequences)):
                    tries += 1
                    candidate_sequences = get_unique_members(candidate_sequences.copy(),
                                                             set(species_of_candidate_sequences), 'Species')
                    species_of_candidate_sequences = set([sequence.give_taxonomy('Species')
                                                          for sequence in candidate_sequences])
                    all_species_represented = species_of_sequences_taken_set.union(species_of_candidate_sequences)
                    species_not_represented = [species.give_taxonomy('Species') for species in remaining_candidates
                                               if species.give_taxonomy('Species') not in all_species_represented]
                    sequences_of_species_not_represented = [sequence for sequence in remaining_candidates
                                                            if sequence.give_taxonomy(
                            'Species') in species_not_represented]
                    if len(sequences_of_species_not_represented) <= num_needed:
                        candidate_sequences.extend(sequences_of_species_not_represented)
                        break
                    new_num_needed = num_of_sequences - len(sequences_taken) - len(candidate_sequences)
                    new_lines_to_take = np.random.randint(0, len(sequences_of_species_not_represented), new_num_needed)

                    candidate_sequences.extend([sequences_of_species_not_represented[num] for num in new_lines_to_take])

                    remaining_candidates = [sequences_of_species_not_represented[num] for num in
                                            range(len(candidate_sequences)) if num not in new_lines_to_take]

                    species_of_candidate_sequences = [sequence.give_taxonomy('Species') for sequence in
                                                      candidate_sequences]
            sequences_taken.extend(candidate_sequences)
            seqs_to_write[j] = np.append(seqs_to_write[j], sequences_taken)
            j += 1
        new_perc = int(i / len(fasta_file_names) * 100)
        if new_perc > old_perc and new_perc % 5 == 0:
            print(f'{new_perc}% of {len(fasta_file_names)} datasets completed...')
        old_perc = new_perc

    print(f'100% of {len(fasta_file_names)} datasets completed.\n{seqs_kept} datasets met inclusion criteria. '
          f'{round(perc_under_seq_num/seqs_kept*100, 2)}% of kept datasets had less than the requested number of '
          f'sequences.\nNow saving in {output_path[:-1]}...')

    # Add to FASTA file and save
    taxonomy_level = silva_by_taxonomy_path.split('/')[-1]
    kingdom_level = silva_by_taxonomy_path.split('/')[1].split('_')[7]
    for i, dataset in enumerate(seqs_to_write):
        file_name = f'{output_path}/{kingdom_level}_by_{taxonomy_level}_{i + 1}.fasta'
        with open(file_name, 'w') as f:
            for sequence in set(dataset):
                f.writelines([sequence.id, sequence.seq, '\n'])

    print('Files saved!\n')
