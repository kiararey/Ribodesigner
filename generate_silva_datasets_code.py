import os
import numpy as np
np.random.seed(1)
# This is a very poorly written script, but it does the job and I am too tired to make it more efficient.

# Determine at what taxonomic level we want the datasets. This must match one of the filenames from SequencePrepper.py:
taxonomy_level = 'Genus'

# Here is the file where the dataset is. I've already run a script on this to separate it by levels of taxonomy.
# To do this, SequencePrepper.py on the SILVA database file.
silva_by_taxonomy_path = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/{taxonomy_level}'

# Here is the file where we want to save the dataset
output_path = 'Datasets_used/SILVA_squished_datasets'

# There are only Bacteria in this dataset, so will only do that:
fasta_file_names = np.array([file_name for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta' in file_name])

# Prepare list to store data in:
num_of_datasets = 5
seqs_to_write = [np.empty(0)]*5


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


def read_fasta_file_full_name(filename):
    target_seqs = []
    with open(filename) as f:
        for line in f:  # read each line
            if line[0] == '>':  # when we reach the ID of a line
                name = line
            else:  # if we are not at an ID line we must be in the next sequ
                target_seqs.append(SilvaSequence(name, line))
    return target_seqs

old_perc = 0
for i, genus in enumerate(fasta_file_names):
    # Pick five unique organisms per genus that are not the same species
    target_seqs_and_names = read_fasta_file_full_name(f'{silva_by_taxonomy_path}/{genus}')
    if len(target_seqs_and_names) <= 5:
        # if there are only five or less sequences, all of the sequences will be in all datasets.
        for row in range(len(seqs_to_write)):
            seqs_to_write[row] = np.append(seqs_to_write[row], target_seqs_and_names)
        continue
    # Now pick five unique organisms per genus five times
    j = 0
    idx_used = []  # sequences we already used, ideally we won't reuse any of these but
    while j < num_of_datasets:
        # take 5 sequences randomly
        lines_to_take = np.random.randint(0, len(target_seqs_and_names), 5)
        sequences_taken = [target_seqs_and_names[num] for num in lines_to_take]
        # check that all sequences are of unique species
        species_of_sequences_taken_set = set(sequence.give_taxonomy('Species') for sequence in sequences_taken)

        if len(species_of_sequences_taken_set) == 5:
            seqs_to_write[j] = np.append(seqs_to_write[j], sequences_taken)
            j += 1
            continue
        else:
            total_species_set = set(sequence.give_taxonomy('Species') for sequence in target_seqs_and_names)
            if len(total_species_set) < 5:
                # if there is no way to get 5 unique species, just take the L:
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
        num_needed = 5 - len(sequences_taken)
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
                                                        if sequence.give_taxonomy('Species') in species_not_represented]
                if len(sequences_of_species_not_represented) <= num_needed:
                    candidate_sequences.extend(sequences_of_species_not_represented)
                    break
                new_num_needed = 5 - len(sequences_taken) - len(candidate_sequences)
                new_lines_to_take = np.random.randint(0, len(sequences_of_species_not_represented), new_num_needed)

                candidate_sequences.extend([sequences_of_species_not_represented[num] for num in new_lines_to_take])

                remaining_candidates = [sequences_of_species_not_represented[num] for num in
                                        range(len(candidate_sequences)) if num not in new_lines_to_take]

                species_of_candidate_sequences = [sequence.give_taxonomy('Species') for sequence in candidate_sequences]
        sequences_taken.extend(candidate_sequences)
        seqs_to_write[j] = np.append(seqs_to_write[j], sequences_taken)
        j += 1
    new_perc = int(i / len(fasta_file_names) * 100)
    if new_perc > old_perc and new_perc % 5 == 0:
        print(f'{new_perc}% of {len(fasta_file_names)} datasets completed...')
    old_perc = new_perc

print(f'100% of {len(fasta_file_names)} datasets completed. Now saving in {output_path[:-1]}...')

# Add to FASTA file and save
for i, dataset in enumerate(seqs_to_write):
    file_name = f'{output_path}Bacteria_by_{taxonomy_level}_{i + 1}.fasta'
    with open(file_name, 'w') as f:
        for sequence in set(dataset):
            f.writelines([sequence.id, sequence.seq, '\n'])

print('Files saved!')
