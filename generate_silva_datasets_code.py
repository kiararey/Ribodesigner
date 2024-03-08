import os
import numpy as np

from silva_sequence_functions import generate_silva_datasets
seed = 1

########################################################
# To generate SILVA squished datasets
# Determine at what taxonomic level we want the datasets. This must match one of the filenames from SequencePrepper.py:
taxonomy_level = 'Genus'

# Here is the file where the dataset is. I've already run a script on this to separate it by levels of taxonomy.
# To do this, SequencePrepper.py on the SILVA database file.
# silva_by_taxonomy_path = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_All_Kingdoms/{taxonomy_level}'
# silva_by_taxonomy_bac = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/{taxonomy_level}'
# silva_by_taxonomy_euk = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Eukaryota_Only/{taxonomy_level}'
# silva_by_taxonomy_arc = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Archaea_Only/{taxonomy_level}'
#
# # Here is the file where we want to save the dataset
# output_path = 'Datasets_used/SILVA_squished_datasets'
#
# # Prepare list to store data in:
# num_of_datasets = 5
# num_of_sequences = 1
#
# # Generate the squished datasets
# generate_silva_datasets(silva_by_taxonomy_bac, output_path, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences, seed=seed)
# generate_silva_datasets(silva_by_taxonomy_euk, output_path, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences, seed=seed)
# generate_silva_datasets(silva_by_taxonomy_arc, output_path, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences, seed=seed)
# generate_silva_datasets(silva_by_taxonomy_path, output_path, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences, seed=seed)


########################################################
# To generate datasets for targeted designs by order level test

# Here is the file where the dataset is. I've already run a script on this to separate it by levels of taxonomy.
# To do this, SequencePrepper.py on the SILVA database file.
silva_by_taxonomy_path = 'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Genus'

# Here is the file where we want to save the dataset
output_path_entero = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Enterobacterales_only_squished'
output_path_pseudo = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Pseudomonadales_only_squished'
output_path_both = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Pseudo_and_entero_only_squished'
output_path_background_no_e = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished_no_entero'
output_path_background_no_p = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished_no_pseudo'
output_path_background_no_p_or_e = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished'
output_path_background_gram_pos = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Gram_positives_only'
output_path_background_no_gram_pos = 'Datasets_used/SILVA_squished_datasets_1_per_genus/No_Gram_positives'

# There are only Bacteria in this dataset, so will only do that:
fasta_file_names = np.array([file_name for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta' in file_name])

# Prepare list to store data in:
num_of_datasets = 5
num_of_sequences_per_genus = 1

include_only_entero = ['Enterobacterales']
include_only_pseudo = ['Pseudomonadales']
include_only_both = ['Enterobacterales', 'Pseudomonadales']
taxonomy_level_of_inclusion = 'Order'
include_only_gram_pos = ['Actinobacteriota', 'Firmicutes']

# Generate the squished datasets
# print('entero only')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_entero, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_entero,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# print('pseudo only')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_pseudo, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_pseudo,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# print('entero or pseudo')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_both, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_both,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# # Generate background datasets (exclude one or the other or both)
# print('all but entero')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_e, num_of_datasets=1,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_entero,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# print('all but pseudo')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_p, num_of_datasets=1,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_pseudo,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# print('all but pseudo or entero')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_p_or_e, num_of_datasets=2,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_both,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# print('gram positives only')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_gram_pos, num_of_datasets=2,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_gram_pos,
#                         exclude_taxonomy_level='Phylum', seed=seed)
#
#
# print('all but gram positives')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_gram_pos, num_of_datasets=2,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_gram_pos,
#                         exclude_taxonomy_level='Phylum', seed=seed)
taxonomy_levels_all = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
to_generate = {'Phylum': ['Proteobacteria', 'Firmicutes'],
               'Class': ['Gammaproteobacteria', 'Bacilli'],
               'Order': ['Enterobacterales', 'Pseudomonadales', 'Bacillales'],
               'Family': ['Enterobacteriaceae', 'Pseudomonadaceae', 'Bacillaceae'],
               'Genus': ['Escherichia-Shigella', 'Pseudomonas', 'Bacillus']}
for taxonomy in taxonomy_levels_all:
    print(f'\nNow generating data for {taxonomy}:')

    for include in to_generate[taxonomy]:
        print(f'\nDataset is now {include}:')
        output_path = f'Datasets_used/SILVA_squished_datasets_1_per_genus/Selective datasets per taxonomy/{taxonomy}_{include}'
        generate_silva_datasets(silva_by_taxonomy_path, output_path + '_included', num_of_datasets=2,
                                num_of_sequences=num_of_sequences_per_genus, include_only=include,
                                exclude_taxonomy_level=taxonomy, seed=seed)
        try:
            generate_silva_datasets(silva_by_taxonomy_path, output_path + '_excluded', num_of_datasets=2,
                                    num_of_sequences=num_of_sequences_per_genus, exclude_only=include,
                                    exclude_taxonomy_level=taxonomy, seed=seed)
        except ZeroDivisionError:
            print(f'Dataset {include} has too few sequences to make exclude data')

print('done!')

