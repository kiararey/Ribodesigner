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
silva_by_taxonomy_path = f'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/{taxonomy_level}'

# Here is the file where we want to save the dataset
output_path = 'Datasets_used/SILVA_squished_datasets'

# There are only Bacteria in this dataset, so will only do that:
fasta_file_names = np.array([file_name for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta' in file_name])

# Prepare list to store data in:
num_of_datasets = 5
num_of_sequences = 5

# Generate the squished datasets
generate_silva_datasets(silva_by_taxonomy_path, output_path, num_of_datasets=num_of_datasets,
                        num_of_sequences=num_of_sequences, seed=seed)


########################################################
# To generate datasets for targeted designs by order level test

# Here is the file where the dataset is. I've already run a script on this to separate it by levels of taxonomy.
# To do this, SequencePrepper.py on the SILVA database file.
silva_by_taxonomy_path = 'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Genus'

# Here is the file where we want to save the dataset
output_path_entero = 'Datasets_used/SILVA_squished_datasets/Enterobacterales_only_squished'
output_path_pseudo = 'Datasets_used/SILVA_squished_datasets/Pseudomonadales_only_squished'
output_path_both = 'Datasets_used/SILVA_squished_datasets/Pseudo_and_entero_only_squished'
output_path_background_no_e = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished_no_entero'
output_path_background_no_p = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished_no_pseudo'
output_path_background_no_p_or_e = 'Datasets_used/SILVA_squished_datasets/Background_Bacteria_squished'

# There are only Bacteria in this dataset, so will only do that:
fasta_file_names = np.array([file_name for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta' in file_name])

# Prepare list to store data in:
num_of_datasets = 5
num_of_sequences_per_genus = 5

include_only_entero = ['Enterobacterales']
include_only_pseudo = ['Pseudomonadales']
include_only_both = ['Pseudomonadales', 'Enterobacterales']
taxonomy_level_of_inclusion = 'Order'

# Generate the squished datasets
# generate_silva_datasets(silva_by_taxonomy_path, output_path_entero, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_entero,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# generate_silva_datasets(silva_by_taxonomy_path, output_path_pseudo, num_of_datasets=num_of_datasets,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_pseudo,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)

generate_silva_datasets(silva_by_taxonomy_path, output_path_both, num_of_datasets=num_of_datasets,
                        num_of_sequences=num_of_sequences_per_genus, include_only=include_only_both,
                        exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)

# Generate background datasets (exclude one or the other or both)
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_e, num_of_datasets=1,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_entero,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)
#
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_p, num_of_datasets=1,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_pseudo,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)

generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_p_or_e, num_of_datasets=1,
                        num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_both,
                        exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed)


