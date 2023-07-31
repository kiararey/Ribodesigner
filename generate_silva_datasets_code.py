import os
import numpy as np
from numpy.random import default_rng
from silva_sequence_functions import generate_silva_datasets
seed = 1


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


