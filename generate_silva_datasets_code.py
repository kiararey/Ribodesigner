from silva_sequence_functions import generate_silva_datasets
seed = 1

# Finally wrote this to be able to just handle the SILVA ref 99 dataset without pre-processing
silva_ref_99_fasta = 'Datasets_used/SILVA_138.1_SSURef_NR99_tax_silva.fasta'
silva_by_taxonomy_path = silva_ref_99_fasta
########################################################
# # To generate SILVA squished datasets
# # Determine at what taxonomic level we want the datasets. This must match one of the filenames from SequencePrepper.py:
# taxonomy_level = 'Order'
#
# # Here is the file where we want to save the dataset
# output_path = 'Datasets_used/SILVA_squished_datasets'
#
# # Prepare list to store data in:
# num_of_sequences_per_order = 10
#
# # Generate the squished datasets
# print('Making general datasets')
# generate_silva_datasets(silva_ref_99_fasta, output_path,
#                         num_of_sequences=num_of_sequences_per_order, seed=seed)
# print('Making bacterial datasets')
# generate_silva_datasets(silva_ref_99_fasta, output_path, include_only=['Bacteria'], exclude_taxonomy_level='Domain',
#                         num_of_sequences=num_of_sequences_per_order, seed=seed)
# print('Making eukaryotic datasets')
# generate_silva_datasets(silva_ref_99_fasta, output_path, include_only=['Eukaryota'], exclude_taxonomy_level='Domain',
#                         num_of_sequences=num_of_sequences_per_order, seed=seed)
# print('Making archaeal datasets')
# generate_silva_datasets(silva_ref_99_fasta, output_path, include_only=['Archaea'], exclude_taxonomy_level='Domain',
#                         num_of_sequences=num_of_sequences_per_order, seed=seed)


########################################################
# # To generate datasets for targeted designs by order level test
#
# # Here is the file where the dataset is. I've already run a script on this to separate it by levels of taxonomy.
# # To do this, SequencePrepper.py on the SILVA database file.
# # silva_by_taxonomy_path = 'Datasets_used/SILVA_Ref_NR_99_dataset_by_taxonomy_Bacteria_Only/Genus'
# silva_by_taxonomy_path = silva_ref_99_fasta
#
# # Here is the file where we want to save the dataset
# output_path_entero = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Enterobacterales_only_squished'
# output_path_pseudo = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Pseudomonadales_only_squished'
# output_path_both = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Pseudo_and_entero_only_squished'
# output_path_background_no_e = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished_no_entero'
# output_path_background_no_p = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished_no_pseudo'
# output_path_background_no_p_or_e = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Background_Bacteria_squished'
# output_path_background_gram_pos = 'Datasets_used/SILVA_squished_datasets_1_per_genus/Gram_positives_only'
# output_path_background_no_gram_pos = 'Datasets_used/SILVA_squished_datasets_1_per_genus/No_Gram_positives'
#
# # # There are only Bacteria in this dataset, so will only do that:
# # fasta_file_names = np.array([file_name for file_name in os.listdir(silva_by_taxonomy_path) if '.fasta' in file_name])
#
# num_of_sequences_per_genus = 1
#
# include_only_entero = ['Enterobacterales']
# include_only_pseudo = ['Pseudomonadales']
# include_only_both = ['Enterobacterales', 'Pseudomonadales']
# taxonomy_level_of_inclusion = 'Order'
# include_only_gram_pos = ['Actinobacteriota', 'Firmicutes']
#
# # Generate the squished datasets
# print('entero only')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_entero,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_entero,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed, divide_by='Genus')
#
# print('pseudo only')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_pseudo,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_pseudo,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed, divide_by='Genus')
#
# print('entero or pseudo')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_both,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_both,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed, divide_by='Genus')
#
# # Generate background datasets (exclude one or the other or both)
# print('all but entero')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_e,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_entero,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed, divide_by='Genus')
#
# print('all but pseudo')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_p,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_pseudo,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed, divide_by='Genus')
#
# print('all but pseudo or entero')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_p_or_e,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_both,
#                         exclude_taxonomy_level=taxonomy_level_of_inclusion, seed=seed, divide_by='Genus')
#
# print('gram positives only')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_gram_pos,
#                         num_of_sequences=num_of_sequences_per_genus, include_only=include_only_gram_pos,
#                         exclude_taxonomy_level='Phylum', seed=seed, divide_by='Genus')
#
#
# print('all but gram positives')
# generate_silva_datasets(silva_by_taxonomy_path, output_path_background_no_gram_pos,
#                         num_of_sequences=num_of_sequences_per_genus, exclude_only=include_only_gram_pos,
#                         exclude_taxonomy_level='Phylum', seed=seed, divide_by='Genus')

########################################################
# taxonomy_levels_all = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
# to_generate = {'Phylum': ['Proteobacteria', 'Firmicutes'],
#                'Class': ['Gammaproteobacteria', 'Bacilli'],
#                'Order': ['Enterobacterales', 'Pseudomonadales', 'Bacillales'],
#                'Family': ['Enterobacteriaceae', 'Pseudomonadaceae', 'Bacillaceae'],
#                'Genus': ['Escherichia-Shigella', 'Pseudomonas', 'Bacillus']}
# for taxonomy in taxonomy_levels_all:
#     print(f'\nNow generating data for {taxonomy}:')
#
#     for include in to_generate[taxonomy]:
#         print(f'\nDataset is now {include}:')
#         output_path = f'Datasets_used/SILVA_squished_datasets_5_per_genus/{taxonomy}_{include}'
#         generate_silva_datasets(silva_by_taxonomy_path, output_path + '_included',
#                                 num_of_sequences=5, include_only=include, unique_at='Species',
#                                 exclude_taxonomy_level=taxonomy, seed=seed, divide_by='Genus')
#         try:
#             generate_silva_datasets(silva_by_taxonomy_path, output_path + '_excluded',
#                                     num_of_sequences=5, exclude_only=include, unique_at='Species',
#                                     exclude_taxonomy_level=taxonomy, seed=seed, divide_by='Genus')
#         except ZeroDivisionError:
#             print(f'Dataset {include} has too few sequences to make exclude data')
# to_generate = ['Enterobacterales', 'Pseudomonadales']
# for include in to_generate:
#         output_path = f'Datasets_used/SILVA_squished_datasets_3000_per_order/Order_{include}'
#         generate_silva_datasets(silva_by_taxonomy_path, output_path + '_included', num_of_sequences=3000,
#                                 unique_at='Species', include_only=[include], exclude_taxonomy_level='Order',
#                                 seed=seed, pick_from_file=True)

print('Now time for fungi hell yeahhhhh ')
output_path = f'Datasets_used/SILVA_squished_datasets_fungi'
generate_silva_datasets(silva_by_taxonomy_path, output_path, num_of_sequences=39, include_only=['Saccharomyces cerevisiae'],
                        exclude_taxonomy_level='any', divide_by='Species', seed=seed, pick_from_file=True,
                        repeat_species=True, only_test=True)

generate_silva_datasets(silva_by_taxonomy_path, output_path, num_of_sequences=3000, divide_by='Family', unique_at='Species',
                        include_only=['Ascomycota', 'Basidiomycota'], exclude_taxonomy_level='Family', seed=seed,
                        pick_from_file=True)

print('done!')

