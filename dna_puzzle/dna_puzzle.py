"""DNA Puzzle - DNA Elements Assembly Script

This script processes DNA elements according to a specified order and creates assembled DNA files.

The script takes two main inputs:
1. An Excel file containing the ordered list of DNA elements to be assembled
2. A CSV source file containing the DNA information for each element (one row per final file)

The script then:
- Reads and validates input files
- Extracts DNA sequences and information for each element
- Assembles the DNA elements in the specified order
- Generates output files in either FASTA+GFF or GenBank format

Usage:
    python dna_puzzle.py <excel_file> <source_csv> --output-format <format>

Args:
    excel_file: Path to Excel file with DNA element order
    source_csv: Path to CSV file with DNA element information
    source_dir: Path to directory containing DNA element files (SnapGene files)
    output_format: Desired output format ('fasta_gff' or 'genbank')

Returns:
    Creates assembled DNA sequence files in the specified format
"""
import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, SimpleLocation
from Bio.Seq import Seq

import helpers

# Argparse
parser = argparse.ArgumentParser(description='DNA Puzzle - DNA Elements Assembly Script')
parser.add_argument(
    '--order_csv',
    help='''CSV with DNA element order. First column "filename", rest are element names.
    Each row is a file to be assembled. The order of the elements is given by the columns.
    Each row can have as many columns as you like. Please leave columns un-named. '''
)
parser.add_argument(
    '--source_csv',
    help='''CSV with DNA info. First column matches element IDs from order_csv.
    Second column can be:
    a) Raw DNA sequence (feature named as first column)
    b) SnapGene file (*.dna) to add as-is
    c) <feature_id>:<snapgene_file.dna> to add specific feature
    d) <snapgene_file.dna>:<start>:<end> to add specific region
    
    In all cases, the feature name in the output column will be the value of the first column.'''
)
parser.add_argument(
    '--source_dir',
    help='Directory containing any SnapGene files referenced in the source_csv',
)
parser.add_argument(
    '--output_format',
    help='Output format (fasta_gff or genbank)',
    required=True
)
parser.add_argument(
    '--output_dir',
    help='Directory to save output files',
    default='output'
)

args = parser.parse_args()

print(args)

order_dict = helpers.read_order_csv(args.order_csv)
sources_dict = helpers.parse_dna_info_file_to_feature_dict(args.source_csv, args.source_dir)

print("heyyy")
# constructed_dna_fragment = None
for key, value in order_dict.items():
    print(f'Processing {key}...')
    this_seqrecord = helpers.construct_order_of_elements_to_dna(order=value, source_dict=sources_dict)
    helpers.write_seqrecord_to_file(this_seqrecord, args.output_format, args.output_dir, key)
        

# source_dict = parse_dna_info(args.source_csv)
# print(order_dict)
# print(source_dict)
# raw_seq_test = create_feature_from_raw_sequence('100k02_LSup', source_dict['100k02_LSup'][1])
# raw_seq_from_snapgene_test = create_feature_from_snapgene_file('100k02_LSup', os.path.join(args.source_dir,"1A2_2025_pSC101_EcTrpRS-h14_colR_VF.dna"))
# chonks_test = source_dict['100k21_chonks']
# raw_seq_from_snapgene_region_test = create_feature_from_snapgene_region(os.path.join(args.source_dir, chonks_test[1][0]), chonks_test[1][1], chonks_test[1][2])
# print(raw_seq_from_snapgene_region_test)
# print(raw_seq_from_snapgene_region_test[0:20].seq)
# print(raw_seq_from_snapgene_region_test[-20::].seq)

# fstest = source_dict["1A1_sC"]
# seq_search_test = create_feature_from_searching_snapgene_with_feature_label(os.path.join(args.source_dir,fstest[1][1]), fstest[1][0], fstest[1][0])
# print(seq_search_test)



print("Done")