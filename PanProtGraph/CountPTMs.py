import csv
import re
import sys
import os
import networkx as nx
import matplotlib.pyplot as plt

def parse_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = []
        name = None

        # Read the sequences
        for line in file:
            # If the line starts with '>', we've reached a new sequence
            if line.startswith(">"):
                if name:  # if there was a previous sequence
                    sequences.append((name, sequence))
                name = line.strip()[1:]  # strip leading '>' and newline character
                sequence = ""
            else:
                # Remove newline characters and concatenate
                sequence += line.strip()

        # Add the last sequence
        if name:
            sequences.append((name, sequence))

    return sequences

import re

def count_ptm_sites(sequences):
  """Counts the number of putative sites for Phosphorylation, acetylation, N-glycosylation and O-glycosylation in protein sequences.

  Args:
    sequences: A list of protein sequences in fasta format.

  Returns:
    A dictionary mapping PTM type to the number of putative sites.
  """

  ptm_sites = {}
  for sequence in sequences:
    seq = sequence.strip()
    if not seq:
      continue

    # Count phosphorylation sites.
    phosphorylation_sites = re.findall(r'[STY]([STK])', seq)
    ptm_sites['phosphorylation'] = len(phosphorylation_sites)

    # Count acetylation sites.
    acetylation_sites = re.findall(r'[K]', seq)
    ptm_sites['acetylation'] = len(acetylation_sites)

    # Count N-glycosylation sites.
    n_glycosylation_sites = re.findall(r'N[^P][^T][^S][^K][^R][^Q][^H]', seq)
    ptm_sites['N-glycosylation'] = len(n_glycosylation_sites)

    # Count O-glycosylation sites.
    o_glycosylation_sites = re.findall(r'[S]', seq)
    ptm_sites['O-glycosylation'] = len(o_glycosylation_sites)

  return ptm_sites

# Use the function
#import sys

#file_name = sys.argv[1]
#sequences = parse_fasta(file_name)
#for name, sequence in sequences:
#    print(f"Name: {name}\nSequence: {sequence}\n")
#
third_column = []

# First command line argument is the file name
file_name = sys.argv[1]

third_column = []

i = 0

G = nx.MultiDiGraph()
ptm_site_counts = {}

with open(file_name, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    print("#Haplotype #phosphorylation #acetylation #N-glycosylation #O-glycosylation")
    for row in reader:
        if row and not row[0].startswith("Assembly") and len(row) >= 3:  # Ensure there is a third column
            i += 1
            protein_seq_file = "Homo_sapiens-" + row[2] + "-2022_07-pep.fa" #generate protein sequence file name
            if not os.path.isfile(protein_seq_file):
                protein_seq_file = "Homo_sapiens-" + row[2] + "-2022_08-pep.fa" #generate protein sequence file name
                if not os.path.isfile(protein_seq_file):
                    print(f"{protein_seq_file} does not exist")
                    exit(0)
            print(protein_seq_file)
            sequences = parse_fasta(protein_seq_file)
            ptm_site_count = count_ptm_sites(sequences)
            ptm_site_counts['phosphorylation'] += ptm_site_count['phosphorylation'];
            ptm_site_counts['acetylation'] += ptm_site_count['acetylation'];
            ptm_site_counts['N-glycosylation'] += ptm_site_count['N-glycosylation'];
            ptm_site_counts['O-glycosylation'] += ptm_site_count['O-glycosylation'];
            print(i, ptm_site_counts['phosphorylation'], ptm_site_counts['acetylation'], ptm_site_counts['N-glycosylation'], ptm_site_counts['O-glycosylation'])

