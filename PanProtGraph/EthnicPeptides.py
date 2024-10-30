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

def tryptic_digestion_mis(protein_sequence, length_threshold, k):
    # Regular expression for Trypsin cleavage rule
    pattern = r"(?<=[KR])(?!P)"

    # Split the protein sequence into peptides
    peptides = re.split(pattern, protein_sequence)

    for i in range(len(peptides)):
        for j in range(i, min(i+k+1, len(peptides))):
            peptides.append(''.join(peptides[i:j+1]))

    #filter short peptides (< length_threshold)
    peptides = [peptide for peptide in peptides if len(peptide) >= length_threshold]

    return peptides;

def tryptic_digestion_single(protein_sequence, length_threshold):
    # Regular expression for Trypsin cleavage rule
    pattern = r"(?<=[KR])(?!P)"

    # Split the protein sequence into peptides
    peptides = re.split(pattern, protein_sequence)

    #filter short peptides (< length_threshold)
    peptides = [peptide for peptide in peptides if len(peptide) >= length_threshold]

    return peptides;

def tryptic_digestion(protein_sequence, length_threshold):
    # Regular expression for Trypsin cleavage rule
    pattern = r"(?<=[KR])(?!P)"

    # Split the protein sequence into peptides
    peptides = re.split(pattern, protein_sequence)

    pairs = []
    for i in range(len(peptides) - 1):
        if len(peptides[i]) >= length_threshold and len(peptides[i + 1]) >= length_threshold:
            pairs.append((peptides[i], peptides[i + 1], []))
        elif len(peptides[i]) >= length_threshold:
            j = i + 1
            discarded_peptides = []
            while j < len(peptides) and len(peptides[j]) < length_threshold:
                discarded_peptides.append(peptides[j])
                j += 1
            if j < len(peptides):
                pairs.append((peptides[i], peptides[j], discarded_peptides))

    return pairs

def write_adjlist_with_attributes(G, filename):
    with open(filename, 'w') as file:
        for node in G.nodes():
            neighbors = G.neighbors(node)
            for neighbor in neighbors:
                edge_dict = G[node][neighbor];
                for edge in edge_dict:
                    single_peptide = edge_dict[edge].get('discarded_peptides')
                    if single_peptide is None:
                        single_peptide = " ";
                    file.write(f"{node} {neighbor} {single_peptide}\n")

def output_peptides_fasta(protein_list, filename):
    with open(filename, 'w') as f:
        for i, protein in enumerate(protein_list):
            # Write the description line
            f.write(f">Protein_{i+1}\n")
            # Write the protein sequence
            f.write(protein + "\n")

# Use the function
#import sys

#file_name = sys.argv[1]
#sequences = parse_fasta(file_name)
#for name, sequence in sequences:
#    print(f"Name: {name}\nSequence: {sequence}\n")
#

def output_protein_path(protein_paths, filename):
    with open(filename, 'w') as file:
        for name, path in protein_paths:
            file.write(f"> {name}\n {' -> '.join(path)}\n")

# First command line argument is the file name
file_name = sys.argv[1]


i = 0
length_threshold = 6 #peptides shorter than this will be discarded

G = nx.MultiDiGraph()
peptides_by_class = {}
ethnic_peptides = []

with open(file_name, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
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
            ethnicgroup = row[1];
            if ethnicgroup not in peptides_by_class:
                peptides_by_class[ethnicgroup] = []
            sequences = parse_fasta(protein_seq_file)
            for name, sequence in sequences:
                peptides = tryptic_digestion_single(sequence, length_threshold)
                for peptide in peptides:
                    ethnic_peptides.append((ethnicgroup, peptide))
                    #if peptide not in peptides_by_class[ethnicgroup]:
                    #    peptides_by_class[ethnicgroup].append(peptide)

for ethnicgroup, peptide_count in peptides_by_class.items():
    print(ethnicgroup, peptide_count)

with open("peptides_by_ethnics.txt", "w") as f:
    for pair in ethnic_peptides:
            f.write(f"{pair[0]};{pair[1]}\n")

#print(peptides_by_class)
