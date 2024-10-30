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
third_column = []

def output_protein_path(protein_paths, filename):
    with open(filename, 'w') as file:
        for name, path in protein_paths:
            file.write(f"> {name}\n {' -> '.join(path)}\n")

# First command line argument is the file name
file_name = sys.argv[1]
# Second command line argument is the output graph file name
graph_file_name = sys.argv[2]
# Third command line argument is the output path file name
path_file_name = sys.argv[3]


third_column = []

i = 0
length_threshold = 6 #peptides shorter than this will be discarded

G = nx.MultiDiGraph()
protein_paths = []
allpeptides = [];
allmispeptides = [];

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
            sequences = parse_fasta(protein_seq_file)
            for name, sequence in sequences:
                #print(f"Name: {name}\nSequence: {sequence}\n")
                #peptides = tryptic_digestion_single(sequence, length_threshold)
                #allpeptides.extend(peptides);
                #mispeptides = tryptic_digestion_mis(sequence, length_threshold, 2)
                #allmispeptides.extend(mispeptides);
                pairs = tryptic_digestion(sequence, length_threshold)
                #print(pairs)
                #print(allpeptides);
                #print(allmispeptides);
                #input()
                # Add edges between adjacent peptides, this will also automatically add the nodes
                protein_path = []
                for peptide1, peptide2, discarded_peptides in pairs:
                    protein_path.append(peptide1)
                    for d_peptides in discarded_peptides:
                        if len(d_peptides) > 0:
                            protein_path.append(d_peptides)
                    single_peptide = "".join(discarded_peptides);
                    #print(single_peptide);
                    label = 0;
                    edge_dict = G.get_edge_data(peptide1, peptide2);
                    if edge_dict is not None:
                        for edge in edge_dict.values():
                            if edge.get('discarded_peptides') == single_peptide:
                                label = 1;
                                break;
                    if label == 0:
                        G.add_edge(peptide1, peptide2, discarded_peptides=single_peptide)
                protein_path.append(peptide2)
                #protein_path.append()
                seqname=protein_seq_file+name;
                protein_paths.append((seqname, protein_path))
                #print(seqname);
            #allpeptides = list(set(allpeptides));
            #num_allpeptides = len(allpeptides);
            #allmispeptides = list(set(allmispeptides));
            #num_allmispeptides = len(allmispeptides);
            #print(i, num_allpeptides, num_allmispeptides);
            #print(f"Num paths: {len(protein_paths)}\n"); 
            #print(f"Num_codes: {len(G.nodes())}\n");
            #input()


# Draw the graph
#nx.draw(G, with_labels=True, node_size=500, node_color="skyblue")
#plt.savefig("protein_graph.pdf", format='pdf')

#nx.write_adjlist(G, "panproteome_graph.adjlist")
write_adjlist_with_attributes(G, graph_file_name)
output_protein_path(protein_paths, path_file_name)
#output_peptides_fasta(allpeptides, 'tryptic_peptides.fa')
#output_peptides_fasta(allmispeptides, 'mis_tryptic_peptides.fa')

# Output the graph is xml format
#nx.write_graphml(G, "panproteome_graph.graphml")
