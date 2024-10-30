import re
import sys
import networkx as nx
import argparse

def generate_peptides(G, k):
    peptides = []

    def dfs(v, path):
        if len(path) > k:
            return
        peptides.append(path)
        for w in G.neighbors(v):
            dfs(w, path + [w])

    for node in G.nodes():
        dfs(node, [node])

    return peptides

def read_adjlist_with_attributes(filename):
    G = nx.DiGraph()
    with open(filename, 'r') as file:
        for line in file:
            node, neighbor, num_proteins = line.strip().split()
            G.add_edge(node, neighbor, num_proteins=int(num_proteins))
    return G

def write_peptides_to_fasta(peptides, filename):
    with open(filename, 'w') as file:
        for i, peptide in enumerate(peptides):
            file.write(f">peptide_{i+1}\n")  # Write the header
            file.write(''.join(peptide) + "\n")  # Write the peptide sequence

# Check that the correct number of arguments was given
if len(sys.argv) != 4:
    print("Usage: python Graph2Pep.py <input_file> <output_file> <k>")
    sys.exit(1)

# Get the arguments
GraphFile = sys.argv[1]
PepFile = sys.argv[2]
k = int(sys.argv[3])

G = read_adjlist_with_attributes(GraphFile);
# Generate the peptides
peptides = generate_peptides(G, k)
# Output the peptides
write_peptides_to_fasta(peptides, PepFile)
