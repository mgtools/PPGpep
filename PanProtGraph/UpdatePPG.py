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

def write_adjlist_with_attributes(G, filename):
    with open(filename, 'w') as file:
        for node in G.nodes():
            neighbors = G.neighbors(node)
            for neighbor in neighbors:
                edge_dict = G[node][neighbor];
                for edge in edge_dict:
                    single_peptide = edge_dict[edge].get('interspersed_peptides')
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

def load_multigraph(filename):
  """Loads a multigraph from a file.

  Args:
    filename: The name of the file to load the graph from.

  Returns:
    A NetworkX graph object.
  """

  with open(filename, "r") as f:
    graph = nx.MultiDiGraph()
    for line in f:
      line = line.strip()
      if not line:
        continue

      tryptic_peptides = line.split()
      if len(tryptic_peptides) == 2:
          single_peptide = "";
      elif len(tryptic_peptides) == 3:
          single_peptide = tryptic_peptides[2];
      graph.add_edge(tryptic_peptides[0], tryptic_peptides[1], interspersed_peptides=single_peptide);

  return graph

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

# First command line argument is the current graph file name
graph_file_name = sys.argv[1]
# Second command line argument is the index file for individual protein sequences
index_file_name = sys.argv[2]
# third command line argument is the output graph file name
outgraph_file_name = sys.argv[3]
# forth command line argument is the output path file name
outpath_file_name = sys.argv[4]

third_column = []

i = 0
length_threshold = 6 #peptides shorter than this will be discarded

G = nx.MultiDiGraph()
protein_paths = []
allpeptides = [];
allmispeptides = [];

#load the current multigraph
G = load_multigraph(graph_file_name)
print(G.number_of_edges())

with open(index_file_name, 'r') as file:
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
                            if edge.get('interspersed_peptides') == single_peptide:
                                label = 1;
                                break;
                    if label == 0:
                        #print(peptide1, peptide2, single_peptide)
                        #input();
                        G.add_edge(peptide1, peptide2, interspersed_peptides=single_peptide)
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
print(G.number_of_edges())
write_adjlist_with_attributes(G, outgraph_file_name)
output_protein_path(protein_paths, outpath_file_name)

# Output the graph is xml format
#nx.write_graphml(G, "panproteome_graph.graphml")
