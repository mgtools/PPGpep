import networkx as nx
import pickle

def extract_all_component_peptides(graph, components):
    """Extract peptides from all components and store them in a dictionary."""
    component_peptides_dict = {}
    for idx, component in enumerate(components):
        component_name = f"Component_{idx + 1}"
        peptides = set()

        for u, v, data in graph.subgraph(component).edges(data=True):
            peptides.add(u.strip().lower())
            peptides.add(v.strip().lower())

            intermediate = data.get('intermediate')
            if intermediate and len(intermediate) > 5:
                peptides.add(intermediate.strip().lower())

        component_peptides_dict[component_name] = peptides
        print(f"Extracted {len(peptides)} peptides for {component_name}")
    return component_peptides_dict

import re

def extract_protein_name(protein_line):
    """Extract a meaningful identifier from a protein line."""
    try:
        ensembl_match = re.search(r"(ENSP|ENSG|ENST)\d+\.\d+", protein_line)
        if ensembl_match:
            return ensembl_match.group(0)

        uniprot_match = re.search(r"\|([A-Z0-9]+_[A-Z]+)\b", protein_line)
        if uniprot_match:
            return uniprot_match.group(1)

        if protein_line.startswith('>'):
            return protein_line.split(' ')[1].strip()  # Second word after '>'

    except Exception as e:
        print(f"Error parsing protein line: {protein_line} - {e}")

    # Return None if no valid identifier is found
    return None

def match_proteins_to_components(protein_file, component_peptides_dict):
    """Match proteins to components and return the mapping."""
    component_to_protein = {name: [] for name in component_peptides_dict.keys()}
    unmatched_proteins = set()

    with open(protein_file, 'r') as infile:
        while True:
            protein_line = infile.readline().strip()
            path_line = infile.readline().strip()

            if not protein_line or not path_line:  # End of file or empty line
                break

            # Extract protein name using enhanced logic
            protein_name = extract_protein_name(protein_line)
            if not protein_name:  # Skip if protein name extraction failed
                print(f"Skipped invalid protein line: {protein_line}")
                continue

            # Normalize peptides from the path line
            peptides = {p.strip().lower() for p in path_line.split(' -> ') if len(p.strip()) > 5}

            matched = False
            for component_name, component_peptides in component_peptides_dict.items():
                if peptides.issubset(component_peptides):
                    component_to_protein[component_name].append(protein_name)
                    print(f"Matched > {protein_name} to {component_name}")
                    matched = True
                    break

            if not matched:
                unmatched_proteins.add(protein_name)

    if unmatched_proteins:
        print(f"Warning: {len(unmatched_proteins)} proteins were not matched.")

    return component_to_protein


def save_data(data, filename):
    """Save data to a pickle file."""
    with open(filename, 'wb') as f:
        pickle.dump(data, f)
    print(f"Data saved to {filename}")

def load_data(filename):
    """Load data from a pickle file."""
    with open(filename, 'rb') as f:
        return pickle.load(f)

def read_adjlist_with_attributes(filename):
    """Reads the adjacency list and builds a directed graph with intermediate peptides as attributes."""
    G = nx.MultiDiGraph()

    with open(filename, 'r') as file:
        for line in file:
            nodes = line.strip().split()
            main_node = nodes[0].strip().lower()
            neighbors = [n.strip().lower() for n in nodes[1:]]

            if len(neighbors) == 1:
                G.add_edge(main_node, neighbors[0], intermediate=None)
            elif len(neighbors) == 2:
                G.add_edge(main_node, neighbors[0], intermediate=neighbors[1])

    return G

graph_file = "PanProteomeGraph_Uniprot.adjlist"
protein_path_file = "protein_uniprot_paths.txt"

G = read_adjlist_with_attributes(graph_file)
components = list(nx.weakly_connected_components(G))

component_peptides_dict = extract_all_component_peptides(G, components)

component_to_protein_mapping = match_proteins_to_components(protein_path_file, component_peptides_dict)

save_data(components, "components.pkl")
save_data(component_peptides_dict, "component_peptides_dict.pkl")
save_data(component_to_protein_mapping, "component_to_protein_mapping.pkl")
