import pickle
import networkx as nx
import re
from collections import defaultdict

UNIPROT_COLORS = {
    1: "red",
    2: "green",
    3: "blue"
}

def load_data(filename):
    """Load data from a pickle file."""
    with open(filename, 'rb') as f:
        return pickle.load(f)

def filter_components_by_size_and_protein(components, protein_mapping, sizes=[96, 192, 288]):
    """Filter components by size and separate them based on UniProt proteins."""
    with_uniprot = {}
    without_uniprot = {}

    for idx, component_nodes in enumerate(components):
        component_name = f"Component_{idx + 1}"

        if component_name in protein_mapping:
            num_proteins = len(protein_mapping[component_name])
            if num_proteins in sizes:
                if any(is_uniprot_protein(p) for p in protein_mapping[component_name]):
                    print(f"Component {component_name} has {num_proteins} proteins with UniProt.")
                    with_uniprot[component_name] = component_nodes
                else:
                    without_uniprot[component_name] = component_nodes

    print(f"Components with UniProt: {len(with_uniprot)}, without UniProt: {len(without_uniprot)}")
    return with_uniprot, without_uniprot

def read_adjlist_with_attributes(filename):
    """Reads the adjacency list and builds a directed graph with attributes."""
    G = nx.MultiDiGraph()

    with open(filename, 'r') as file:
        for line in file:
            nodes = line.strip().split()
            main_node = nodes[0].strip().lower()
            neighbors = [n.strip().lower() for n in nodes[1:]]

            if len(neighbors) == 1:
                G.add_edge(main_node, neighbors[0], intermediate="None")
            elif len(neighbors) == 2:
                G.add_edge(main_node, neighbors[0], intermediate=neighbors[1])

    return G

def is_uniprot_protein(protein_name):
    """Check if the protein is from UniProt."""
    return bool(re.match(r"^[A-Z0-9]+_[A-Z]+$", protein_name))

def load_selected_protein_data(filename, relevant_proteins):
    """Load only relevant UniProt protein data from the protein path file."""
    peptide_to_proteins = defaultdict(set)

    with open(filename, 'r') as file:
        while True:
            protein_line = file.readline().strip()
            path_line = file.readline().strip()
            if not protein_line or not path_line:
                break

            # Extract UniProt protein name
            uniprot_match = re.search(r"\|([A-Z0-9]+_[A-Z]+)\b", protein_line)
            protein_name = uniprot_match.group(1) if uniprot_match else None

            # Map peptides only if the protein is relevant
            if protein_name and protein_name in relevant_proteins:
                print(f"Found relevant protein: {protein_name}")
                peptides = [p.strip().lower() for p in path_line.split(' -> ')]
                for peptide in peptides:
                    peptide_to_proteins[peptide].add(protein_name)

    return peptide_to_proteins
def assign_protein_colors_for_component(proteins):
    """Assign colors only to UniProt proteins within a component."""
    # Define up to 7 distinct colors
    colors = ["red", "blue", "green", "purple", "orange", "yellow", "cyan"]
    protein_to_color = {}

    # Filter out non-UniProt proteins (e.g., ENSP) before assigning colors
    uniprot_proteins = [p for p in proteins if re.match(r"^[A-Z0-9]+_[A-Z]+$", p)]

    # Assign colors to the first 7 UniProt proteins
    for i, protein in enumerate(sorted(uniprot_proteins)):
        if i < 7:
            protein_to_color[protein] = colors[i]

    print(f"Assigned colors to UniProt proteins: {protein_to_color}")
    return protein_to_color

def set_node_colors_for_pie(G, peptide_to_proteins, protein_to_color):
    """Assign pie chart data to nodes based on UniProt proteins."""
    # Initialize pie chart colors with 0 for all 7 colors
    color_template = {
        "red": 0.0, "blue": 0.0, "green": 0.0, 
        "purple": 0.0, "orange": 0.0, "yellow": 0.0, "cyan": 0.0, 
        "black": 0.0  # Default color if no UniProt proteins are present
    }

    for node in G.nodes:
        proteins = peptide_to_proteins.get(node, set())

        # Start with the template for each node
        color_distribution = color_template.copy()

        # Assign colors only if the proteins exist in the color mapping
        valid_proteins = [p for p in proteins if p in protein_to_color]
        if valid_proteins:
            for protein in valid_proteins:
                color = protein_to_color[protein]
                color_distribution[color] = 1.0  # Assign the color
        else:
            # If no UniProt protein, assign black = 1
            color_distribution["black"] = 1.0

        # Store the pie chart data in node attributes
        for color, value in color_distribution.items():
            G.nodes[node][f"pie_{color}"] = value

def recreate_and_export_subgraphs(components, original_graph, peptide_to_proteins, protein_mapping, output_prefix="component"):
    """Recreate subgraphs, assign colors, and export them to GraphML."""
    for component_name, nodes in components.items():
        subgraph = original_graph.subgraph(nodes).copy()

        # Get the proteins associated with this component
        component_proteins = protein_mapping.get(component_name, [])

        # Assign colors only to UniProt proteins within the component
        protein_to_color = assign_protein_colors_for_component(component_proteins)

        # Set pie chart colors for the subgraph nodes
        set_node_colors_for_pie(subgraph, peptide_to_proteins, protein_to_color)

        # Replace None values in edge attributes with "None"
        for u, v, data in subgraph.edges(data=True):
            for key, value in data.items():
                if value is None:
                    data[key] = "None"

        # Export the subgraph as a directed GraphML file
        output_file = f"{output_prefix}_{component_name}.graphml"
        nx.write_graphml(subgraph, output_file, encoding='utf-8')
        print(f"Exported {component_name} to {output_file}")


def find_relevant_proteins(protein_mapping, selected_components):
    """Extract relevant protein names for the selected components."""
    relevant_proteins = set()
    for component_name in selected_components:
        proteins = set(protein_mapping[component_name])
        relevant_proteins.update(proteins)
    return relevant_proteins

# Load components and protein mappings
components = load_data("components.pkl")
protein_mapping = load_data("component_to_protein_mapping.pkl")

# Load the original graph
graph_file = "PanProteomeGraph_Uniprot.adjlist"
G = read_adjlist_with_attributes(graph_file)

# Select components by size (e.g., 384, 672 proteins)
with_uniprot, without_uniprot = filter_components_by_size_and_protein(
    components, protein_mapping, sizes=[ 96*4,96*6]
)

# Combine a few components from both categories
selected_components = {}
selected_components.update(dict(list(with_uniprot.items())[:20]))  # 20 with UniProt

print(selected_components)
# Find relevant proteins for the selected components
relevant_proteins = find_relevant_proteins(protein_mapping, selected_components)
# Load only relevant peptide-to-protein mappings
protein_path_file = "protein_uniprot_paths.txt"
peptide_to_proteins = load_selected_protein_data(protein_path_file, relevant_proteins)
print(peptide_to_proteins)

# Recreate and export the selected components to GraphML
recreate_and_export_subgraphs(selected_components, G, peptide_to_proteins, protein_mapping)
