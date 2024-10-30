import networkx as nx
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def read_adjlist_with_attributes(filename):
    """Reads the adjacency list and builds a directed graph with intermediate peptides as edge attributes."""
    G = nx.MultiDiGraph()  # Use MultiDiGraph to support multiple edges

    with open(filename, 'r') as file:
        for line in file:
            nodes = line.strip().split()
            main_node = nodes[0]
            neighbors = nodes[1:]

            # Add edges with the intermediate peptide stored as an attribute if it exists
            if len(neighbors) == 1:
                G.add_edge(main_node.strip(), neighbors[0].strip(), intermediate=None)
            elif len(neighbors) == 2:
                G.add_edge(main_node.strip(), neighbors[0].strip(), intermediate=neighbors[1].strip())

    return G

def get_custom_components(G):
    """Extracts components by starting from each main node."""
    components = [list(c) for c in nx.weakly_connected_components(G)]
    return components

# Load the graph with intermediate peptides as edge labels
file_path = "PanProteomeGraph_Uniprot.adjlist"
G = read_adjlist_with_attributes(file_path)

# Extract components
components = get_custom_components(G)

# Calculate component sizes
component_sizes = [len(c) for c in components]

# Set thresholds for small, middle, and large components
small_threshold = 2000
large_threshold = 2000

# Separate components into small, middle, and large categories
small_components = [size for size in component_sizes if size < small_threshold]
middle_components = [size for size in component_sizes if small_threshold <= size < large_threshold]
large_components = [size for size in component_sizes if size >= large_threshold]

# Calculate statistics
stats = {
    "Total Nodes": G.number_of_nodes(),
    "Total Edges": G.number_of_edges(),
    "Average Degree": sum(dict(G.degree()).values()) / G.number_of_nodes(),
    "Graph Density": nx.density(G),
    "Number of Custom Components": len(components),
    "Small Components (size < 100)": len(small_components),
    "Middle Components (100 <= size < 2000)": len(middle_components),
    "Large Components (size >= 2000)": len(large_components),
    "Singleton Nodes": sum(1 for size in small_components if size == 1),
}

# Display statistics
df_statistics = pd.DataFrame(list(stats.items()), columns=["Statistic", "Value"])
print(df_statistics)

# Visualization: Histograms of small, middle, and large components
fig, axs = plt.subplots(1, 2, figsize=(18, 5))
axs[0].hist(small_components, bins=20, edgecolor='black')
axs[0].set_title("Small Component Sizes")
axs[0].set_xlabel("Size")
axs[0].set_ylabel("Frequency")

# axs[1].hist(middle_components, bins=20, edgecolor='black')
# axs[1].set_title("Middle Component Sizes")
# axs[1].set_xlabel("Size")
# axs[1].set_ylabel("Frequency")

axs[1].hist(large_components, bins=20, edgecolor='black')
axs[1].set_title("Large Component Sizes")
axs[1].set_xlabel("Size")
axs[1].set_ylabel("Frequency")
axs[1].xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
plt.show()

# Visualization: Box plot of component sizes
plt.figure(figsize=(10, 5))
plt.boxplot([small_components, middle_components, large_components], vert=False, patch_artist=True, showfliers=True)
plt.yticks([1, 2, 3], ['Small', 'Middle', 'Large'])
plt.title("Box Plot of Component Sizes")
plt.xlabel("Component Size")
plt.show()

# Define the degree threshold
threshold = 20

# Get the degree of each node and count the frequency of each degree
degree_sequence = [d for n, d in G.degree()]
degree_count = Counter(degree_sequence)

# Separate degrees below and above the threshold
low_degrees = {k: v for k, v in degree_count.items() if k < threshold}
high_degrees = {k: v for k, v in degree_count.items() if k >= threshold}

# Plot degrees below the threshold (Bar plot)
plt.figure(figsize=(10, 5))
plt.bar(low_degrees.keys(), low_degrees.values(), color='blue', edgecolor='black')
plt.xlabel("Degree", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.show()

# Plot degrees above the threshold (Scatter plot)
plt.figure(figsize=(12, 6))
plt.scatter(high_degrees.keys(), high_degrees.values(), marker='o', color='red')
plt.xlabel("Degree", fontsize=20)
plt.ylabel("Frequency", fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()
# Extract nodes with a degree over 300
high_degree_nodes = [(n, d) for n, d in G.degree() if d > 300]
nodes_degree_2 = [node for node, degree in G.degree() if degree == 2]
print(f"Number of nodes with degree = 2: {len(nodes_degree_2)}")
# Display the nodes and their degrees
print(f"Nodes with degree > 300 (Total: {len(high_degree_nodes)}):")
for node in high_degree_nodes:
    in_edges = G.in_degree(node)
    out_edges = G.out_degree(node)
    print(f"Node: {node}, In-degree: {in_edges}, Out-degree: {out_edges}")

# Visualization: Draw the largest custom component (if it exists)
if components:
    largest_component = max(components, key=len)
    H = G.subgraph(largest_component).copy()

    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(H)  # Use a layout for better visualization
    nx.draw(H, pos, with_labels=True, node_size=50, font_size=8, edge_color='gray')

    # Draw edge labels (intermediate peptides) if available
    edge_labels = nx.get_edge_attributes(H, 'label')
    nx.draw_networkx_edge_labels(H, pos, edge_labels=edge_labels, font_size=6)

    plt.title("Visualization of the Largest Custom Component")
    plt.show()
