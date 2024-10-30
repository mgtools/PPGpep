import pickle
import matplotlib.pyplot as plt
from collections import Counter

def load_data(filename):
    """Efficiently load data from a pickle file."""
    with open(filename, 'rb') as f:
        return pickle.load(f)

def split_clusters(mapping_file, small_threshold=200, medium_threshold=5000):
    """Split components into small, medium, and large clusters."""
    component_to_protein_mapping = load_data(mapping_file)

    # Initialize lists for the three categories
    small_clusters = []
    medium_clusters = []
    large_clusters = []

    # Categorize components based on the number of mapped proteins
    for proteins in component_to_protein_mapping.values():
        count = len(proteins)
        if count <= small_threshold:
            small_clusters.append(count)
        elif count <= medium_threshold:
            medium_clusters.append(count)
        else:
            large_clusters.append(count)

    print(f"Small clusters: {len(small_clusters)}")
    print(f"Medium clusters: {len(medium_clusters)}")
    print(f"Large clusters: {len(large_clusters)}")

    return small_clusters, medium_clusters, large_clusters

def plot_cluster_distribution(cluster_data, title, output_file):
    """Plot the distribution of clusters."""
    frequency = Counter(cluster_data)
    counts, frequencies = zip(*sorted(frequency.items()))

    # Plot the distribution
    plt.figure(figsize=(10, 5))
    plt.bar(counts, frequencies, width=0.8, align='center', edgecolor='black')
    plt.xlabel("Number of Proteins", fontsize=16)
    plt.ylabel("Frequency", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Save the plot as an image
    plt.savefig(output_file)
    print(f"Plot saved as '{output_file}'")

    # Display the plot
    plt.show()

# Load the mapping data and split clusters based on the thresholds
mapping_file = "component_to_protein_mapping.pkl"
small_clusters, medium_clusters, large_clusters = split_clusters(
    mapping_file, small_threshold=1000, medium_threshold=2000
)

# Plot and save the distributions for each category
plot_cluster_distribution(
    small_clusters,
    "Distribution of Small Clusters (<= 50 Proteins)",
    "small_clusters_distribution.png"
)
plot_cluster_distribution(
    medium_clusters,
    "Distribution of Medium Clusters (70 - 130 Proteins)",
    "medium_clusters_distribution.png"
)
plot_cluster_distribution(
    large_clusters,
    "Distribution of Large Clusters (> 20000 Proteins)",
    "large_clusters_distribution.png"
)
