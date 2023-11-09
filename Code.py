from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage


# Read in gene sequences and their corresponding names from file
with open('1249 DNA seq.fasta', 'r') as f:
    gene_data = [line.strip() for line in f if line.strip()]

# Split gene data into sequence names and sequences
sequence_names = [line[1:] for i, line in enumerate(gene_data) if i % 2 == 0]
gene_sequences = [line for i, line in enumerate(gene_data) if i % 2 == 1]

# Calculate pairwise similarities between all sequences
similarities = np.zeros((len(gene_sequences), len(gene_sequences)))
for i, j in combinations(range(len(gene_sequences)), 2):
    seq1, seq2 = gene_sequences[i], gene_sequences[j]
    if len(seq1) > 0:
        similarity = sum(1 for a, b in zip(seq1, seq2) if a == b) / len(seq1)
    else:
        similarity = 0.0
    similarities[i, j] = similarity
    similarities[j, i] = similarity

# Perform hierarchical clustering on the similarities
Z = linkage(similarities, method='complete')

# Plot dendrogram with sequence names
plt.figure(figsize=(10, 5))
plt.title('RGAs sequence pairwise similarities-hierarchical clustering', fontsize =15)
plt.xlabel('RGAs', fontsize= 15)
plt.ylabel('Distance', fontsize=15)
dendrogram(Z, leaf_rotation=90., leaf_font_size=8., labels=sequence_names)
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import pandas as pd

# Read in gene sequences and their corresponding names from a file
with open('1249 DNA seq.fasta', 'r') as f:
    gene_data = [line.strip() for line in f if line.strip()]

# Split gene data into sequence names and sequences
sequence_names = [line[1:] for i, line in enumerate(gene_data) if i % 2 == 0]
gene_sequences = [line for i, line in enumerate(gene_data) if i % 2 == 1]

# Calculate pairwise similarities between all sequences
similarities = np.zeros((len(gene_sequences), len(gene_sequences)))
for i in range(len(gene_sequences)):
    for j in range(i + 1, len(gene_sequences)):
        seq1, seq2 = gene_sequences[i], gene_sequences[j]
        if len(seq1) > 0:
            similarity = sum(1 for a, b in zip(seq1, seq2) if a == b) / len(seq1)
        else:
            similarity = 0.0
        similarities[i, j] = similarity
        similarities[j, i] = similarity

# Perform hierarchical clustering on the similarities
Z = linkage(similarities, method='complete')

# Set the number of clusters
num_clusters = 5   # Set the cluster number based on the number of clusters  from the above code

# Get cluster labels
cluster_labels = fcluster(Z, num_clusters, criterion='maxclust')

# Create a dictionary to store clusters
gene_clusters = {cluster_id: [] for cluster_id in range(1, num_clusters + 1)}

# Assign genes to clusters based on cluster labels
for i, cluster_id in enumerate(cluster_labels):
    gene_clusters[cluster_id].append((sequence_names[i], gene_sequences[i]))

# Export genes of each cluster to separate CSV files
for cluster_id, genes in gene_clusters.items():
    cluster_csv = pd.DataFrame(genes, columns=['Sequence Name', 'Gene Sequence'])
    cluster_csv.to_csv(f'cluster_{cluster_id}.csv', index=False)



