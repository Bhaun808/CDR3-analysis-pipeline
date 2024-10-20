# Combined code for distance matrix calculation, clustering, Pearson correlation analysis, heatmap visualization, and network analysis of CDR3 sequences

# This script performs multiple stages of analysis for CDR3 amino acid sequences, starting with the computation of a distance matrix using Levenshtein distances.
# The sequences are then clustered to explore relationships between them, and a Pearson correlation analysis is performed to understand the co-variation among sequences.
# The script also visualizes these correlations in a heatmap and uses network analysis techniques to explore connections between sequences, including their association with different ethnicities.
# This combined workflow aims to identify key patterns in immune receptor diversity, cluster similarities, and analyze how genetic relationships may correlate with ethnicity.

# Example Input File Structure
# The input CSV file ('Filtered_Sequences_90_Percent_Identical.csv') should be structured as follows:
# - It should contain at least three columns:
#   1. 'seqID': A unique identifier for each sequence (e.g., Seq1, Seq2, ...).
#   2. 'CDR3aa': The CDR3 amino acid sequence for each entry (e.g., CASSPGGTDTQYF).
#   3. 'Ethnicity': The ethnicity group to which the sequence belongs (e.g., BBR, ASI, LAT, BLA, CAU).
#
# Example CSV file (CDR3_AA_w_Ethnicity.csv):
#
# seqID,CDR3aa,Ethnicity
# Seq1,CASSPGGTDTQYF,BBR
# Seq2,CASSLGGADTQYF,ASI
# Seq3,CASRRGGDTQYF,LAT
# Seq4,CASSQGGTETQYF,BLA
# Seq5,CASSPGGTDTQYF,CAU
#
# Make sure the CSV file follows this structure to ensure the script functions correctly.

from scipy.cluster.hierarchy import linkage, leaves_list
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from scipy.spatial.distance import pdist, squareform
import Levenshtein as lev
import networkx as nx
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from joblib import Parallel, delayed
import multiprocessing

# Step 1: Create Distance Matrix

# The Levenshtein distance is used here to measure the similarity between CDR3 amino acid sequences.
# It calculates the minimum number of single-character edits (insertions, deletions, substitutions) required to change one sequence into another.
# This metric is particularly useful for comparing immune receptor sequences, as small changes in these sequences can significantly impact their function.

# Read sequences and ethnicity from CSV file
data = pd.read_csv('CDR3_AA_w_Ethnicity.csv')

# Extract sequences and ethnicity
sequences = data['CDR3aa'].tolist()
ethnicity = data['Ethnicity'].tolist()
sequence_ids = [f"Seq{i}" for i in range(len(sequences))]  # Generating sequence IDs

# Create an empty matrix for storing the distances
dist_matrix = np.zeros((len(sequences), len(sequences)))

# Calculate the pairwise Levenshtein distance with parallel processing
# The nested loop for calculating pairwise distances can be computationally intensive for large datasets.
# To optimize this process, we use parallelization with joblib, which can significantly speed up the computation.

def calculate_distance(i, j):
    return lev.distance(sequences[i], sequences[j])

num_cores = multiprocessing.cpu_count()

distance_results = Parallel(n_jobs=num_cores)(
    delayed(calculate_distance)(i, j)
    for i in range(len(sequences))
    for j in range(i + 1, len(sequences))
)

# Fill in the distance matrix with the calculated values
k = 0
for i in range(len(sequences)):
    for j in range(i + 1, len(sequences)):
        dist_matrix[i, j] = distance_results[k]
        dist_matrix[j, i] = distance_results[k]  # The distance is symmetric
        k += 1

# Cluster the sequences
try:
    condensed_dist_matrix = squareform(dist_matrix)
    # The 'average' linkage method was chosen because it tends to provide a balanced clustering, reducing the influence of outliers that might otherwise dominate the results.
    # Unlike 'single' linkage, which can create long chains of loosely connected sequences, or 'complete' linkage, which can be overly influenced by the most distant points, 'average' linkage offers a compromise that often leads to more meaningful clusters in biological data.
    linkage_matrix = linkage(condensed_dist_matrix, method='average')
    dendro_order = leaves_list(linkage_matrix)

    # Reorder the matrix and sequence ids
    clustered_matrix = dist_matrix[dendro_order, :][:, dendro_order]
    clustered_ids = [sequence_ids[i] for i in dendro_order]

    # Convert the matrix to a DataFrame for easier handling
    df = pd.DataFrame(clustered_matrix, index=clustered_ids, columns=clustered_ids)
    df.to_csv("Distmatrix_scores.csv")

    # Define a custom colormap for the heatmap
    colors = ["#D218D9", "#FFFFFF", "#16B8D9"]  # Pink, white, blue
    cmap = LinearSegmentedColormap.from_list("mycmap", colors)
    ethnicity_colors = {'BBR': '#AB9493', 'ASI': '#9493AB', 'LAT': '#9FCABF', 'BLA': '#8DB3CE', 'CAU': '#9CADA9'}
    default_color = '#808080'
    row_colors = pd.Series([ethnicity_colors.get(e, default_color) for e in ethnicity], index=clustered_ids)

    # Plot a clustered heatmap
    g = sns.clustermap(df, cmap=cmap, center=0, robust=True, row_cluster=True, col_cluster=False,
                       cbar_pos=(0.02, 0.8, .02, .2), linewidths=0,
                       cbar_kws={"orientation": "vertical", "shrink": .5, "drawedges": False},
                       figsize=(5, 5),
                       row_colors=row_colors,
                       xticklabels=False, yticklabels=False)

    g.cax.set_position([0.01, .85, .02, .1])
    g.cax.set_yticks([0, 25, 50, 75, 100])
    g.cax.set_yticklabels(['0', '25', '50', '75', '100'])
    g.ax_heatmap.tick_params(left=False, bottom=False, right=False)
    g.fig.suptitle("CDR3 Distance", y=0.92)
    plt.show()

except ValueError as e:
    print(f"An error occurred during clustering: {e}")

# Step 2: Generate Pearson Correlation

# Pearson correlation is used here to assess the linear relationship between different sequences based on their distance matrix values.
# This analysis helps identify pairs of sequences that change in similar ways, indicating possible functional similarities or shared features.
# The resulting correlation values can reveal patterns of similarity that are not immediately obvious from the distance matrix alone.

# Prepare a list to store the correlations and p-values
results = []

# Compute the Pearson correlation for each pair of sequences (columns)
for i, col1 in enumerate(df.columns):
    for col2 in df.columns[i + 1:]:  # Avoid repeating pairs and self-comparison
        corr, p_val = pearsonr(df[col1], df[col2])
        results.append({
            'Sequence_1': col1,
            'Sequence_2': col2,
            'Pearson_Correlation': corr,
            'P_value': p_val
        })

# Convert the results list to a DataFrame
results_df = pd.DataFrame(results)
results_df.to_csv('pearson_correlations_pvalues.csv', index=False)

# Step 3: Extract Top Values and Plot Heatmap

# The top 200 correlations are extracted to focus on the most significant relationships between sequences.
# These top correlations can help identify clusters of sequences that exhibit strong similarities, which may indicate common functionality or evolutionary relationships.
# By visualizing these in a heatmap, we can quickly identify patterns and subgroups within the dataset, providing insight for further analysis.

results_df_top = results_df.nlargest(200, 'Pearson_Correlation')
results_df_top.to_csv('top_200_r_values.csv', index=False)

pivot_df = results_df_top.pivot(index='Sequence_1', columns='Sequence_2', values='Pearson_Correlation')
plt.figure(figsize=(10, 8))
sns.heatmap(pivot_df, annot=True, fmt=".2f", cmap='coolwarm', cbar_kws={'label': 'Pearson Correlation'})
plt.title('Heatmap of Top 200 Pearson Correlations between CDR3 Sequences')
plt.tight_layout()
plt.show()

# Step 4: Network Analysis with Ethnicities

# Including ethnicity in the network analysis allows us to explore whether specific genetic relationships or sequence similarities are associated with particular ethnic groups.
# This can help reveal how genetic diversity is distributed across different populations, providing insights into evolutionary pressures or the impact of certain diseases in different ethnicities.
# By visualizing the network with ethnicities, we can identify patterns that may indicate population-specific immune responses or other unique characteristics that are relevant to immune function.

seq_id_to_ethnicity = dict(zip(data['seqID'], data['Ethnicity']))

G = nx.Graph()

min_pval = results_df['P_value'].min()
max_neg_log_pval = -np.log10(min_pval if min_pval > 0 else 1e-100)  # avoid log(0)
edge_width_scale = 5 / max_neg_log_pval

for _, row in results_df_top.iterrows():
    weight = -np.log10(row['P_value'] if row['P_value'] > 0 else 1e-100) * edge_width_scale
    G.add_edge(row['Sequence_1'], row['Sequence_2'], weight=weight, r_value=row['Pearson_Correlation'])

node_r_values = {node: np.mean([G[node][neigh]['r_value'] for neigh in G.neighbors(node)]) for node in G.nodes()}
node_colors = [(node_r_values[node] - 0.9) * 10 for node in G.nodes()]
ethnicity_colors = {'BBR': '#AB9493', 'ASI': '#9493AB', 'LAT': '#9FCABF', 'BLA': '#8DB3CE', 'CAU': '#9CADA9'}
node_ethnicity_colors = [ethnicity_colors.get(seq_id_to_ethnicity.get(node, None), '#808080') for node in G.nodes()]

filtered_nodes = [node for node in G.nodes() if node in seq_id_to_ethnicity]
filtered_node_colors = [color for node, color in zip(G.nodes(), node_colors) if node in seq_id_to_ethnicity]
filtered_node_ethnicity_colors = [color for node, color in zip(G.nodes(), node_ethnicity_colors) if node in seq_id_to_ethnicity]
node_sizes = [G.degree(node) * 600 for node in filtered_nodes]

pos = nx.spring_layout(G, k=1.5, iterations=100)
nx.draw_networkx_edges(G, pos, edge_color='#D3D3D3', width=[G[u][v]['weight'] for u, v in G.edges()])

# Visualization explanation: 
# - Node colors represent the average Pearson correlation (R-value) with neighboring nodes. Darker colors indicate higher similarity.
# - Edge colors represent ethnicity information, allowing us to see how nodes (sequences) group by population.
# - Node size is determined by the degree (number of connections), indicating how central a node is within the network.

nodes = nx.draw_networkx_nodes(G, pos, nodelist=filtered_nodes, node_color=filtered_node_colors, cmap=node_cmap, edgecolors=filtered_node_ethnicity_colors, linewidths=6, node_size=node_sizes)
nx.draw_networkx_labels(G, pos, font_size=8)

sm = plt.cm.ScalarMappable(cmap=node_cmap, norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = plt.colorbar(sm)
cbar.set_label('Average R-value')
