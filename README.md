CDR3 Sequence Analysis Pipeline
This repository contains a Python script for analyzing CDR3 amino acid sequences, including:

Distance Matrix Calculation: Computes pairwise Levenshtein distances between sequences to assess their similarities.
Clustering: Uses hierarchical clustering to explore relationships among sequences.
Pearson Correlation Analysis: Evaluates co-variation between sequences.
Heatmap Visualization: Provides a clustered heatmap to visualize sequence similarities and correlations.
Network Analysis: Generates a network to visualize the relationships between sequences, including an exploration of associations with different ethnicities.
Features
Comprehensive Analysis: From sequence similarity calculation to network analysis, this script provides a full workflow for understanding relationships between immune sequences.
Parallel Processing: Efficient calculation of Levenshtein distances using parallel processing to handle large datasets.
Custom Visualizations: Heatmaps and network graphs to highlight key patterns in immune receptor diversity.
Prerequisites
Python 3.x
Required packages: scipy, pandas, seaborn, matplotlib, networkx, joblib, Levenshtein, numpy

Example Input File Structure
The input CSV file (example: CDR3_AA_w_Ethnicity.csv) should contain the following columns:

seqID: A unique identifier for each sequence (e.g., Seq1, Seq2, ...).
CDR3aa: The CDR3 amino acid sequence for each entry (e.g., CASSPGGTDTQYF).
Ethnicity: The ethnicity group to which the sequence belongs (e.g., BBR, ASI, LAT, BLA, CAU)
*modify to fit your needs

Usage
To run the analysis:

Prepare Input File: Ensure your CSV file is structured as shown above.
Run the Script: Execute the Python script to perform the analysis.
Output: The script will generate:
Distance Matrix CSV: Distmatrix_scores.csv
Top Correlations CSV: pearson_correlations_pvalues.csv, top_200_r_values.csv
Visualizations: Heatmaps and network plots.
Output Files
Distmatrix_scores.csv: A matrix of pairwise distances between sequences.
pearson_correlations_pvalues.csv: Pearson correlations and p-values for each sequence pair.
top_200_r_values.csv: Top 200 sequence correlations for further analysis.
Heatmap and Network Visualizations: Visual representations of sequence relationships and ethnicity associations.
Key Visualizations
Clustered Heatmap: Visualizes sequence similarities and clusters.
Network Graph: Shows relationships between sequences, with nodes colored by ethnicity and edges weighted by correlation significance.
Notes
The average linkage method is used for clustering to balance between extremes (e.g., single vs. complete linkage).
Ethnicity information is used to explore genetic diversity among different populations, highlighting population-specific immune responses.

Contributions are welcome! Feel free to submit issues or pull requests to enhance the features or fix bugs.
