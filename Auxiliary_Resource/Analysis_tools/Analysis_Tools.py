import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy.plot import barplot, dotplot

def read_gene_expression_data(file_path):

#     Reads gene expression data from a CSV file.
    
#     Parameters:
#     file_path (str): Path to the CSV file containing the gene expression data.
    
#     Returns:
#     DataFrame: Gene expression matrix.

    gene_expression_matrix = pd.read_csv(file_path, index_col=0)
    return gene_expression_matrix

def perform_differential_expression_analysis(data, group1, group2):

#     Performs differential gene expression analysis between two groups.
    
#     Parameters:
#     data (DataFrame): Gene expression matrix.
#     group1 (list): List of sample names for group 1.
#     group2 (list): List of sample names for group 2.
    
#     Returns:
#     DataFrame: DataFrame with log2 fold change and adjusted p-values for each gene.

    group1_data = data[group1]
    group2_data = data[group2]

    # Perform t-tests
    t_stat, p_values = ttest_ind(group1_data, group2_data, axis=1)

    # Calculate log2 fold change
    log2_fc = np.log2(group1_data.mean(axis=1) / group2_data.mean(axis=1))
    
    # Adjust p-values using Benjamini-Hochberg correction
    from statsmodels.stats.multitest import multipletests
    _, adj_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    results = pd.DataFrame({
        'log2FC': log2_fc,
        'p-value': p_values,
        'adj_p-value': adj_p_values
    }, index=data.index)
    
    return results

def plot_volcano(results, p_value_threshold=0.05, log2fc_threshold=1):

#     Plots a volcano plot from the differential expression analysis results.
    
#     Parameters:
#     results (DataFrame): DataFrame with log2 fold change and adjusted p-values for each gene.
#     p_value_threshold (float): p-value threshold for significance.
#     log2fc_threshold (float): log2 fold change threshold for significance.

    plt.figure(figsize=(10, 8))
    
    # Scatter plot
    plt.scatter(results['log2FC'], -np.log10(results['adj_p-value']), color='grey')
    
    # Highlight significant points
    significant_up = (results['adj_p-value'] < p_value_threshold) & (results['log2FC'] > log2fc_threshold)
    significant_down = (results['adj_p-value'] < p_value_threshold) & (results['log2FC'] < -log2fc_threshold)
    
    plt.scatter(results[significant_up]['log2FC'], -np.log10(results[significant_up]['adj_p-value']), color='red', label='Upregulated')
    plt.scatter(results[significant_down]['log2FC'], -np.log10(results[significant_down]['adj_p-value']), color='blue', label='Downregulated')
    
    # Add lines for thresholds
    plt.axhline(y=-np.log10(p_value_threshold), color='black', linestyle='--')
    plt.axvline(x=log2fc_threshold, color='black', linestyle='--')
    plt.axvline(x=-log2fc_threshold, color='black', linestyle='--')
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted p-value')
    plt.title('Volcano Plot')
    plt.legend()
    plt.show()

def main(input_file, group1, group2):
    # Read the gene expression data
    data = read_gene_expression_data(input_file)
    
    # Perform differential gene expression analysis
    results = perform_differential_expression_analysis(data, group1, group2)
    
    # Plot the volcano plot
    plot_volcano(results)

if __name__ == "__main__":
    # Example usage
    input_file = 'path_data.csv'
    group1 = ['sample1', 'sample2', 'sample3'] 
    group2 = ['sample4', 'sample5', 'sample6']  
    
    main(input_file, group1, group2)
    
def perform_kegg_pathway_analysis(results, organism='Human', gene_column='log2FC', p_value_cutoff=0.05):

    # Filter significant genes
    significant_genes = results[results['adj_p-value'] < p_value_cutoff]
    
    # Create a dictionary for gene scores
    gene_scores = significant_genes[gene_column].to_dict()
    
    # Convert keys to a list
    gene_list = list(gene_scores.keys())
    
    # Perform KEGG pathway analysis
    enrichment_results = gp.enrichr(gene_list=gene_list,
                                    gene_sets='KEGG_2021_Human',
                                    organism=organism,
                                    outdir=None)
    
    # Convert results to DataFrame
    enrichment_results_df = enrichment_results.results
    
    return enrichment_results_df

def main(input_file, group1, group2):
    # Read the gene expression data
    data = read_gene_expression_data(input_file)
    
    # Perform differential gene expression analysis
    results = perform_differential_expression_analysis(data, group1, group2)
    
    # Plot the volcano plot
    plot_volcano(results)
    
    # Perform KEGG pathway analysis
    kegg_results = perform_kegg_pathway_analysis(results)
    
    # Plot KEGG pathway enrichment results
    try:
        barplot(kegg_results, title='KEGG Pathway Enrichment')
        plt.show()
    except ValueError as e:
        print(f"Warning: {e}")

if __name__ == "__main__":
    # Example usage
    input_file = 'path_data.csv'
    group1 = ['sample1', 'sample2', 'sample3']  
    group2 = ['sample4', 'sample5', 'sample6']  
    
    main(input_file, group1, group2)