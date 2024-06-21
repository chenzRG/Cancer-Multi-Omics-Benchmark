
import pandas as pd
from dimensionKnockOutExperiments import nullingDimensions
from omiShapExplainer import omiShapExplainer
import numpy as np
import torch


if __name__ == "__main__":
    expr_path = 'data/GDC-PANCAN_htseq_fpkm_'
    input_path = 'DataSources/GDC-PANCAN_'

    print('Loading data...')
    expr_df = pd.read_pickle("expr_df.pkl")
    sample_id = np.loadtxt(input_path + 'both_samples.tsv', delimiter='\t', dtype='str')
    label = pd.read_csv(input_path + 'both_samples_tumour_type_digit.tsv', sep='\t', header=0, index_col=0)
    label_array = label['tumour_type'].to_numpy()


    #explain interim layer. Need to pass in tumour name.
    #omiShapExplainer(sample_id, label_array, expr_df, tumourName="TCGA-BRCA", NormalvsTumourInterimExplain=True)

    #Example of knocking out dimension. Pass in chosen tumour tissue and ID to evaluate. Fill in dimensions within the TestingNullingDimensions code.
    #e.g. BRCA tumour ID=3, lUAD=17
    omiShapExplainer(sample_id, label_array, expr_df, tumourID=3, tumourName="TCGA-BRCA", TestingNullingDimensions=True)

    # example of explaining the most important genes for a tissue
    #omiShapExplainer(sample_id, label_array, expr_df, NormalvsTumourExplain=True, tumourName="TCGA-LUAD")
  

    # explain the most important dimension in the supervised part of the model. Pass in dimension number and tumour name.
    #omiShapExplainer(sample_id, label_array, expr_df, tumourName="TCGA-HNSC", dimension=42,NormalvsTumourDimensionExplain=True)

