import torch
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import torch.nn as nn
from torch.nn import functional as F
import sys
import shap
from sklearn import metrics
from scipy.stats import skew
#import gffpandas.gffpandas as gffpd

def loadData(input_path):
    # Sample ID and order that has both gene expression and DNA methylation data
    sample_id = np.loadtxt(input_path + 'both_samples.tsv', delimiter='\t', dtype='str')
    # Loading label
    label = pd.read_csv(input_path + 'both_samples_tumour_type_digit.tsv', sep='\t', header=0, index_col=0)
    class_num = len(label.tumour_type.unique())
    label_array = label['tumour_type'].to_numpy()
    label_df = label['tumour_type']
    return class_num, label_array, label_df, sample_id

def processPhenotypeDataForSamples(clinical_df, sample_id, le):
    phenotype = clinical_df
    phenotype = phenotype.loc[sample_id, :]
    return phenotype



def booleanConditional(sampleOne, sampleTwo, sampleThree, sampleFour, sampleFive):
    condition = np.logical_or.reduce(sampleOne, sampleTwo,sampleThree,sampleFour,sampleFive)
    return condition


def randomTrainingSample(expr,sampleSize):
    randomTrainingSampleexpr = expr.sample(n=sampleSize, axis=0)
    return randomTrainingSampleexpr




def splitExprandSample(condition, sampleSize, expr):
    expr_df_T = expr
    split_expr = expr_df_T[condition]
    split_expr = split_expr.sample(n=sampleSize, axis=0)
    return split_expr

def splitForGenders(sample_id):
    phenotype = pd.read_csv('DataSources/GDC-PANCAN.basic_phenotype.tsv', sep='\t', header=0, index_col=0)
    phenotype = phenotype.T
    phenotype = phenotype[sample_id].T
    female = phenotype['Gender'] == "Female"
    male = phenotype['Gender'] == "Male"
    return female, male

def printConditionalSelection(conditional,label_array):
    malecounts = label_array[conditional]
    unique, counts = np.unique(malecounts.iloc[:, 0], return_counts=True)




def addToTensor(expr_selection,device):
    selection = expr_selection.values.astype(dtype='float32')
    selection = torch.Tensor(selection).to(device)
    return selection

class ModelWrapper(nn.Module):
    def __init__(self, vae_model, source):
        super(ModelWrapper, self).__init__()
        self.vae_model = vae_model
        self.source = source

    def forward(self, input):
        return self.vae_model.source_predict(input, self.source)

def explain_vae(sample_id, vae_model, expr_df, clinical_df, source, subtype, device='cpu', le=None):
    """
    :param sample_id: List of samplesid to consider for explanation.
    :param vae_model: The CustOmics model to explain, the output should be a 1xn_class tensor.
    :param expr_df: DataFrame of data to explain, the input has to match the source selected for the explanation.
    :param clinical_df: DataFrame containing the clinical data.
    :param source: Omic source to explain.
    :param subtype: Subtype that needs explanation, shap values will be computed for this subtypes against the others.
    :param le: LabelEncoder to revert back from numerical encoding to original names.
    :return:None, plots of the top 10 genes and their shap values
    """
    #This class has combined the different analysis' of the Deep SHAP values we conducted.
    #SHAP reference: Lundberg et al., 2017: http://papers.nips.cc/paper/7062-a-unified-approach-to-interpreting-model-predictions.pdf

    sample_id = list(set(sample_id).intersection(set(expr_df.index)))
    phenotype = processPhenotypeDataForSamples(clinical_df, sample_id, le)
    print("BRCA subtype is")
    print(subtype)

    conditionaltumour=phenotype.loc[:, 'PAM50'] == subtype

        
    print('Printing the dataframe')
    print(expr_df)
    expr_df = expr_df.loc[sample_id,:]
    normal_expr = randomTrainingSample(expr_df, 10)
    tumour_expr = splitExprandSample(condition=conditionaltumour, sampleSize=10, expr=expr_df)
    # put on device as correct datatype
    background = addToTensor(expr_selection=normal_expr, device=device)
    print("Background Shape : {}".format(background.shape))
    male_expr_tensor = addToTensor(expr_selection=tumour_expr, device=device)
    print("Tumor Shape : {}".format(male_expr_tensor.shape))


    e = shap.DeepExplainer(ModelWrapper(vae_model, source=source), background)
    print("calculating shap values")
    shap_values_female = e.shap_values(male_expr_tensor, ranked_outputs=None)

    shap.summary_plot(shap_values_female[0],features=tumour_expr,feature_names=list(tumour_expr.columns), show=False, plot_type="violin", max_display=10, plot_size=[4,6])
    plt.savefig('shap_{}_{}.png'.format(source, subtype), bbox_inches='tight')
    print("calculated shap values")