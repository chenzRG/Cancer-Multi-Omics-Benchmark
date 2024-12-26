import torch
import numpy as np
from scipy import stats
import pandas as pd

def dataPrepForDeepShap(sample_ids, label_arrays, exprdf, number_of_samples, label_id,device):
    sample_id_with_chosen_labels = sample_ids[np.nonzero(label_arrays == label_id)]
    chosen_expr = exprdf[sample_id_with_chosen_labels]
    chosen_expr_feature_labels = chosen_expr.iloc[:, 0:5].index
    chosen_expr_feature_labels_list = chosen_expr_feature_labels.values.tolist()
    # find in gencode and replace with gene acronym
    chosen_expr_tensor = chosen_expr.iloc[:, 0:number_of_samples].values
    chosen_expr_tensor = chosen_expr_tensor.astype(dtype='float32')
    chosen_expr_tensor = torch.Tensor(chosen_expr_tensor)
    chosen_expr_tensor = chosen_expr_tensor.to(device)
    chosen_expr_tensor = torch.transpose(chosen_expr_tensor, 0, 1)
    return chosen_expr_feature_labels_list, chosen_expr_tensor, sample_id_with_chosen_labels


def sampleSameAmount(label_array, label_dig, no_labels, gender, phenotype):
    countlabel = 0
    target_array = np.zeros(label_array.shape)
    for x in range(0, len(label_array)):
        if label_array[x] == label_dig and countlabel < no_labels and phenotype['Gender'].iloc[x] == gender:
            countlabel += 1
            target_array[x] = 1
        else:
            target_array[x] = 0
    target_array = target_array.astype(bool)
    print("end count")
    print(countlabel)
    return target_array

def booleanConditional(sampleOne, sampleTwo, sampleThree, sampleFour, sampleFive):
    condition = np.logical_or.reduce(sampleOne, sampleTwo,sampleThree,sampleFour,sampleFive)
    return condition


def saveMostStatisticallySignificantIndex(conditionOne, conditionTwo, fileNameOne, fileNameTwo, z):
    female_z = z[conditionOne]
    male_z = z[conditionTwo]
    statistics = stats.ttest_ind(female_z, male_z, axis=0, equal_var=False, nan_policy='propagate')
    stat_lowest_index = np.argmin(statistics.pvalue)
    female_z_column = female_z[:, stat_lowest_index]
    male_z_column = male_z[:, stat_lowest_index]
    file_path_one = 'data/' + fileNameOne + '.csv'
    file_path_two = 'data/' + fileNameTwo + '.csv'
    np.savetxt(file_path_one, female_z_column)
    np.savetxt(file_path_two, male_z_column)

def randomTrainingSample(expr,sampleSize):
    randomTrainingSampleexpr = expr.sample(n=sampleSize, axis=1)
    return randomTrainingSampleexpr


def pathwayComparison(sample_id, pathway):
    pathways = pd.read_csv('data/PathwayAnalysis.csv', sep=',', header=0, index_col=0)
    relevantpathways = pathways.loc[sample_id]
    logicalMutant = relevantpathways[pathway] == 1
    logicalNonMutant = relevantpathways[pathway] == 0
    return logicalMutant,logicalNonMutant

def multipleSampling(label_array,phenotype, tumour_id_one, tumour_id_two, tumour_id_three, tumour_id_four, tumour_id_five,sample_number):
    female_thirteen = sampleSameAmount(label_array, tumour_id_one, sample_number, "Female", phenotype)
    female_eighteen = sampleSameAmount(label_array, tumour_id_two, sample_number, "Female", phenotype)
    female_four = sampleSameAmount(label_array, tumour_id_three, sample_number, "Female", phenotype)
    female_five = sampleSameAmount(label_array, tumour_id_four, sample_number, "Female", phenotype)
    female_six = sampleSameAmount(label_array, tumour_id_five, sample_number, "Female", phenotype)

    male_thirteen = sampleSameAmount(label_array, tumour_id_one, sample_number, "Male", phenotype)
    male_eighteen = sampleSameAmount(label_array, tumour_id_two, sample_number, "Male", phenotype)
    male_four = sampleSameAmount(label_array, tumour_id_three, sample_number, "Male", phenotype)
    #changed as only 25 to sample from
    male_five = sampleSameAmount(label_array, tumour_id_four, 25, "Male", phenotype)
    male_six = sampleSameAmount(label_array, tumour_id_five, sample_number, "Male", phenotype)

    maleCondition = booleanConditional(male_thirteen, male_eighteen, male_four, male_five, male_six)
    femaleCondition = booleanConditional(female_thirteen, female_eighteen, female_four, female_five,female_six)
    return femaleCondition, maleCondition

def splitExprandSample(condition, sampleSize, expr):
    expr_df_T = expr.T
    split_expr = expr_df_T[condition].T
    split_expr = split_expr.sample(n=sampleSize, axis=1)
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
    print("male sample id counts")
    print(np.asarray((unique, counts)).T)

#when not using ranked output i.e. not explaining the outputs (therefore exlaining the z dimension or mu)
def saveShapValues(shap_vals, gene, chrom, ranked_output=True):
    # if not using ranked outputs (i.e. wanting to explain just the top predicted label... then use vals = np.abs(shap_vals).mean(0)
    # if passing in shap_value[0] as ranked_outputs then again need to use shap_vals[0] here
    #if want to understand the positive and negative SHAP values then change absolute SHAP value here to separating between positive and negative values
    if ranked_output==True:
        vals = np.abs(shap_vals[0]).mean(0)
    else:
        vals = np.abs(shap_vals).mean(0)
    feature_importance = pd.DataFrame(list(zip(gene, chrom, vals)),
                                      columns=['gene', 'chrom', 'feature_importance_vals'])
    feature_importance.sort_values(by=['feature_importance_vals'], ascending=False, inplace=True)
    feature_importance.to_csv('data/shapValues.csv')
    return feature_importance

def getGenes(expr_df):
    # get genes and chromosomes
    gencode_ids = pd.read_csv('DataSources/gencodev22annotationgeneCOPY.tsv', sep='\t', header=0, index_col=0)
    new = expr_df.merge(gencode_ids, left_index=True, right_index=True, how='left')
    print("new")
    print(new)
    genes = new.iloc[:, -5]
    chrom = new.iloc[:, -4]
    ensg = new.index

    return genes, chrom, ensg

def getTopShapValues(shap_vals, numberOfTopFeatures, expr_df, ranked_output=True, cancerType="TCGA-BRCA",absolute=True):
    gene, chrom, ensg = getGenes(expr_df)
    if absolute:
        if ranked_output==True:
            print("here absolute and ranked")
            shap_value=shap_vals[0]
            vals = np.abs(shap_value[0]).mean(0)
        else:
            vals = np.abs(shap_vals).mean(0)
    else:
        print("should not print here")
        if ranked_output==True:
            shap_value=shap_vals[0]
            vals = shap_value[0].mean(0)
        else:
            print("should not print here")
            vals = shap_vals.mean(0)

    #feature_importance = pd.DataFrame(list(zip(gene, chrom, ensg, vals)),columns=['gene', 'chrom', 'id', 'feature_importance_vals'])
    feature_importance = pd.DataFrame(list(zip(gene, chrom, ensg, vals)),
                                      columns=['gene', 'chrom', 'id', 'feature_importance_vals'])
    feature_importance.sort_values(by=['feature_importance_vals'], ascending=False, inplace=True)
    #plotGenes(cancerType,feature_importance)

    mostImp_shap_values = feature_importance.head(numberOfTopFeatures)
    print(mostImp_shap_values)

    feature_importance.to_csv(cancerType +"_feature_imp_secondtry")
    """
    print(mostImp_shap_values)
    print("least importance absolute values")
    feature_importance.sort_values(by=['feature_importance_vals'], ascending=True, inplace=True)
    leastImp_shap_values = feature_importance.head(numberOfTopFeatures)
    print(leastImp_shap_values)
    """
    return mostImp_shap_values

def plotGenes(cancertype,shap_values):
    import matplotlib.pyplot as plt
    import seaborn as sns

    shap_values = shap_values.head(10)
    fig, axs = plt.subplots(ncols=1)
    sns.set_style("white")
    fig.tight_layout(pad=6.0)
    sns.barplot(x="feature_importance_vals", y="gene", data=shap_values, color="skyblue", ax=axs)
    axs.set(xlabel="Mean |SHAP value| ", ylabel="Gene")
    axs.set_title(cancertype, pad=15)
    plt.savefig(cancertype+".png", dpi=1500)
    #plt.show()




