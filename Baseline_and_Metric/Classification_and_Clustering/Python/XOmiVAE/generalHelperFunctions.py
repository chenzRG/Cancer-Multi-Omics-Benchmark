import torch
import numpy as np
from torch.nn import functional as F
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def loadData(input_path):
    # Sample ID and order that has both gene expression and DNA methylation data
    sample_id = np.loadtxt(input_path + 'both_samples.tsv', delimiter='\t', dtype='str')
    # Loading label
    label = pd.read_csv(input_path + 'both_samples_tumour_type_digit.tsv', sep='\t', header=0, index_col=0)
    class_num = len(label.tumour_type.unique())
    label_array = label['tumour_type'].to_numpy()
    return class_num, label_array, sample_id

def processPhenotypeDataForSamples(sample_id):
    phenotype = pd.read_csv('DataSources/GDC-PANCAN.basic_phenotype.tsv', sep='\t', header=0, index_col=0)
    phenotype = phenotype.T
    phenotype = phenotype[sample_id]
    phenotype = phenotype.T
    return phenotype

def processDataForDEA(sample_id,expr_df,label_array,tumourType):
    phenotype = processPhenotypeDataForSamples(sample_id)
    conditionone = phenotype['sample_type'] == "Solid Tissue Normal"
    conditiontwo = phenotype['project_id'] == tumourType
    conditioncombined = np.logical_and(conditionone, conditiontwo)
    expr_df_T = expr_df.T
    normal_BRCA = expr_df_T[conditioncombined]
    # np.savetxt('data/normal_BRCA.csv', normal_BRCA.index, delimiter=",", fmt="%s")
    import csv
    data = ["%s" % i for i in normal_BRCA.index]
    out = csv.writer(open("data/normal_BRCA.csv", "w"), delimiter=',', quoting=csv.QUOTE_ALL)
    out.writerow(data)
    sample_id_with_chosen_labels = sample_id[np.nonzero(label_array == 3)]
    brca_tumour = expr_df[sample_id_with_chosen_labels]
    brca_tumour = brca_tumour.T
    data = ["%s" % i for i in brca_tumour.index]
    out = csv.writer(open("../data/tumour_brca.csv", "w"), delimiter=',', quoting=csv.QUOTE_ALL)
    out.writerow(data)


def preprocessExpr_df(expr_df):
    expr_df.info()
    expr_df = expr_df.apply(pd.to_numeric, errors='coerce')
    expr_df.info()
    check = expr_df.isnull().sum().sum()
    expr_df = expr_df.fillna(0)
    checkafter = expr_df.isnull().values.any()
    return expr_df


def checkPredictions(expr_tensor,vae_model):
    vae_model.eval()
    with torch.no_grad():
        _, recon_data, mean, log_var, pred_y = vae_model(expr_tensor)
        print("predicted y is")
        pred_y_softmax = F.softmax(pred_y, dim=1)
        _, predicted = torch.max(pred_y_softmax, 1)
        print(predicted)

def addToTensor(expr_selection,device):
    selection = expr_selection.values.astype(dtype='float32')
    selection = torch.Tensor(selection).to(device)
    selection = torch.transpose(selection, 0, 1)
    return selection

def saveLatentSpace(vae_model,expr_tensor):
    vae_model.eval
    with torch.no_grad():
        z, recon_data, mean, log_var, gpred_y = vae_model(expr_tensor)
    z = z.detach().numpy()
    np.savetxt('save_z.csv', z, delimiter=',')


def plotLatentSpaceTCGATsne(latentSpace, labels):
    from sklearn.manifold import TSNE
    tsne = TSNE(n_components=2, random_state=0)
    # data X should be the z values
    print("1")
    tsne_obj = tsne.fit_transform(latentSpace)
    print("2")
    tsne_df = pd.DataFrame({'X': tsne_obj[:, 0],
                            'Y': tsne_obj[:, 1],
                            'digit': labels})
    print("3")
    print(tsne_df.head())
    colours = ['dimgray', 'indianred', 'darkred', 'salmon', 'orangered',
               'chocolate', 'burlywood', 'khaki', 'y', 'darkolivegreen', 'springgreen', 'darkslategrey', 'aqua',
               'powderblue', 'mediumblue', 'darkslateblue', 'blueviolet', 'purple', 'magenta', 'pink', 'k',
               'yellow', 'darkgray', 'lightgrey', 'aquamarine', 'pink', 'lemonchiffon', 'rosybrown', 'c',
               'greenyellow', 'sienna', 'deeppink', 'teal', 'mediumslateblue']
    sns.scatterplot(x="X", y="Y",
                    hue="digit",
                    palette=sns.color_palette(colours, n_colors=tsne_df.digit.unique().shape[0]),
                    legend='full',
                    data=tsne_df,
                    s=15);
    print("4")
    plt.show()

def processSubtypeSamples(sample_id,subtypeOne,subtypeTwo):
    subtypes = pd.read_csv('data/TCGASubtypeAnalysis.csv', sep=',', header=0, index_col=1)
    relevantsubtypes = subtypes.loc[sample_id]
    conditionalOne = relevantsubtypes['SUBTYPE'] == subtypeOne
    conditionalTwo = relevantsubtypes['SUBTYPE'] == subtypeTwo
    return conditionalOne, conditionalTwo

def sampleSameAmount(label_array, label_dig, no_labels):
    countlabel = 0
    target_array = np.zeros(label_array.shape)
    for x in range(0, len(label_array)):
        if label_array[x] == label_dig and countlabel < no_labels:
            countlabel += 1
            target_array[x] = 1
        else:
            target_array[x] = 0
    target_array = target_array.astype(bool)
    return target_array


def plotThreeBarCharts():
    KIRC_shap_values = pd.read_csv('data/NormalvsKIRCasBackground_feature_importance_absolute.csv', header=0,
                                   index_col=0)
    KIRC_shap_values = KIRC_shap_values.head(15)
    KIRP_shap_values = pd.read_csv('data/NormalvsKIRPasBackground_feature_importance_absolute.csv', header=0,
                                   index_col=0)
    KIRP_shap_values = KIRP_shap_values.head(15)
    KICH_shap_values = pd.read_csv('data/NormalvsKICHasBackground_feature_importance_absolute.csv', header=0,
                                   index_col=0)
    KICH_shap_values = KICH_shap_values.head(15)
    fig, axs = plt.subplots(ncols=3)
    sns.set_style("white")
    fig.tight_layout(pad=6.0)
    sns.barplot(x="feature_importance_vals", y="gene", data=KIRC_shap_values, color="skyblue", ax=axs[0])
    sns.barplot(x="feature_importance_vals", y="gene", data=KIRP_shap_values, color="skyblue", ax=axs[1])
    sns.barplot(x="feature_importance_vals", y="gene", data=KICH_shap_values, color="skyblue", ax=axs[2])
    axs[0].set(xlabel="Mean |SHAP value| ", ylabel="Gene")
    axs[1].set(xlabel="Mean |SHAP value| ", ylabel="Gene")
    axs[2].set(xlabel="Mean |SHAP value| ", ylabel="Gene")
    axs[0].set_title('KIRC vs. Training Sample', pad=15)
    axs[1].set_title('KIRP vs. Training Sample', pad=15)
    axs[2].set_title('KICH vs. Training Sample', pad=15)
    plt.savefig("ksNormalvstissue.png", dpi=1500)
    plt.show()






