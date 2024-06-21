import torch
import generalHelperFunctions as helper
import numpy as np

from exprVAEwithAdditionalFeatures import ExprOmiVAE
import generalHelperFunctions as GeneralHelper
import shapExplainerHelper as ShapHelper
import disentanglingBarChartPlots as barChartHelper
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import LogNorm


def disentanglingFunctions(sample_id,label_array,fullTsnePlot=False, subtypeTsne=False,
               dimensionstandardDevs=False, heatmaps=False,brcaBarCharts=True,traverseLatentSpace=False):
    vae_model = ExprOmiVAE.VAE()

    if fullTsnePlot:
        """
        #method to load and save latent space if required
        expr_tensor = GeneralHelper.addToTensor(expr_df, device)
        vae_model.load_state_dict(torch.load('data/beta15.pt', map_location=torch.device('cpu')))
        ShapHelper.saveLatentSpace(vae_model, expr_tensor)
        """
        z = np.genfromtxt('data/z_before_supervised_normalvae_32.csv')
        y = label_array
        GeneralHelper.plotLatentSpaceTCGATsne(latentSpace=z, labels=y)

    if subtypeTsne:
        Basal, LumB = GeneralHelper.processSubtypeSamples(sample_id, subtypeOne="Basal",subtypeTwo="LumB")
        LumA, Her2 = GeneralHelper.processSubtypeSamples(sample_id,subtypeOne="LumA",subtypeTwo="Her2")
        #Get normal BRCA tissue for comaprison in t-SNE
        phenotype = GeneralHelper.processPhenotypeDataForSamples(sample_id)
        conditionone = phenotype['sample_type'] == "Solid Tissue Normal"
        conditiontwo = phenotype['project_id'] == "TCGA-BRCA"
        normal = np.logical_and(conditionone, conditiontwo)

        combined = np.logical_or.reduce(Basal, LumB, LumA, Her2,normal)

        z = np.genfromtxt('data/mean_unsup_32_BETA1pt2_KLanneal_klstart30_annealtime15.csv')
        subtype_z = z[combined]

        # Loading label
        input_path = "input_path = 'data/GDC-PANCAN_"
        label = pd.read_csv(input_path + 'both_samples_tumour_type_digit.tsv', sep='\t', header=0, index_col=0)
        label_array = label['tumour_type'].to_numpy()
        #Number the subtypes (different to normal labels)
        label_array[LumB] = 100
        label_array[Basal] = 101
        label_array[LumA] = 102
        label_array[Her2] = 103
        label_array = label_array[combined]
        y = label_array
        GeneralHelper.plotLatentSpaceTCGATsne(latentSpace=subtype_z, labels=y)

    #Get the samples either side of the distribution as a measure of disentanglement
    if dimensionstandardDevs:
        label = pd.read_csv(input_path + 'both_samples_tumour_type_digit.tsv', sep='\t', header=0, index_col=0)
        label_array = label['tumour_type'].to_numpy()
        z = np.genfromtxt('data/z_before_supervised_32_KLaneaupto2.3KLLOSS14(notopt).csv')
        phenotype = GeneralHelper.processPhenotypeDataForSamples(sample_id)
        conditionone = phenotype['sample_type'] == "Solid Tissue Normal"
        conditiontwo = phenotype['project_id'] == "TCGA-BRCA"
        conditionthree = phenotype['project_id'] == "Primary Tumor"
        combined=np.logical_and.reduce(conditionone, conditiontwo,conditionthree)

        selected_z = z[normal]
        label_array[normal] = 4
        label_array = label_array[combined]

        barChartHelper.compute_diff_capacity_latent_dim(selected_z, subtype_z, label_array,
                                         "CellDifferentiationKidneyKLaneaupto2.pt3NOAbs")

    # Disentangling assessment method
    if heatmaps:
        cancer = label_array == 3
        allElse = label_array != 3
        z = np.genfromtxt('data/z_before_supervised_dim32_beta10.csv')
        cancer_z = z[cancer]
        allElse_z = z[allElse]
        statistics = stats.ttest_ind(cancer_z, allElse_z, axis=0, equal_var=False, nan_policy='propagate')
        pvalue = statistics.pvalue
        # make sure there are no zero values in array- here just set to average
        pvalue[pvalue == 0] = 1.62916828e-001
        sns.set()
        # 4.55993791e-270, 0.001
        log_norm = LogNorm(vmin=4.6e-265, vmax=1.0e-001)
        ax = sns.heatmap([pvalue], norm=log_norm, cbar_kws={"ticks": [0, 1, 10, 1e2, 1e3, 1e4, 1e5]})
        plt.show()

        ShapHelper.saveMostStatisticallySignificantIndex(cancer, allElse, 'data/cancer_z.csv', 'data/normal_z.csv', z)

    if brcaBarCharts:
        label = pd.read_csv(input_path + 'both_samples_tumour_type_digit.tsv', sep='\t', header=0, index_col=0)
        label_array = label['tumour_type'].to_numpy()
        all_z = np.genfromtxt('data/z_before_supervised_loss_32_cap20_beta10_lossclimbedto5noAnneal.csv')
        barChartHelper.compute_diff_capacity_latent_dimBRCA(all_z, all_z, label_array, "Bargraphs_BRCA_vs_NormalTest")
