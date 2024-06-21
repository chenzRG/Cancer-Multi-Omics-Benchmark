import torch
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from exprVAEwithAdditionalFeatures import ExprOmiVAE
import generalHelperFunctions as GeneralHelper
import shapExplainerHelper as ShapHelper
from dimensionKnockOutExperiments import nullingDimensions
from torch.nn import functional as F
import sys
sys.path.insert(0, '~/PycharmProjects/shapLundberg')
import shapLundberg
from shapLundberg import shap
from sklearn import metrics
from scipy.stats import skew
#import gffpandas.gffpandas as gffpd

def omiShapExplainer(sample_id,label_array,expr_df,dimension=0,tumourID=0,device='cpu',LatentSpaceExplain=False, statisticsTest=False, NormalvsTumourExplain=False,
               NormalvsTumourInterimExplain=False, TestingNullingDimensions=False,
               saveLatentSpace=False, statisticsTestEvenGender=False,histogram=False,NormalvsTumourDimensionExplain=False, knockoutGenes=False,idProvided=True,tumourName="TCGA-KIRC",file_name="place_holder",knockOutNumber=50,mean=False,toZero=True,TopGenes=True,skewValues=False,knockOutAcrossTumours=False,measureGenesinCommon=False):
    """
    :param tumourName: String version of TCGA tumour e.g. "TCGA-KIRC" we would like to explain
    :param NormalvsTumourExplain: Explain predictions (from supervised OmiVAE). Pass in tumourName to explain.
    :param NormalvsTumourInterimExplain: Explain most important dimensions (from supervised OmiVAE). Pass in tumourName to explain.
    :param NormalvsTumourDimensionExplain: Explain a specified dimension (either after supervised or unsupervised training). Pass in dimension and tumourName.
    :param statisticsTest: Find the most important dimension in the unsupervised part of the model. Adjust within the code for which latent space (z) file to load in and which factrs we are investigating.
    :param LatentSpaceExplain: Explain the mean latent variables. Need to specifiy dimension. Dimension found by running statisticsTest. Adjust within section
    :param tumourID: TCGA tumour code (see tumour type digit file)  that we want to explain prediction. Required for testing nulling dimensions.
    :return:
    """
    #This class has combined the different analysis' of the Deep SHAP values we conducted.
    #SHAP reference: Lundberg et al., 2017: http://papers.nips.cc/paper/7062-a-unified-approach-to-interpreting-model-predictions.pdf


    input_path = 'DataSources/GDC-PANCAN_'
    vae_model = ExprOmiVAE(input_path, expr_df)

    #load vae trained model
    vae_model.load_state_dict(torch.load('DataSources/vae_saved_model(original).pt',map_location=torch.device('cpu')))

    #Get the most statistically significant dimensions: use for explaining the UNSUPERVISED model (just VAE)
    if NormalvsTumourExplain:
        """
            Explain the most important genes for a specific tumour, when using the normal tissue as a background
            Note: adjusting to use a random training sample as the background is simple, use the 
            ShapHelper.randomTrainingSample(expr) function.
            """
        phenotype = GeneralHelper.processPhenotypeDataForSamples(sample_id)
        print("cancer type is")
        print(tumourName)

        conditionone=phenotype['sample_type'] == "Primary Tumor"
        conditiontwo=phenotype['project_id'] == tumourName
        conditionthree=phenotype['sample_type'] == "Solid Tissue Normal"

        conditionaltumour=np.logical_and(conditionone,conditiontwo)

        conditionalnormal = np.logical_and(conditiontwo, conditionthree)
        """
        #recommended to have 100 samples, but in some cases it is not possible to sample 100 for a tissue type, especially for normal tissue
        print("count of normal tissue")
        print(np.count_nonzero(conditionalnormal))
        if np.count_nonzero(conditionalnormal) < 10:
            normal_number=np.count_nonzero(conditionalnormal)
        else:
            normal_number=10
        
        normal_expr = ShapHelper.splitExprandSample(condition=conditionalnormal, sampleSize=10, expr=expr_df)
    
        if np.count_nonzero(conditionaltumour) < 10:
            conditionalcount = np.count_nonzero(conditionaltumour)
        else:
            conditionalcount = 10
        """

        normal_expr = ShapHelper.randomTrainingSample(expr_df, 10)
        tumour_expr = ShapHelper.splitExprandSample(condition=conditionaltumour, sampleSize=10, expr=expr_df)

        # put on device as correct datatype
        background = GeneralHelper.addToTensor(expr_selection=normal_expr, device=device)
        male_expr_tensor = GeneralHelper.addToTensor(expr_selection=tumour_expr, device=device)
        GeneralHelper.checkPredictions(male_expr_tensor, vae_model)

        #output number is based from OmiVAE forward function. 4= predicted y value.
        e = shap.DeepExplainer(vae_model, background, outputNumber=4)
        print("calculating shap values")
        shap_values_female = e.shap_values(male_expr_tensor, ranked_outputs=1)
        genes,chrom, ensg =ShapHelper.getGenes(expr_df)

        shap.summary_plot(shap_values_female[0][0],features=tumour_expr.T,feature_names=genes)
        print("calculated shap values")
        most_imp  = ShapHelper.getTopShapValues(shap_vals=shap_values_female, numberOfTopFeatures=50,
                                                          expr_df=expr_df, ranked_output=True, cancerType=tumourName,absolute=True)

    if NormalvsTumourDimensionExplain:
        """
            Explain the most important genes for a specific dimension for a tumour.            
            Note: adjusting to use a random training sample as the background is simple, use the 
            ShapHelper.randomTrainingSample(expr) function.
            Here we are explaining the mean output (outputNumber=2).
            """
        phenotype = GeneralHelper.processPhenotypeDataForSamples(sample_id)
        conditionone=phenotype['sample_type'] == "Primary Tumor"
        conditiontwo=phenotype['project_id'] == tumourName
        conditionthree=phenotype['sample_type'] == "Solid Tissue Normal"

        conditionaltumour=np.logical_and(conditionone,conditiontwo)
        conditionalnormal = np.logical_and(conditiontwo, conditionthree)
        #recommended to have 100 samples, but in some cases it is not possible to sample 100 for a tissue type, especially for normal tissue
        normal_expr = ShapHelper.splitExprandSample(condition=conditionalnormal, sampleSize=2, expr=expr_df)
        tumour_expr = ShapHelper.splitExprandSample(condition=conditionaltumour, sampleSize=2, expr=expr_df)
        # put on device as correct datatype
        background = GeneralHelper.addToTensor(expr_selection=normal_expr, device=device)
        male_expr_tensor = GeneralHelper.addToTensor(expr_selection=tumour_expr, device=device)
        GeneralHelper.checkPredictions(background, vae_model)

        #output number 2 is the mean (latent space)
        #I edited the SHAP library to allow a dimension to be passed in and edited the output value. important that explainLatentSpace=True.
        e = shap.DeepExplainer(vae_model, background, dim=dimension, outputNumber=2, explainLatentSpace=True)
        #e = shap.DeepExplainer(vae_model, background)
        #as we are explaining the genes for a dimension (unlike previously where we had to choose which prediction we were explaining) we do not need to specifiy to rank the outputs
        shap_values_female = e.shap_values(male_expr_tensor)

        most_imp = ShapHelper.getTopShapValues(shap_vals=shap_values_female, numberOfTopFeatures=50,
                                                          expr_df=expr_df, ranked_output=False)

    #Note: this can be easily adjusted to explain any two samples; simply change 'normal_expr' and 'tumour_expr' to the chosen samples
    if NormalvsTumourInterimExplain:
        """
        Explain the interim layer; pass in the interim layer and then it explains the top dimensions.
        """
        phenotype = GeneralHelper.processPhenotypeDataForSamples(sample_id)
        conditionone = phenotype['sample_type'] == "Primary Tumor"
        conditiontwo = phenotype['project_id'] == tumourName
        conditionthree = phenotype['sample_type'] == "Solid Tissue Normal"
        conditionaltumour = np.logical_and(conditionone, conditiontwo)
        conditionalnormal = np.logical_and(conditionthree, conditiontwo)


        # Recommended to have 100 samples, but in some cases it is not possible to sample 100 for a tissue type, especially for normal tissue
        #change background here; e.g. random training sample or normal tissue
        normal_expr = ShapHelper.randomTrainingSample(expr_df, 50)
        #normal_expr = ShapHelper.splitExprandSample(condition=conditionone, sampleSize=50, expr=expr_df)
        tumour_expr = ShapHelper.splitExprandSample(condition=conditionaltumour, sampleSize=50, expr=expr_df)
        # Ensure tensor is on device and has the correct data type
        background = GeneralHelper.addToTensor(expr_selection=normal_expr, device=device)
        male_expr_tensor = GeneralHelper.addToTensor(expr_selection=tumour_expr, device=device)

        GeneralHelper.checkPredictions(background, vae_model)
        GeneralHelper.checkPredictions(male_expr_tensor, vae_model)
        #output number is based from OmiVAE forward function. 4= predicted y value
        #this should be c_fc1 for autoencoded layer
        e = shap.DeepExplainer((vae_model, vae_model.c_fc1), background, outputNumber=4)

        #explains layer before mean (512 dimensions)
        #e = shap.DeepExplainer((vae_model, vae_model.e_fc4_mean), background)
        shap_values_female = e.shap_values(male_expr_tensor, ranked_outputs=1)
        # Here look at the numbers to left (they should range 1 to no. of dimensions) and this is the most important dimension
        most_imp, least_imp = ShapHelper.getTopShapValues(shap_vals=shap_values_female, numberOfTopFeatures=50,
                                                          expr_df=expr_df, ranked_output=True)


    if LatentSpaceExplain:
        """
            This requires a known dimension to explain (found using the statistical significance above),
            This is for gender, subtype and pathway latent space explanations for the unsupervised part of the model.
            """
        #Example of how to obtain the relevant samples for the subtypes we would like to analyse (source: https://www.sciencedirect.com/science/article/pii/S0092867418303593)
        subtypeConditionalOne, subtypeConditionalTwo = GeneralHelper.processSubtypeSamples(sample_id,subtypeOne="Basal",subtypeTwo="LumB")

        #Example of how to obtain the relevant samples for the pathways we would like to analyse (source: https://www.sciencedirect.com/science/article/pii/S0092867418303593)
        logicalOne, logicalTwo = ShapHelper.pathwayComparison(sample_id=sample_id,pathway='RTK RAS')

        #Example of how to obtain the relevant samples for the different genders to explain the latent space
        femaleCondition, maleCondition=ShapHelper.splitForGenders(sample_id=sample_id,)
        ShapHelper.printConditionalSelection(maleCondition)

        #Change condition to one of the relevant conditions above
        female_expr = ShapHelper.splitExprandSample(condition=femaleCondition,sampleSize=50,expr=expr_df)
        male_expr = ShapHelper.splitExprandSample(condition=maleCondition, sampleSize=50, expr=expr_df)

        #put on device as correct datatype
        background= GeneralHelper.addToTensor(expr_selection=female_expr,device=device)
        male_expr_tensor = GeneralHelper.addToTensor(expr_selection=male_expr, device=device)
        GeneralHelper.checkPredictions(background, vae_model)
        e = shap.DeepExplainer(vae_model, background,outputNumber=2,dim=dimension)
        # If explaining the z/mu dimension then don't need ranked outputs like we used before (as only one output from the model)
        shap_values_female = e.shap_values(male_expr_tensor)
        most_imp, least_imp=ShapHelper.getTopShapValues(shap_vals=shap_values_female, numberOfTopFeatures=50, expr_df=expr_df, ranked_output=False)

    if saveLatentSpace:
        expr_tensor=GeneralHelper.addToTensor(expr_df,device)
        vae_model.load_state_dict(torch.load('data/beta15.pt', map_location=torch.device('cpu')))
        GeneralHelper.saveLatentSpace(vae_model,expr_tensor)

    if statisticsTest:
        # z latent space; explain beta-vae and normal vae latent space (can also feed in the mean here)
        z = np.genfromtxt('data/z_before_supervised_loss32_ForBeta10AnnealingKlLoss8.csv')
        #Split for gender here but can also split for subtype if we want to find the most important dimension for a subtype
        femaleCondition, maleCondition = ShapHelper.splitForGenders(sample_id=sample_id, )

        # Example for adapting this code to measure the statistically significant dimensions in a subtype
        # subtypeConditionalOne, subtypeConditionalTwo = GeneralHelper.processSubtypeSamples(sample_id, subtypeOne="Basal",
        #                                                                                        subtypeTwo="LumB")

        female_z = z[femaleCondition]
        male_z = z[maleCondition]
        statistics = stats.ttest_ind(female_z, male_z, axis=0, equal_var=False, nan_policy='propagate')
        stat_lowest_index = np.argmin(statistics.pvalue)
        female_z_column = female_z[:, stat_lowest_index]
        male_z_column = male_z[:, stat_lowest_index]
        np.savetxt('data/brca_z.csv', female_z_column)
        np.savetxt('data/normal_z.csv', male_z_column)

    # Method to sample the same amount to ensure gender split is even
    if statisticsTestEvenGender:
        # Get the z vector we would like to explain
        z = np.genfromtxt('data/z_before_supervised_loss32_ForBeta10AnnealingKlLoss8.csv')
        phenotype = GeneralHelper.processPhenotypeDataForSamples(sample_id)
        # tumour id's were chosen as they were typically non-gender specific cancers
        femaleCondition, maleCondition = ShapHelper.multipleSampling(label_array=label_array, phenotype=phenotype,
                                                                     tumour_id_one=13,
                                                                     tumour_id_two=18, tumour_id_three=10,
                                                                     tumour_id_four=5, tumour_id_five=6,
                                                                     sample_number=50)
        ShapHelper.saveMostStatisticallySignificantIndex(conditionOne=femaleCondition, conditionTwo=maleCondition,
                                                     fileNameOne="female_z", fileNameTwo="male_z",z=z)

    if histogram:
        feature_importance = pd.read_csv('data/luaddimensions.csv', header=0, index_col=0)
        plt.style.use('seaborn-white')
        hist = feature_importance.hist(column='feature_importance_vals',normed=False, alpha=0.5,
         histtype='stepfilled', color='skyblue',
         edgecolor='none')
        plt.title('SHAP value histogram')
        plt.xlabel('SHAP_value')
        plt.ylabel('Frequency')
        plt.yscale('log')
        plt.savefig("SHAPvalueLoghistogram.png", dpi=1500)
        plt.show()

    if TestingNullingDimensions:
        #Example of testing four different dimensions
        print("dim 42")
        nullingDimensions(sample_id=sample_id,expr_df=expr_df,diseaseCode=tumourID, chosenTissue=tumourName, dimension=42)
        print("dim 125")
        nullingDimensions(sample_id=sample_id,expr_df=expr_df,diseaseCode=tumourID, chosenTissue=tumourName, dimension=125)

    def runThroughModel(vae_mod, exprTensor, disease_code):
        with torch.no_grad():
            _, recon_data, mean, log_var, pred_y = vae_mod(exprTensor)
            pred_y_softmax = F.softmax(pred_y, dim=1)
            _, predicted = torch.max(pred_y_softmax, 1)
            correct_vals = (predicted == disease_code).sum().item()
            accuracy_tumour = correct_vals / len(pred_y) * 100
            return accuracy_tumour, predicted

    def runThroughModelFullExp(vae_mod, exprTensor):
        with torch.no_grad():
            _, recon_data, mean, log_var, pred_y = vae_mod(exprTensor)
            pred_y_softmax = F.softmax(pred_y, dim=1)
            _, predicted = torch.max(pred_y_softmax, 1)
            return predicted

    def getAccuracyOfModel(vae_mod, disease_code,norm_expr_tensor, tumour_expr_tensor):
        accuracy, predictedTumour = runThroughModel(vae_mod, tumour_expr_tensor, disease_code)
        print('Cancer tissue accuracy is: {}'.format(accuracy))
        accuracy, predictedNormal = runThroughModel(vae_mod, norm_expr_tensor, 0)
        print('Normal tissue accuracy is: {}'.format(accuracy))
        unique, counts = np.unique(predictedTumour, return_counts=True)

        print("Predictions for tumour")
        print(np.asarray((unique, counts)).T)


    def extract_relevantExpr(expr_df, condition):
        expr_df_T = expr_df.T
        relevantExpr = expr_df_T[condition].T
        # tumour_expr = tumour_expr.sample(n=100, axis=1)
        relevantExpr_tensor=GeneralHelper.addToTensor(relevantExpr,device)
        return relevantExpr_tensor

    def geneKnockOutExpr(expr_df,conditionaltumour,file_name,knockOutNumber,idProvided,TopGenes):
        shap_genes=pd.read_csv(file_name, sep=',', header=0, index_col=1)
        #knock out the bottom genes- sort by lowests to highest
        if not TopGenes:
            shap_genes.sort_values(by=['feature_importance_vals'], ascending=True, inplace=True)
            print("rearranges bottom to top")
        #get top number of genes to knock out
        gencode_ids = pd.read_csv('DataSources/gencodev22annotationgene.tsv', sep='\t', header=0, index_col=1)
        # merge includes duplicate genes!
        if idProvided:
            shapwithID = pd.merge(shap_genes, gencode_ids, on='id', how='inner')
            mostImp_shap_values = shapwithID.head(knockOutNumber)
            ensembleCodes = mostImp_shap_values.iloc[:, 2]
        else:
            shapwithID = pd.merge(shap_genes, gencode_ids, on='gene', how='inner')
            mostImp_shap_values = shapwithID.head(knockOutNumber)
            ensembleCodes = mostImp_shap_values.iloc[:, 3]

        expr_df=expr_df.T
        """
        if mean:
            newValuesForGenes = expr_df.loc[expr_df.index.intersection(ensembleCodes)].mean(axis=1)
            newvalues = np.resize(newValuesForGenes, (np.count_nonzero(conditionaltumour), newValuesForGenes.size))
            expr_df.loc[expr_df[conditionaltumour].index, expr_df.columns.intersection(ensembleCodes)] = newvalues
        """
        expr_df.loc[expr_df[conditionaltumour].index, expr_df.columns.intersection(ensembleCodes)] = 0
        print("values after")
        print(expr_df.loc[expr_df[conditionaltumour].index, expr_df.columns.intersection(ensembleCodes)])

        relevantExpr = expr_df[conditionaltumour].T
        # tumour_expr = tumour_expr.sample(n=100, axis=1)
        relevantExpr_tensor = GeneralHelper.addToTensor(relevantExpr, device)
        expr_df = expr_df.T
        return relevantExpr_tensor,expr_df

    if knockoutGenes:
        vae_model.eval()
        phenotype = pd.read_csv('DataSources/GDC-PANCAN.basic_phenotype.tsv', sep='\t', header=0, index_col=0)
        phenotype = phenotype.T[sample_id]
        phenotype = phenotype.T
        conditionone = phenotype['sample_type'] == "Primary Tumor"
        conditiontwo = phenotype['project_id'] == tumourName
        conditionthree = phenotype['sample_type'] == "Solid Tissue Normal"
        conditionaltumour = np.logical_and(conditionone, conditiontwo)
        conditionalnormal = np.logical_and(conditionthree, conditiontwo)
        norm_expr_tensor = extract_relevantExpr(expr_df, conditionalnormal)
        tumour_expr_tensor = extract_relevantExpr(expr_df, conditionaltumour)

        print("accuracy before")
        ### accuracy of just tumour
        getAccuracyOfModel(vae_model,disease_code=tumourID, norm_expr_tensor=norm_expr_tensor,
                           tumour_expr_tensor=tumour_expr_tensor)
        ##total accuracy
        expr_tensor = GeneralHelper.addToTensor(expr_df, device)
        predicted = runThroughModelFullExp(vae_model, expr_tensor)
        unique, counts = np.unique(predicted, return_counts=True)

        print("Predictions for ALL")
        print(np.asarray((unique, counts)).T)
        #accuracy = metrics.accuracy_score(label_array, predicted, average='weighted')
        #print('accuracy: {}'.format(accuracy))

        print("accuracy after")
        tumour_expr_tensor=geneKnockOutExpr(expr_df,conditionaltumour,file_name,knockOutNumber,idProvided)

        getAccuracyOfModel(vae_model, disease_code=tumourID, norm_expr_tensor=norm_expr_tensor,
                           tumour_expr_tensor=tumour_expr_tensor)
        ##total accuracy
        expr_tensor = GeneralHelper.addToTensor(expr_df, device)
        predicted = runThroughModelFullExp(vae_model, expr_tensor)
        unique, counts = np.unique(predicted, return_counts=True)
        print("Predictions for ALL")
        print(np.asarray((unique, counts)).T)

        #accuracy = metrics.accuracy_score(label_array, predicted, average='weighted')
        #print('accuracy: {}'.format(accuracy))
    def getPredictions(expr_df,device,vae_model):
        expr_tensor = GeneralHelper.addToTensor(expr_df, device)
        predicted = runThroughModelFullExp(vae_model, expr_tensor)
        unique, countsBefore = np.unique(predicted, return_counts=True)
        print(np.asarray((unique, countsBefore)).T)
        return unique, countsBefore

    def getConditionalTumour(sample_id,tumourName):
        tumourName= "TCGA-"+tumourName
        phenotype = pd.read_csv('DataSources/GDC-PANCAN.basic_phenotype.tsv', sep='\t', header=0, index_col=0)
        phenotype = phenotype.T[sample_id]
        phenotype = phenotype.T
        conditionone = phenotype['sample_type'] == "Primary Tumor"
        conditiontwo = phenotype['project_id'] == tumourName
        conditionaltumour = np.logical_and(conditionone, conditiontwo)
        return conditionaltumour

    def runThroughKnockOutAndSave(expr_df,device,vae_model,conditionaltumour,tumourType,count,TopGenes):
         uniqueBefore, countsBefore = getPredictions(expr_df, device, vae_model)
         tumour_expr_tensor, expr_df = geneKnockOutExpr(expr_df, conditionaltumour, file_name, count, idProvided,TopGenes)
         uniqueAfter, countsAfter = getPredictions(expr_df, device, vae_model)

         transcriptCountsBefore = np.asarray((uniqueBefore, countsBefore)).T
         transcriptCountsAfter =   np.asarray((uniqueAfter, countsAfter)).T

         transcriptCountsBeforeDF = pd.DataFrame(transcriptCountsBefore,
                           columns=['unique', 'countsBefore'])
         transcriptCountsAfterDF = pd.DataFrame(transcriptCountsAfter,
                                                 columns=['unique', 'countsAfter'])
         #left: useo nly keys from left frame
         countsMerges = pd.merge(transcriptCountsBeforeDF, transcriptCountsAfterDF[['unique','countsAfter']], on='unique', how='left')

         countsMerges=countsMerges.fillna(0)
         subtract=(countsMerges['countsAfter'] - countsMerges['countsBefore']).values
         countsMerges.insert(3, "Subtraction", subtract, True)
         accuracy=((countsMerges['countsAfter'] / countsMerges['countsBefore'])*100).values
         countsMerges.insert(4,"Accuracy", accuracy, True)


         if TopGenes:
             fileAddition = "_top_SHAP_values"
         else:
             fileAddition = "_everything_but_top_SHAP_values"
         countsMerges.to_csv("TCGA-" + tumourType + "knockout_"+str(count)+fileAddition)

    if knockOutAcrossTumours:
        TumourID = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                    'KIRP', 'LAML','LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD',
                    'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        TumourID = ['ACC']
        vae_model.eval()
        print("tumour is")
        print(TumourID)
        #knockOutNumbers = [0.2,0.5,1,2,5,6]
        knockOutNumbers = [6]
        for tumourType in TumourID:
            conditionaltumour = getConditionalTumour(sample_id, tumourType)
            file_name = 'SHAP_with_ID/TCGA-' + tumourType + '_feature_imp_secondtry'
            for numberOfGenes in knockOutNumbers:
                if numberOfGenes==6:
                    #newknockOutNumbers=[0.2,0.5,1,5]
                    newknockOutNumbers = [0.2]
                    for i in newknockOutNumbers:
                        print("Knocking out everything but the amount of genes numbering")
                        count = round(numberOfGenes / 100 * expr_df.shape[0])
                        count=expr_df.shape[0]-count
                        print("count is")
                        print(count)
                        TopGenes=False
                        runThroughKnockOutAndSave(expr_df, device, vae_model, conditionaltumour, tumourType,count,TopGenes)
                print("count is")
                count=round(numberOfGenes/100*expr_df.shape[0])

                print(count)
                TopGenes = True
                runThroughKnockOutAndSave(expr_df,device,vae_model,conditionaltumour,tumourType,count,TopGenes)

    if measureGenesinCommon:
    #go through Top 1% of genes in SHAP file
        TumourID = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                    'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD',
                    'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        ShapVals=[]
        print("count is")
        #top 0.1%
        count = round(0.1 / 100 * expr_df.shape[0])
        for x in TumourID:
            feature_importance = pd.read_csv('SHAP_with_ID/TCGA-' + x + '_feature_imp_secondtry', header=0, index_col=0)
            feature_importance = feature_importance.head(count)
            ###ACTUALLY WANT TOP 1% INDEX OF ID HERE!!!!!
            ShapVals.append(feature_importance['id'].values)
        #list_of_tuples = list(zip(TumourID, ShapVals))
        #vdf = pd.DataFrame(list_of_tuples,columns=['TumourID', 'SHAP_vals'])
        flattened_list = []
        # flatten the lis
        for x in ShapVals:
            for y in x:
                flattened_list.append(y)
        unique, counts = np.unique(flattened_list, return_counts=True)

        print("Predictions for tumour")
        transcriptCounts=np.asarray((unique, counts)).T
        #transcriptCountsPd=pd.DataFrame(data=transcriptCounts[:,1],index=transcriptCounts[:,0],columns=['count'])
        print(transcriptCounts)
        gencode_ids = pd.read_csv('DataSources/gencodev22annotationgene.tsv', sep='\t', header=0)
        transcriptCounts=pd.DataFrame(data=transcriptCounts[:, :], columns=['id', 'count'])
        # merge includes duplicate genes!
        flattened_list_with_gene = pd.merge(transcriptCounts, gencode_ids[['gene','id']], on='id', how='inner')
        flattened_list_with_gene['count'] = pd.to_numeric(flattened_list_with_gene['count'])
        flattened_list_with_gene.sort_values(by=['count'], ascending=False, inplace=True)
        print(flattened_list_with_gene)
        flattened_list_with_gene.to_csv("topGenesAllTumoursForTop_"+str(count)+"genes.csv",columns = ['gene', 'count'])






    if skewValues:
        TumourID=['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML',
                  'LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM']
        Skew=[]
        for x in TumourID:
            feature_importance = pd.read_csv('SHAP_feature_importance/TCGA-'+x+'_feature_imp', header=0, index_col=0)
            Skew.append(skew(feature_importance['feature_importance_vals']))
            print(x)
        list_of_tuples = list(zip(TumourID, Skew))
        df = pd.DataFrame(list_of_tuples,
                          columns=['TumourID', 'Skew'])
        df.to_csv("SHAP_value_skew_by_tumourID")