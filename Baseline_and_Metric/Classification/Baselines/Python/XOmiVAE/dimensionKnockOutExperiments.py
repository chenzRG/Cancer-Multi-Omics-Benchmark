import torch
import numpy as np
from torch.nn import functional as F
from sklearn import metrics
import pandas as pd
from exprVAEwithAdditionalFeatures import ExprOmiVAE

import shapExplainerHelper as ShapHelper
import generalHelperFunctions as GeneralHelper


def nullingDimensions(sample_id,expr_df,diseaseCode,chosenTissue,dimension,device='cpu'):
    """
    Tests four different methods to perturb the dimension in the network.
    :param diseaseCode: Required to measure the change in accuracy.
    :param chosenTissue: Required to measure the change in accuracy.
    :param dimension: Chosen dimension to perturb
    """
    input_path = 'DataSources/GDC-PANCAN_'
    vae_model=ExprOmiVAE(input_path,expr_df)


    def runThroughModel(vae_mod,exprTensor,disease_code,dim_number):
        with torch.no_grad():
            _, recon_data, mean, log_var, pred_y = vae_mod(exprTensor)
            #print("mean is")
            #print(mean[:, dim_number])
            pred_y_softmax = F.softmax(pred_y, dim=1)
            _, predicted = torch.max(pred_y_softmax, 1)
            correct_vals = (predicted == disease_code).sum().item()
            accuracy_tumour = correct_vals / len(pred_y) * 100
            return accuracy_tumour,predicted

    def getAccuracyOfModel(disease_code, vae_mod, norm_expr_tensor, tumour_expr_tensor,dim_number):
            accuracy,predictedTumour=runThroughModel(vae_mod,tumour_expr_tensor,disease_code,dim_number)
            print('Cancer tissue accuracy is: {}'.format(accuracy))
            accuracy,predictedNormal = runThroughModel(vae_mod, norm_expr_tensor, 0,dim_number)
            print('Normal tissue accuracy is: {}'.format(accuracy))
            """
            print("Other metrics")
            output_y_normal = np.zeros(predictedNormal.shape)
            output_y_tumour = np.zeros(predictedTumour.shape)
            output_y_tumour[output_y_tumour == 0] = disease_code
            output_y = np.concatenate([output_y_normal, output_y_tumour])
            unique, counts = np.unique(predictedTumour, return_counts=True)
            
            print("Predictions for tumour")
            print(np.asarray((unique, counts)).T)
            unique, counts = np.unique(predictedNormal, return_counts=True)
            print("Predictions for normal tissue")
            print(np.asarray((unique, counts)).T)
            prediction = torch.cat((predictedNormal, predictedTumour), 0)
            precision = metrics.precision_score(output_y, prediction, average='weighted')
            print('Precision: {}'.format(precision))
            recall = metrics.recall_score(output_y, prediction, average='weighted')
            print('Recall: {}'.format(recall))
            f1Score = metrics.f1_score(output_y, prediction, average='weighted')
            print('f1Score: {}'.format(f1Score))
            """


    def changeWeightandBias(weightZero=False, weightOne=False, weightReverse=False, changeOutput=True, dim_number=0, vae_mod=vae_model,
                            disease_code=1, norm_expr_tensor=0, tumour_expr_tensor=0):
        vae_mod.load_state_dict(torch.load('DataSources/vae_saved_model(original).pt', map_location=torch.device('cpu')))
        vae_mod.eval()
        with torch.no_grad():

            if weightZero == True:
                vae_mod.e_fc4_mean[0].weight[dim_number] = 0.
                #vae_mod.e_fc4_mean[0].weight[dim_number, :] = 0.

            elif weightOne == True:
                vae_mod.e_fc4_mean[0].weight[dim_number, :] = 1.
                vae_mod.e_fc4_mean[0].bias[dim_number] = 200.
            elif weightReverse == True:
                vae_mod.e_fc4_mean[0].weight[dim_number, :] = (vae_mod.e_fc4_mean[0].weight[dim_number, :]) * -1
            # this is for changing the output of the dimension

            elif changeOutput == True:
                vae_mod.load_state_dict(
                    torch.load('DataSources/vae_saved_model(original).pt',
                               map_location=torch.device('cpu')))
                vae_mod.eval()
        getAccuracyOfModel(disease_code, vae_mod, norm_expr_tensor, tumour_expr_tensor,dim_number)


    def extract_relevantExpr(expr_df, condition):
        expr_df_T = expr_df.T
        relevantExpr = expr_df_T[condition].T
        # tumour_expr = tumour_expr.sample(n=100, axis=1)
        relevantExpr_tensor=GeneralHelper.addToTensor(relevantExpr,device)
        return relevantExpr_tensor


    def nullingDimensions(disease_code, disease_string, dim_number):
        phenotype = pd.read_csv('DataSources/GDC-PANCAN.basic_phenotype.tsv', sep='\t', header=0, index_col=0)
        phenotype = phenotype.T[sample_id]
        phenotype = phenotype.T
        conditionone = phenotype['sample_type'] == "Primary Tumor"
        conditiontwo = phenotype['project_id'] == disease_string
        conditionthree = phenotype['sample_type'] == "Solid Tissue Normal"
        conditionaltumour = np.logical_and(conditionone, conditiontwo)
        conditionalnormal = np.logical_and(conditionthree, conditiontwo)
        norm_expr_tensor = extract_relevantExpr(expr_df, conditionalnormal)
        tumour_expr_tensor = extract_relevantExpr(expr_df, conditionaltumour)
        #saved pickle file of omiVAE we are using
        vae_model.load_state_dict(torch.load('DataSources/vae_saved_model(original).pt', map_location=torch.device('cpu')))
        vae_model.eval()

        print("normal before")
        getAccuracyOfModel(disease_code, vae_model, norm_expr_tensor=norm_expr_tensor,
                           tumour_expr_tensor=tumour_expr_tensor,dim_number=dim_number)
        #print("weight of model before")
        #print(vae_model.e_fc4_mean[0].weight[dim_number, :])
        
        print("weight set to zero")
        changeWeightandBias(weightZero=True, dim_number=dim_number, disease_code=disease_code,
                            norm_expr_tensor=norm_expr_tensor, tumour_expr_tensor=tumour_expr_tensor)
       
        print("weight set to 1")
        changeWeightandBias(weightOne=True, dim_number=dim_number, vae_mod=vae_model, disease_code=disease_code,
                            norm_expr_tensor=norm_expr_tensor, tumour_expr_tensor=tumour_expr_tensor)
       
        print("weight reversed")
        changeWeightandBias(weightReverse=True, dim_number=dim_number, vae_mod=vae_model, disease_code=disease_code,
                            norm_expr_tensor=norm_expr_tensor, tumour_expr_tensor=tumour_expr_tensor)

        print("change output of model")
        from exprVAEwithAdditionalFeatures import ExprOmiVAE
        input_path="DataSources/GDC-PANCAN_"
        vae_model_two = ExprOmiVAE(input_path=input_path, expr_df=expr_df, latent_dim=128, knockingOut=True, dimension_number=dim_number)

        vae_model_two.load_state_dict(
            torch.load('DataSources/vae_saved_model(original).pt',
                       map_location=torch.device('cpu')))

        vae_model_two.eval()
        changeWeightandBias(changeOutput=True, dim_number=dim_number, vae_mod=vae_model_two, disease_code=disease_code,
                            norm_expr_tensor=norm_expr_tensor, tumour_expr_tensor=tumour_expr_tensor)


    nullingDimensions(disease_code=diseaseCode, disease_string=chosenTissue, dim_number=dimension)
