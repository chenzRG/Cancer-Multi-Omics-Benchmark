from random import sample
import numpy as np
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from src.tools.utils import get_samples, get_sub_omics_df, read_data, save_splits, get_splits, extract_tumour_type, read_data_pancan
from src.datasets.multi_omics_dataset import MultiOmicsDataset


def prepare_dataset(cohort, sources, n_split, save_split=True, ruche=False, label='PAM50'):
    if cohort == 'PANCAN':
        omics_df, clinical_df, data_y, lt_samples = read_data_pancan(sources)
        data_y = extract_tumour_type(data_y)
        clinical_df.loc[:, 'tumor_type'] = data_y[:,0]
        le = LabelEncoder().fit(clinical_df.loc[:, 'tumor_type'].values)
        clinical_df.loc[:, 'tumor_type'] = le.transform(clinical_df.loc[:, 'tumor_type'].values)
        ohe = OneHotEncoder(sparse=False).fit(clinical_df.loc[:, 'tumor_type'].values.reshape(-1,1))
    else:
        if cohort == 'TCGA-BRCA':
            omics_df, clinical_df , lt_samples = read_data(cohort, sources, label)
        else:
            omics_df, clinical_df , lt_samples = read_data(cohort, sources, label)
            ohe = None
            le = None
    if save_split:
        save_splits(lt_samples, cohort)
    samples_train, samples_val, samples_test = get_splits(cohort, n_split)
    
    omics_train = get_sub_omics_df(omics_df, samples_train)
    omics_val = get_sub_omics_df(omics_df, samples_val)
    omics_test = get_sub_omics_df(omics_df, samples_test)


    x_dim = [omics_df[omic_source].shape[1] for omic_source in omics_df.keys()]

    return omics_df, clinical_df, omics_train, omics_val, omics_test, lt_samples, x_dim
