from random import sample
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from src.tools.utils import get_samples, get_sub_omics_df, read_data, save_splits, get_splits, extract_tumour_type, read_data_pancan, read_data_mlomics
from src.datasets.multi_omics_dataset import MultiOmicsDataset


def prepare_dataset(cohort, sources, n_split, save_split=True, ruche=False, label='PAM50',
                    data_root=None, version='Top', result_directory='results/'):
    if cohort == 'PANCAN':
        omics_df, clinical_df, data_y, lt_samples = read_data_pancan(sources)
        data_y = extract_tumour_type(data_y)
        clinical_df.loc[:, 'tumor_type'] = data_y[:,0]
        le = LabelEncoder().fit(clinical_df.loc[:, 'tumor_type'].values)
        clinical_df.loc[:, 'tumor_type'] = le.transform(clinical_df.loc[:, 'tumor_type'].values)
        ohe = OneHotEncoder(sparse_output=False).fit(clinical_df.loc[:, 'tumor_type'].values.reshape(-1,1))
    elif cohort.startswith('GS-'):
        omics_df, labels = read_data_mlomics(cohort, version=version, data_root=data_root)
        # Keep only requested sources, in CLI order
        missing = [s for s in sources if s not in omics_df]
        if missing:
            raise KeyError(f'Omics sources not found in dataset: {missing}')
        omics_df = {s: omics_df[s] for s in sources}
        lt_samples = list(list(omics_df.values())[0].index)
        if len(labels) != len(lt_samples):
            raise ValueError(
                f'Label/sample size mismatch for {cohort}: '
                f'{len(labels)} labels vs {len(lt_samples)} samples'
            )
        clinical_df = pd.DataFrame({
            label: labels,
            'status': 0,
            'overall_survival': 0.0,
            'OS': 0,
            'OS.time': 0.0,
        }, index=pd.Index(lt_samples))
        ohe = None
        le = None
    else:
        if cohort == 'TCGA-BRCA':
            omics_df, clinical_df , lt_samples = read_data(cohort, sources, label)
        else:
            omics_df, clinical_df , lt_samples = read_data(cohort, sources, label)
            ohe = None
            le = None

    # Write/read splits under result_directory so packaged splits/ are never overwritten
    split_root = os.path.join(result_directory, 'splits')
    os.makedirs(os.path.join(split_root, cohort), exist_ok=True)
    if save_split:
        save_splits(lt_samples, cohort, split_root=split_root)
    samples_train, samples_val, samples_test = get_splits(cohort, n_split, split_root=split_root)

    omics_train = get_sub_omics_df(omics_df, samples_train)
    omics_val = get_sub_omics_df(omics_df, samples_val)
    omics_test = get_sub_omics_df(omics_df, samples_test)

    x_dim = [omics_df[s].shape[1] for s in sources]

    return omics_df, clinical_df, omics_train, omics_val, omics_test, lt_samples, x_dim
