import torch
from torch.utils.data import Dataset
import numpy as np


class MultiOmicsDataset(Dataset):
    """
    Load multi-omics data
    """

    def __init__(self, omics_df, clinical_df, lt_samples, label, event, surv_time):
        self.omics_df = omics_df
        self.clinical_df = clinical_df
        self.lt_samples = lt_samples
        self.label = label
        self.event = event
        self.surv_time = surv_time

    def __len__(self):
        return len(self.lt_samples)

    def __getitem__(self, index):
        sample = self.lt_samples[index]
        omics_data = []
        for source, omic_df in zip(self.omics_df.keys(), self.omics_df.values()):
            omic_line = omic_df.loc[sample, :].values
            omic_line = omic_line.astype(np.float32)
            omic_line_tensor = torch.Tensor(omic_line)
            omics_data.append(omic_line_tensor)
            if self.label:
                label = self.clinical_df.loc[sample, self.label]
            else:
                label = 0
            os_time = int(self.clinical_df.loc[sample, self.surv_time])
            os_event = int(self.clinical_df.loc[sample, self.event])  
        return omics_data, label, os_time, os_event

    def return_samples(self):
        return self.lt_samples


class PANCANDataset(Dataset):
    """
    Load multi-omics data
    """

    def __init__(self, omics_df, labels, cohort):
        self.omics_df = omics_df
        self.labels = labels
        self.cohort = cohort

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, index):
        omics_data = []
        for source, omic_df in zip(self.omics_df.keys(), self.omics_df.values()):
            omic_line = omic_df.iloc[index, :].values
            omic_line = omic_line.astype(np.float32)
            omic_line_tensor = torch.Tensor(omic_line)
            omics_data.append(omic_line_tensor)
        label = self.labels[index]
        return omics_data, label