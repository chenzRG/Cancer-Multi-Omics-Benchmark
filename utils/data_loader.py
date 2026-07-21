"""
MLOmics unified data loader.

Usage:
    from utils.data_loader import MLOmicsLoader

    loader = MLOmicsLoader('/path/to/Main_Dataset')

    # Classification
    omics = loader.load_omics('GS-BRCA', version='Top')
    X     = loader.get_X('GS-BRCA', version='Top')          # (samples × features)
    y     = loader.load_labels('GS-BRCA')

    # Clustering
    omics  = loader.load_omics('ACC', version='Aligned')
    surv   = loader.load_survival('ACC')

    # Imputation
    omics  = loader.load_omics('Imp-BRCA', version='Top')
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple


OMICS_TYPES = ['mRNA', 'miRNA', 'CNV', 'Methy']

CLASSIFICATION_DATASETS = ['BRCA', 'COAD', 'GBM', 'LGG', 'OV', 'Pan-cancer']
GS_DATASETS             = ['BRCA', 'COAD', 'GBM', 'LGG', 'OV']
CLUSTERING_DATASETS     = ['ACC', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'THCA', 'THYM']
IMPUTATION_DATASETS     = ['Imp-BRCA', 'Imp-COAD', 'Imp-GBM', 'Imp-LGG', 'Imp-OV']


class MLOmicsLoader:
    """
    Load MLOmics datasets in a unified format.

    All omics files are returned as pandas DataFrames with shape
    (features × samples).  Call .T or use get_X() to get (samples × features).

    Parameters
    ----------
    data_root : str
        Path to the Main_Dataset directory.
        e.g. '/path/to/MLOmics/Main_Dataset'
    """

    def __init__(self, data_root: str):
        if not os.path.isdir(data_root):
            raise FileNotFoundError(f'data_root not found: {data_root}')
        self.data_root = data_root

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def load_omics(
        self,
        dataset: str,
        version: str = 'Top',
        omics_types: Optional[List[str]] = None,
    ) -> Dict[str, pd.DataFrame]:
        """
        Load one or more omics modalities for a dataset.

        Returns
        -------
        dict
            {omics_name: DataFrame(features × samples)}
        """
        if omics_types is None:
            omics_types = OMICS_TYPES
        folder = self._resolve_folder(dataset, version)
        result = {}
        for omics in omics_types:
            fname = self._omics_filename(dataset, omics, version)
            path = os.path.join(folder, fname)
            if not os.path.exists(path):
                raise FileNotFoundError(f'Omics file not found: {path}')
            result[omics] = pd.read_csv(path, index_col=0, low_memory=False)
        return result

    def load_labels(self, dataset: str, version: str = 'Original') -> np.ndarray:
        """
        Load integer class labels for a classification dataset.

        Returns
        -------
        np.ndarray, shape (n_samples,)
        """
        folder = self._resolve_folder(dataset, version)
        cancer = self._canonical_cancer(dataset)
        if dataset == 'Pan-cancer':
            label_fname = 'Pan-cancer_label_num.csv'
        else:
            label_fname = f'{cancer}_label_num.csv'
        path = os.path.join(folder, label_fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f'Label file not found: {path}')
        df = pd.read_csv(path)
        return df.iloc[:, 0].values.astype(int)

    def load_survival(self, dataset: str, version: str = 'Original') -> pd.DataFrame:
        """
        Load survival data for a clustering dataset.

        Returns
        -------
        pd.DataFrame with columns: sample_name, event_observed, survival_times
        """
        cancer = self._canonical_cancer(dataset)
        folder = self._resolve_folder(dataset, version)
        path = os.path.join(folder, f'survival_{cancer}.csv')
        if not os.path.exists(path):
            raise FileNotFoundError(f'Survival file not found: {path}')
        return pd.read_csv(path)

    def get_X(
        self,
        dataset: str,
        version: str = 'Top',
        omics_types: Optional[List[str]] = None,
        concatenate: bool = True,
    ):
        """
        Load omics data as a sample-first matrix ready for ML models.

        Parameters
        ----------
        concatenate : bool
            If True, return a single np.ndarray of shape (n_samples, total_features)
            by concatenating all omics types.
            If False, return a dict {omics: np.ndarray(n_samples, n_features)}.

        Returns
        -------
        np.ndarray or dict of np.ndarray
        """
        omics_dict = self.load_omics(dataset, version, omics_types)
        # Each df is (features × samples); transpose → (samples × features)
        transposed = {k: df.T for k, df in omics_dict.items()}
        if not concatenate:
            return {k: df.values for k, df in transposed.items()}
        # Ensure sample order is consistent across omics before concatenating
        ref_index = list(transposed.values())[0].index
        arrays = []
        for omics, df in transposed.items():
            if not df.index.equals(ref_index):
                df = df.loc[ref_index]
            arrays.append(df.values)
        return np.concatenate(arrays, axis=1)

    def get_sample_names(self, dataset: str, version: str = 'Top') -> List[str]:
        """Return sample IDs for a dataset/version."""
        folder = self._resolve_folder(dataset, version)
        cancer = self._canonical_cancer(dataset)
        fname  = self._omics_filename(dataset, 'mRNA', version)
        path   = os.path.join(folder, fname)
        return pd.read_csv(path, index_col=0, nrows=0).columns.tolist()

    def get_feature_names(
        self, dataset: str, version: str = 'Top', omics: str = 'mRNA'
    ) -> List[str]:
        """Return feature (gene/probe) names for a given omics modality."""
        folder = self._resolve_folder(dataset, version)
        fname  = self._omics_filename(dataset, omics, version)
        path   = os.path.join(folder, fname)
        df = pd.read_csv(path, index_col=0, nrows=0)
        # index holds feature names after index_col=0
        return pd.read_csv(path, usecols=[0], header=0).iloc[:, 0].tolist()

    def dataset_info(self, dataset: str) -> dict:
        """Return a summary dict of available versions and shapes."""
        info = {'dataset': dataset, 'versions': {}}
        for version in ['Original', 'Top', 'Aligned']:
            try:
                folder = self._resolve_folder(dataset, version)
                shapes = {}
                for omics in OMICS_TYPES:
                    fname = self._omics_filename(dataset, omics, version)
                    path = os.path.join(folder, fname)
                    if os.path.exists(path):
                        df = pd.read_csv(path, index_col=0, nrows=0)
                        with open(path) as f:
                            n_feat = sum(1 for _ in f) - 1
                        shapes[omics] = (n_feat, len(df.columns))
                if shapes:
                    info['versions'][version] = shapes
            except (FileNotFoundError, KeyError):
                pass
        return info

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _canonical_cancer(self, dataset: str) -> str:
        """Strip 'GS-' prefix and 'Imp-' prefix to get the TCGA cancer code."""
        if dataset.startswith('GS-'):
            return dataset[3:]
        if dataset.startswith('Imp-'):
            return dataset[4:]
        return dataset

    def _resolve_folder(self, dataset: str, version: str) -> str:
        """Map (dataset, version) to the absolute folder path."""
        cancer = self._canonical_cancer(dataset)

        if dataset == 'Pan-cancer':
            return os.path.join(
                self.data_root, 'Classification_datasets', 'Pan-cancer', 'Original'
            )

        # Check Imp- prefix first (before GS classification lookup)
        if dataset.startswith('Imp-'):
            folder = os.path.join(
                self.data_root, 'Imputation_datasets', f'Imp-{cancer}', 'Top'
            )
            if os.path.isdir(folder):
                return folder

        if dataset.startswith('GS-'):
            folder = os.path.join(
                self.data_root, 'Classification_datasets', f'GS-{cancer}', version
            )
            if os.path.isdir(folder):
                return folder

        if cancer in CLUSTERING_DATASETS:
            folder = os.path.join(
                self.data_root, 'Clustering_datasets', cancer, version
            )
            if os.path.isdir(folder):
                return folder

        raise FileNotFoundError(
            f'Cannot resolve folder for dataset={dataset!r}, version={version!r}'
        )

    def _omics_filename(self, dataset: str, omics: str, version: str) -> str:
        """Return the expected filename for an omics file."""
        cancer = self._canonical_cancer(dataset)

        if dataset == 'Pan-cancer':
            return f'Pan-cancer_{omics}.csv'

        suffix_map = {'Original': '', 'Top': '_top', 'Aligned': '_aligned'}

        # Imputation datasets: Top version uses no suffix
        if dataset.startswith('Imp-'):
            return f'{cancer}_{omics}.csv'

        suffix = suffix_map.get(version, '')
        return f'{cancer}_{omics}{suffix}.csv'
