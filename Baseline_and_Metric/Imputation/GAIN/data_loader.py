# coding=utf-8
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''Data loader for UCI letter, spam and MNIST datasets.
'''

# Necessary packages
import numpy as np
from utils import binary_sampler
from keras.datasets import mnist


def data_loader (data_name, miss_rate, data_path=None):  # MLOmics: added data_path argument
  '''Loads datasets and introduce missingness.

  Args:
    - data_name: letter, spam, or mnist
    - miss_rate: the probability of missing components
    - data_path: (MLOmics) explicit path to the CSV file; if provided,
                 overrides the default relative path construction.
                 Use this to load MLOmics imputation datasets directly.

  Returns:
    data_x: original data
    miss_data_x: data with missing values
    data_m: indicator matrix for missing components
  '''

  # Load data
  # MLOmics: if data_path is supplied, use it directly; otherwise fall back to
  #          the original relative path (backward-compatible with letter/spam).
  if data_path is not None:  # MLOmics: explicit path
    file_name = data_path  # MLOmics: e.g. '/path/to/Main_Dataset/.../BRCA_mRNA.csv'
  else:
    file_name = '../../../Main_Dataset/Imputation_datasets/Imp-'+data_name.split("_")[0]+'/Top/'+data_name+'.csv'
  data_x = np.genfromtxt(file_name, delimiter=',', dtype=str)[1:,1:].astype(float)

  # Parameters
  no, dim = data_x.shape
  
  # Introduce missing data
  data_m = binary_sampler(1-miss_rate, no, dim)
  miss_data_x = data_x.copy()
  miss_data_x[data_m == 0] = np.nan
      
  return data_x, miss_data_x, data_m
