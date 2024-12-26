import numpy as np
import pickle
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from sklearn.manifold import TSNE
import matplotlib.patheffects as PathEffects


sns.set_style('darkgrid')
sns.set_palette('muted')
sns.set_context("notebook", font_scale=1.5,
                rc={"lines.linewidth": 2.5})

def get_common_samples(dfs):
    lt_indices = []
    for df in dfs:
        lt_indices.append(list(df.index))
    common_indices = set(lt_indices[0])
    for i in range(1, len(lt_indices)):
        common_indices = common_indices & set(lt_indices[i])
    return list(common_indices)

def get_sub_omics_df(omics_df, lt_samples):
    return {key: value.loc[lt_samples, :] for key, value in omics_df.items()}


def read_data(cohort, omic_sources, label):
    omics_df = {}
    cnv_path = '../TCGA/{}/Omics/CNV/{}.gistic.tsv'.format(cohort,cohort)
    rna_path = '../TCGA/{}/Omics/RNAseq/{}.htseq_fpkm.tsv'.format(cohort,cohort)
    methyl_path = '../TCGA/{}/Omics/Methyl450/{}.methylation450.tsv'.format(cohort,cohort)
    clinical_path = '../TCGA/{}/Clinical/{}.GDC_phenotype.tsv'.format(cohort,cohort)
    if 'CNV' in omic_sources:
        omics_df['CNV'] = pd.read_csv(cnv_path, sep='\t', index_col=0).T.dropna(axis=1)
    if 'RNAseq' in omic_sources:
        omics_df['RNAseq'] = pd.read_csv(rna_path, sep='\t', index_col=0).T.dropna(axis=1)
    if 'methyl' in omic_sources:
        omics_df['methyl'] = pd.read_csv(methyl_path, sep='\t', index_col=0).T.dropna(axis=1)
    clinical_df = pd.read_csv(clinical_path, sep='\t', index_col=0).T
    clinical_df = clinical_df[clinical_df['overall_survival'].notna()]
    if cohort == 'TCGA-BRCA':
        clinical_df = clinical_df[clinical_df[label].notna()]
    lt_samples = get_common_samples([df for df in omics_df.values()] + [clinical_df])
    for source in omics_df.keys():
        omics_df[source] = omics_df[source].loc[lt_samples, :]
    return omics_df, clinical_df.loc[lt_samples,:], lt_samples

    


def omics_train_test_split(omics_df, data_y, split_test, seed=2):
    n_samples = len(data_y)
    lt_samples = list(range(n_samples))
    i_train, i_test = train_test_split(lt_samples, test_size=split_test, random_state=seed)
    y_train, y_test = data_y[i_train], data_y[i_test]
    omics_train = {}
    omics_test = {}
    for name in omics_df.keys():
        omics_train[name] = omics_df[name].iloc[i_train, :]
        omics_test[name] = omics_df[name].iloc[i_test, :]
    return omics_train, omics_test, y_train, y_test

def to_categorical(y_train, y_val, y_test):
    le = LabelEncoder()
    le.fit(y_train)
    y_train_le = le.transform(y_train)
    y_val_le = le.transform(y_val)
    y_test_le = le.transform(y_test)
    return y_train_le, y_val_le ,y_test_le, le


def from_categorical(y_train, y_test, le):
    return le.inverse_transform(y_train), le.inverse_transform(y_test)



#####PANCAN

def get_samples():
    lt = []
    with open('src/config/samples.txt', 'r') as sample_file:
        for line in sample_file.readlines():
            lt.append(line.rstrip('\n'))
    return lt

def read_data_pancan(omic_sources):
    omics_df = {}
    for omic_source in omic_sources:
        directory = '../TCGA/Pancancer/processed_pancan/{}'.format(omic_source)
        chr = 0
        for filename in os.listdir(directory):
            if omic_source == 'Methyl':
                with open(directory + '/' + filename, 'rb') as omic_file:
                    df = pickle.load(omic_file)
                    chr += 1
                omics_df[omic_source + str(chr)] = df
            else:
                with open(directory + '/' + filename, 'rb') as omic_file:
                    df = pickle.load(omic_file)
                if omic_source == 'RNAseq':
                    omics_df[omic_source] = df.T
                else:
                    omics_df[omic_source] = df
    with open('../TCGA/Pancancer/processed_pancan/clinical/clinical.pkl', 'rb') as clinical_file:
        df_clinical = pickle.load(clinical_file)
    lt_samples = get_common_samples([df for df in omics_df.values()] + [df_clinical])
    for source in omics_df.keys():
        omics_df[source] = omics_df[source].loc[lt_samples, :]
    return omics_df, df_clinical.loc[lt_samples, :], df_clinical.loc[lt_samples, :].values, lt_samples

def extract_tumour_type(data):
    def get_tumour_name_from_id(project_id):
        pos = 0
        for i in range(len(project_id)):
            if project_id[i] == '-':
                pos = i
        return project_id[pos+1:]
    data_tumour = []
    for sample in data:
        if sample [3] == 'Solid Tissue Normal':
            data_tumour.append(['Normal', int(sample[-2]), float(sample[-1])])
        else:
            data_tumour.append([get_tumour_name_from_id(sample[-6]), float(sample[4]) ,int(sample[-1]), float(sample[-3])])
    return np.array(data_tumour)


#######


def save_splits(lt_samples, cohort):
    kf = KFold(n_splits=5)
    for i, (train_index, test_index) in enumerate(kf.split(lt_samples)):
        train_index, val_index = train_test_split(train_index, test_size=0.15)
        samples_train = [lt_samples[i] for i in train_index]
        samples_val = [lt_samples[i] for i in val_index]
        samples_test = [lt_samples[i] for i in test_index]
        split_dir = 'splits/{}/'.format(cohort)
        if not os.path.isdir(split_dir):
            os.mkdir(split_dir)
        with open(split_dir + 'split_train_{}.txt'.format(i+1), 'w+') as split_file:
            for sample in samples_train:
                split_file.write(sample + '\n')
        with open(split_dir + 'split_val_{}.txt'.format(i+1), 'w+') as split_file:
            for sample in samples_val:
                split_file.write(sample + '\n')
        with open(split_dir + 'split_test_{}.txt'.format(i+1), 'w+') as split_file:
            for sample in samples_test:
                split_file.write(sample + '\n')

def get_splits(cohort, split):
    split_dir = 'splits/{}/'.format(cohort)
    with open(split_dir + 'split_train_{}.txt'.format(split), 'r') as split_file:
        samples_train = [sample.rstrip() for sample in split_file.readlines()]
    with open(split_dir + 'split_val_{}.txt'.format(split), 'r') as split_file:
        samples_val = [sample.rstrip() for sample in split_file.readlines()]
    with open(split_dir + 'split_test_{}.txt'.format(split), 'r') as split_file:
        samples_test = [sample.rstrip() for sample in split_file.readlines()]
    return samples_train, samples_val, samples_test


def fashion_scatter(x, colors):
    # choose a color palette with seaborn.
    num_classes = len(np.unique(colors))
    palette = np.array(sns.color_palette("hls", num_classes))

    # create a scatter plot.
    f = plt.figure(figsize=(8, 8))
    ax = plt.subplot(aspect='equal')
    sc = ax.scatter(x[:,0], x[:,1], lw=0, s=40, c=palette[colors.astype(np.int)])
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)
    ax.axis('off')
    ax.axis('tight')

    # add the labels for each digit corresponding to the label
    txts = []

    for i in range(num_classes):

        # Position of each label at median of data points.

        xtext, ytext = np.median(x[colors == i, :], axis=0)
        txt = ax.text(xtext, ytext, str(i), fontsize=24)
        txt.set_path_effects([
            PathEffects.Stroke(linewidth=5, foreground="w"),
            PathEffects.Normal()])
        txts.append(txt)

    return f, ax, sc, txts

def save_plot_score(filename, z, y, title, show=False):
    tsne = TSNE(n_components=2, verbose=0, perplexity=40, n_iter=300)
    tsne_embed = tsne.fit_transform(z)
    df = pd.DataFrame()
    df['targets'] = y
    df['x-axis'] = tsne_embed[:,0]
    df['y-axis'] = tsne_embed[:,1]
    #fashion_scatter(tsne_embed, y)
    sns_plot = sns.scatterplot(x='x-axis', y='y-axis', hue=df.targets.tolist(),
                    palette=sns.color_palette('hls', len(np.unique(y))),data=df).set(title=title)
    plt.legend(bbox_to_anchor=(1.5, 1.1), loc=2, borderaxespad=0.)
    plt.savefig(filename +  '.png', bbox_inches='tight')
    if show:
        plt.show()
    # Put the legend out of the figure
    plt.clf()


