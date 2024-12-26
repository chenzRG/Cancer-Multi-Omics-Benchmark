
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#Bar chart code for disentanglement experiments
#Reference: @DOI={10.1101/801605}
def compute_diff_capacity_latent_dim(latent_dim_data, subtype_z, labels, model):
    # Compute the percentage distribution for the tissue more than a standard deviation away from the mean
    # change depending on number of columns
    subtype_latent = np.zeros(shape=(subtype_z.shape[1], 4))
    print(subtype_z.shape[1])
    for latent_dim in range(latent_dim_data.shape[1]):
        latent_dim_across_subtype = subtype_z[:, latent_dim]
        latent_dim_across_cells = latent_dim_data[:, latent_dim]
        # 1x1
        latent_dim_mean = np.mean(latent_dim_across_cells)
        latent_dim_std = np.std(latent_dim_across_cells)
        # Return the indices of the elements that are non-zero
        variable_cells = np.where(latent_dim_across_subtype > latent_dim_mean + latent_dim_std)
        variable_cells_left = np.where(latent_dim_across_subtype < latent_dim_mean - latent_dim_std)
        variable_labels = labels[variable_cells]
        variable_two_labels = labels[variable_cells_left]
        variable_cells = variable_labels.tolist()
        variable_cells_left = variable_two_labels.tolist()
        # Add one to end of range.. these are the label codes
        counter_dict = {x: float(variable_cells.count(x)) for x in range(11, 14)}
        counter_dict_two = {x: float(variable_cells_left.count(x)) for x in range(11, 14)}

        counter = (np.array(list(counter_dict.values())) - np.array(list(counter_dict_two.values()))) * -1
        # counter=np.abs(counter)
        # counter =  np.array(list(counter_dict.values()))/float(len(latent_dim_across_cells))
        # counter = np.around(counter * 100.0, decimals=2)
        subtype_latent[latent_dim][1:] = counter
        subtype_latent[latent_dim][0] = int(latent_dim)

    subtype_latent = pd.DataFrame(subtype_latent, columns=['Dim.', 'KIRC', 'KIRP', 'KICH'])
    subtype_latent['Dim.'] = subtype_latent['Dim.'].astype(int)

    subtype_latent = subtype_latent.melt(id_vars=['Dim.'], value_vars=['KIRC', 'KIRP', 'KICH'],
                                         var_name='Tissue type', value_name='Percentage')
    print(subtype_latent)
    sns.set(font_scale=2.5)
    # flatui = ["#9b59b6", "#2ecc71", "#95a5a6", "#e74c3c", "#3498db" ]
    flatui = ["#9b59b6", "#2ecc71", "#95a5a6"]
    sns.set_palette(sns.color_palette(flatui))
    sns.set_style('darkgrid')
    g = sns.factorplot(x='Tissue type', y='Percentage', col='Dim.', data=subtype_latent, saturation=.5,
                       col_wrap=5,
                       kind="bar", ci=None, aspect=1.3, legend_out=True)

    g.set_xticklabels(rotation=70)

    plt.show()
    g.savefig("data/CellDifferentiation" + model + ".pdf")


def compute_diff_capacity_latent_dimBRCA(latent_dim_data, subtype_z, labels, model):
    # Compute the percentage distribution for the tissue more than a standard deviation away from the mean
    # change depending on number of columns e.g. dim, brca, normal then=3
    subtype_latent = np.zeros(shape=(subtype_z.shape[1], 35))

    print(subtype_z.shape[1])
    for latent_dim in range(latent_dim_data.shape[1]):
        latent_dim_across_subtype = subtype_z[:, latent_dim]
        latent_dim_across_cells = latent_dim_data[:, latent_dim]
        # 1x1
        latent_dim_mean = np.mean(latent_dim_across_cells)
        latent_dim_std = np.std(latent_dim_across_cells)
        # Return the indices of the elements that are non-zero... note this does +/- I believe
        variable_cells = np.where(latent_dim_across_subtype > latent_dim_mean + latent_dim_std)
        variable_cells_left = np.where(latent_dim_across_subtype < latent_dim_mean - latent_dim_std)
        variable_labels = labels[variable_cells]
        variable_two_labels = labels[variable_cells_left]
        variable_cells = variable_labels.tolist()
        variable_cells_left = variable_two_labels.tolist()
        # I think add one to end of range
        counter_dict = {x: float(variable_cells.count(x)) for x in range(0, 34)}
        print("counter to left")
        print(counter_dict)
        counter_dict_two = {x: float(variable_cells_left.count(x)) for x in range(0, 34)}
        print("counter to right")
        print(counter_dict_two)
        # more interesting if length is all cells.. do that i think ere
        print("length")
        print(len(latent_dim_across_cells))
        # trying just counterof values?
        counter = (np.array(list(counter_dict.values())) - np.array(list(counter_dict_two.values())))
        # counter=np.abs(counter)
        # counter =  np.array(list(counter_dict.values()))/float(len(latent_dim_across_cells))
        # counter = np.around(counter * 100.0, decimals=2)
        subtype_latent[latent_dim][1:] = counter
        subtype_latent[latent_dim][0] = int(latent_dim)

    subtype_latent = pd.DataFrame(subtype_latent, columns=['Dim.', 'BRCA', 'normal'])
    subtype_latent['Dim.'] = subtype_latent['Dim.'].astype(int)

    subtype_latent = subtype_latent.melt(id_vars=['Dim.'], value_vars=['BRCA', 'normal'],
                                         var_name='Tissue type', value_name='Count')
    print(subtype_latent)
    sns.set(font_scale=2.5)
    # flatui = ["#9b59b6", "#2ecc71", "#95a5a6", "#e74c3c", "#3498db" ]
    flatui = ["#3498db", "#e74c3c", "#95a5a6"]
    sns.set_palette(sns.color_palette(flatui))
    sns.set_style('darkgrid')
    g = sns.factorplot(x='Tissue type', y='Count', col='Dim.', data=subtype_latent, saturation=.5,
                       col_wrap=5, kind="bar", ci=None, aspect=1.3, legend_out=True)

    g.set_xticklabels(rotation=70)

    plt.show()
    g.savefig("data/" + model + ".pdf")