import numpy as np

from src.network.customics import CustOMICS
from src.tools.prepare_dataset import prepare_dataset

from torch.utils.data import DataLoader


def train(task, cohorts, sources, split, device, num_classes=4,
          batch_size=32, n_epochs=10, beta=1, lr=1e-3,
          hidden_dim=[512, 256], central_dim=[512, 256], latent_dim=128, dropout=0.2,
          classifier_dim=[128, 64], survival_dim=[64, 32], lambda_classif=5,
          lambda_survival=5, explain=False, explained_source='RNAseq', explained_class='Her2'):

    if cohorts == 'PANCAN':
        label = 'tumor_type'
        surv_time = 'OS.time'
        event = 'OS'
    else:
        label = 'pathology_T_stage'
        event = 'status'
        surv_time = 'overall_survival'

    omics_df, clinical_df, omics_train, omics_val, omics_test, lt_samples, x_dim = prepare_dataset(
        cohorts, sources, split, label)
    print(clinical_df[label].values)
    num_classes = len(np.unique(clinical_df[label].values))

    hidden_dim = [512, 256]
    central_dim = [512, 256]
    rep_dim = 128
    latent_dim = 128

    source_params = {}
    central_params = {'hidden_dim': central_dim, 'latent_dim': latent_dim,
                      'norm': True, 'dropout': dropout, 'beta': beta}
    classif_params = {'n_class': num_classes, 'lambda': lambda_classif,
                      'hidden_layers': classifier_dim, 'dropout': dropout}
    surv_params = {'lambda': lambda_survival, 'dims': survival_dim,
                   'activation': 'SELU', 'l2_reg': 1e-2, 'norm': True, 'dropout': dropout}
    for i, source in enumerate(sources):
        source_params[source] = {'input_dim': x_dim[i], 'hidden_dim': hidden_dim,
                                 'latent_dim': rep_dim, 'norm': True, 'dropout': 0.2}
    train_params = {'switch': 5, 'lr': 1e-3}

    model = CustOMICS(source_params=source_params, central_params=central_params, classif_params=classif_params,
                      surv_params=surv_params, train_params=train_params, device=device).to(device)
    model.get_number_parameters()
    model.fit(omics_train=omics_train, clinical_df=clinical_df, label=label, event=event, surv_time=surv_time,
              omics_val=omics_val, batch_size=batch_size, n_epochs=n_epochs, verbose=True)
    metric = model.evaluate(omics_test=omics_test, clinical_df=clinical_df, label=label, event=event, surv_time=surv_time,
                            task=task, batch_size=1024, plot_roc=False)
    model.plot_loss()
    model.plot_representation(omics_train, clinical_df, label,
                              'plot_representation', 'Representation of the latent space')
    if cohorts != 'PANCAN':
        model.explain(lt_samples, omics_df, clinical_df,
                      'RNAseq', 'Her2', label, device)
    if task == 'survival':
        print(model.predict_survival(omics_test))
        model.stratify(omics_df=omics_train, clinical_df=clinical_df,
                       event='status', surv_time='overall_survival')
    return metric
