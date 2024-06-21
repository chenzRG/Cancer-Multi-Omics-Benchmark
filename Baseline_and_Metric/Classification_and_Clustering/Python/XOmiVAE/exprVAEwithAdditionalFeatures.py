import torch
import generalHelperFunctions as helper
import numpy as np

from sklearn.model_selection import train_test_split
from torch import nn, optim
from torch.utils.data import Dataset, DataLoader
from torch.nn import functional as F
from torch.utils.tensorboard import SummaryWriter
from earlystoping import Earlystopping
from sklearn import metrics
import sys
import pandas as pd
sys.path.insert(0, '~/PycharmProjects/shapLundberg')


def ExprOmiVAE(input_path, expr_df, random_seed=42, no_cuda=False, model_parallelism=True,
               separate_testing=True, batch_size=32, latent_dim=128, learning_rate=1e-300, p1_epoch_num=50,
               p2_epoch_num=100, output_loss_record=True, classifier=True, early_stopping=True, Running=False,
               Pickle=True, kl_annealing=False, beta_vae=False,beta=1,knockingOut=False,dimension_number=16):
    """
    This altered version is primarily used for the disentangling analysis. Additionally, it includes the 'knock-out' test of altering the output of the model.
    :param Running: Set to true if we want to train VAE, otherwise if we just want the class structure, or other functions
    in this file, set to False.
    :param Pickle: Set to true if we have a saved file of the model.
    :param kl_annealing: Set to true if we want to attempt disentangling using KL annealing.
    :param beta_vae: Set to true if we want to attempt disentangling using a beta weight.
    :param beta: Specificy beta-weight here. If beta=1 then model acts as a normal VAE.
    :param knockingOut: Set to true if we want to reverse the output of the model
    :param dimension_number: Specify dimension output to reverse here.
    """

    torch.manual_seed(random_seed)
    torch.cuda.manual_seed_all(random_seed)
    device = torch.device('cuda:0' if not no_cuda and torch.cuda.is_available() else 'cpu')
    parallel = torch.cuda.device_count() > 1 and model_parallelism

    class_num, label_array, sample_id = helper.loadData(input_path)

    # cleans expr_df if pickle not provided
    if not Pickle:
        expr_df=helper.preprocessExpr_df(expr_df)
        # save cleaned data frame to a pickle to be used later
        expr_df.to_pickle("data/expr_df.pkl")

    if separate_testing:
        # Get testing set index and training set index
        # Separate according to different tumour types
        testset_ratio = 0.2
        valset_ratio = 0.5
        train_index, test_index, train_label, test_label = train_test_split(sample_id, label_array,
                                                                            test_size=testset_ratio,
                                                                            random_state=random_seed,
                                                                            stratify=label_array)
        (uniquetwo, countstwo) = np.unique(test_label, return_counts=True)
        frequenciesq = np.asarray((uniquetwo, countstwo)).T


        val_index, test_index, val_label, test_label = train_test_split(test_index, test_label, test_size=valset_ratio,
                                                                        random_state=random_seed)
        expr_df_test = expr_df[test_index]
        expr_df_val = expr_df[val_index]
        expr_df_train = expr_df[train_index]


    # Get multi-omics dataset information
    sample_num = len(sample_id)
    expr_feature_num = expr_df.shape[0]
    """
    print('\nNumber of samples: {}'.format(sample_num))
    print('Number of gene expression features: {}'.format(expr_feature_num))
    """
    if classifier:
        print('Number of classes: {}'.format(class_num))

    class ExprOmiDataset(Dataset):
        def __init__(self, expr_df, labels):
            self.expr_df = expr_df
            self.labels = labels

        def __len__(self):
            return self.expr_df.shape[1]

        def __getitem__(self, index):
            expr_line = self.expr_df.iloc[:, index].values

            expr_line = expr_line.astype(dtype='float32')
            expr_line_tensor = torch.Tensor(expr_line)
            label = self.labels[index]
            return [expr_line_tensor, label]

    # DataSets and DataLoaders
    if separate_testing:
        train_dataset = ExprOmiDataset(expr_df=expr_df_train, labels=train_label)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=6)
        val_dataset = ExprOmiDataset(expr_df=expr_df_val, labels=val_label)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=True, num_workers=6)
        test_dataset = ExprOmiDataset(expr_df=expr_df_test, labels=test_label)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True, num_workers=6)
    else:
        train_dataset = ExprOmiDataset(expr_df=expr_df, labels=label_array)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=6)

    full_dataset = ExprOmiDataset(expr_df=expr_df, labels=label_array)
    full_loader = DataLoader(full_dataset, batch_size=batch_size, num_workers=6)

    # Setting dimensions
    latent_space_dim = latent_dim
    input_dim_expr = expr_feature_num
    level_2_dim_expr = 4096
    level_3_dim_expr = 1024
    level_4_dim = 512
    classifier_1_dim = 128
    classifier_2_dim = 64
    classifier_out_dim = class_num

    class VAE(nn.Module):
        def __init__(self):
            super(VAE, self).__init__()
            # ENCODER fc layers
            # level 1
            # Expr input
            self.e_fc1_expr = self.fc_layer(input_dim_expr, level_2_dim_expr)

            # Level 2
            self.e_fc2_expr = self.fc_layer(level_2_dim_expr, level_3_dim_expr)
            # self.e_fc2_expr = self.fc_layer(level_2_dim_expr, level_3_dim_expr, dropout=True)

            # Level 3
            self.e_fc3 = self.fc_layer(level_3_dim_expr, level_4_dim)
            # self.e_fc3 = self.fc_layer(level_3_dim_expr, level_4_dim, dropout=True)

            # Level 4
            self.e_fc4_mean = self.fc_layer(level_4_dim, latent_space_dim, activation=0)
            self.e_fc4_log_var = self.fc_layer(level_4_dim, latent_space_dim, activation=0)
            self.device = torch.device('cuda:0' if not no_cuda and torch.cuda.is_available() else 'cpu')

            # model parallelism
            if parallel:
                self.e_fc1_expr.to('cuda:0')
                self.e_fc2_expr.to('cuda:0')
                self.e_fc3.to('cuda:0')
                self.e_fc4_mean.to('cuda:0')
                self.e_fc4_log_var.to('cuda:0')

            # DECODER fc layers
            # Level 4
            self.d_fc4 = self.fc_layer(latent_space_dim, level_4_dim)

            # Level 3
            self.d_fc3 = self.fc_layer(level_4_dim, level_3_dim_expr)
            # self.d_fc3 = self.fc_layer(level_4_dim, level_3_dim_expr, dropout=True)

            # Level 2
            self.d_fc2_expr = self.fc_layer(level_3_dim_expr, level_2_dim_expr)
            # self.d_fc2_expr = self.fc_layer(level_3_dim_expr, level_2_dim_expr, dropout=True)

            # level 1
            # Expr output
            self.d_fc1_expr = self.fc_layer(level_2_dim_expr, input_dim_expr, activation=2)

            # model parallelism
            if parallel:
                self.d_fc4.to('cuda:1')
                self.d_fc3.to('cuda:1')
                self.d_fc2_expr.to('cuda:1')
                self.d_fc1_expr.to('cuda:1')

            # CLASSIFIER fc layers
            self.c_fc1 = self.fc_layer(latent_space_dim, classifier_1_dim)
            self.c_fc2 = self.fc_layer(classifier_1_dim, classifier_2_dim)
            # self.c_fc2 = self.fc_layer(classifier_1_dim, classifier_2_dim, dropout=True)
            self.c_fc3 = self.fc_layer(classifier_2_dim, classifier_out_dim, activation=0)

            # model parallelism
            if parallel:
                self.c_fc1.to('cuda:1')
                self.c_fc2.to('cuda:1')
                self.c_fc3.to('cuda:1')

        # Activation - 0: no activation, 1: ReLU, 2: Sigmoid
        def fc_layer(self, in_dim, out_dim, activation=1, dropout=False, dropout_p=0.5):
            if activation == 0:
                layer = nn.Sequential(
                    nn.Linear(in_dim, out_dim),
                    nn.BatchNorm1d(out_dim))
            elif activation == 2:
                layer = nn.Sequential(
                    nn.Linear(in_dim, out_dim),
                    nn.BatchNorm1d(out_dim),
                    nn.Sigmoid())
            else:
                if dropout:
                    layer = nn.Sequential(
                        nn.Linear(in_dim, out_dim),
                        nn.BatchNorm1d(out_dim),
                        nn.ReLU(),
                        nn.Dropout(p=dropout_p))
                else:
                    layer = nn.Sequential(
                        nn.Linear(in_dim, out_dim),
                        nn.BatchNorm1d(out_dim),
                        nn.ReLU())

            return layer

        def encode(self, x):

            expr_level2_layer = self.e_fc1_expr(x)

            level_3_layer = self.e_fc2_expr(expr_level2_layer)

            level_4_layer = self.e_fc3(level_3_layer)

            latent_mean = self.e_fc4_mean(level_4_layer)
            latent_log_var = self.e_fc4_log_var(level_4_layer)

            return latent_mean, latent_log_var

        def reparameterize(self, mean, log_var):

            sigma = torch.exp(0.5 * log_var)
            eps = torch.randn_like(sigma)
            return mean + eps * sigma

        def decode(self, z):

            level_4_layer = self.d_fc4(z)

            level_3_layer = self.d_fc3(level_4_layer)

            expr_level2_layer = self.d_fc2_expr(level_3_layer)

            recon_x = self.d_fc1_expr(expr_level2_layer)

            return recon_x

        def classifier(self, mean):
            #Output change (4th in figure) variant of OmiVAE knock-out experiments
            #size of input x mean
            if knockingOut:
                mean=mean.numpy()
                dimension_number=85
                mean_overZero = mean[:, dimension_number] > 0
                mean_Zero = mean[:, dimension_number] <= 0
                mean[:, dimension_number][mean_overZero] = 0.
                mean[:, dimension_number][mean_Zero] = 1.
                mean=torch.from_numpy(mean)
            level_1_layer = self.c_fc1(mean)
            level_2_layer = self.c_fc2(level_1_layer)
            output_layer = self.c_fc3(level_2_layer)
            return output_layer

        def forward(self, x):

            mean, log_var = self.encode(x)

            z = self.reparameterize(mean, log_var)

            classifier_x = mean
            if parallel:
                z = z.to('cuda:1')

                classifier_x = classifier_x.to('cuda:1')

            recon_x = self.decode(z)
            pred_y = self.classifier(classifier_x)
            z = z.to(device)
            recon_x = recon_x.to(device)
            mean = mean.to(device)
            log_var = log_var.to(device)
            pred_y = pred_y.to(device)
            """
            with torch.no_grad():
                pred_y_softmax = F.softmax(pred_y, dim=1)
                _, predicted = torch.max(pred_y_softmax, 1)
                print("predicted is")
                print(predicted)
            """



            return [z, recon_x, mean, log_var, pred_y]

    # Instantiate VAE
    if parallel:
        vae_model = VAE()
    else:
        vae_model = VAE().to(device)

    # Early Stopping
    if early_stopping:
        early_stop_ob = Earlystopping()

    # Tensorboard writer
    train_writer = SummaryWriter(log_dir='../logs/train')
    val_writer = SummaryWriter(log_dir='../logs/val')

    # print the model information
    # print('\nModel information:')
    # print(vae_model)
    total_params = sum(params.numel() for params in vae_model.parameters())
    print('Number of parameters: {}'.format(total_params))

    optimizer = optim.Adam(vae_model.parameters(), lr=learning_rate)

    def expr_recon_loss(recon_x, x):
        loss = F.binary_cross_entropy(recon_x, x, reduction='sum')
        return loss

    def kl_loss(mean, log_var):
        loss = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())
        return loss

    def classifier_loss(pred_y, y):
        loss = F.cross_entropy(pred_y, y, reduction='sum')
        return loss

    # loss record
    loss_array = np.zeros(shape=(9, p1_epoch_num + p2_epoch_num + 1))
    # performance metrics
    metrics_array = np.zeros(4)

    def train(e_index, e_num, k_expr_recon, k_kl, k_c,K_weight=1):
        vae_model.train()
        train_expr_recon = 0
        train_kl = 0
        train_classifier = 0
        train_correct_num = 0
        train_total_loss = 0
        print("enumerate loader is")
        print(enumerate(train_loader))
        for batch_index, sample in enumerate(train_loader):
            data = sample[0]
            y = sample[1]
            data = data.to(device)
            y.to(device)

            y = y.to(device)
            optimizer.zero_grad()
            _, recon_data, mean, log_var, pred_y = vae_model(data)
            if parallel:
                recon_data = recon_data.to('cuda:0')
                pred_y = pred_y.to('cuda:0')
            expr_recon = expr_recon_loss(recon_data, data)
            kl = kl_loss(mean, log_var)
            class_loss = classifier_loss(pred_y, y)
            loss = k_expr_recon * expr_recon + k_kl * kl * K_weight + k_c * class_loss

            loss.backward()

            with torch.no_grad():
                pred_y_softmax = F.softmax(pred_y, dim=1)
                _, predicted = torch.max(pred_y_softmax, 1)
                correct = (predicted == y).sum().item()
                train_expr_recon += expr_recon.item()
                train_kl += kl.item()
                train_classifier += class_loss.item()
                train_correct_num += correct
                train_total_loss += loss.item()

            optimizer.step()

            # if batch_index % 15 == 0:
            #     print('Epoch {:3d}/{:3d}  ---  [{:5d}/{:5d}] ({:2d}%)\n'
            #           'Expr Recon Loss: {:.2f}   KL Loss: {:.2f}   '
            #           'Classification Loss: {:.2f}\nACC: {:.2f}%'.format(
            #         e_index + 1, e_num, batch_index * len(data), len(train_dataset),
            #         round(100. * batch_index / len(train_loader)),
            #         expr_recon.item() / len(data), kl.item() / len(data), class_loss.item() / len(data),
            #         correct / len(data) * 100))

        train_expr_recon_ave = train_expr_recon / len(train_dataset)
        train_kl_ave = train_kl / len(train_dataset)
        train_classifier_ave = train_classifier / len(train_dataset)
        train_accuracy = train_correct_num / len(train_dataset) * 100
        train_total_loss_ave = train_total_loss / len(train_dataset)

        print('Epoch {:3d}/{:3d}\n'
              'Training\n'
              'Expr Recon Loss: {:.2f}   KL Loss: {:.2f}   '
              'Classification Loss: {:.2f}\nACC: {:.2f}%'.
              format(e_index + 1, e_num, train_expr_recon_ave, train_kl_ave, train_classifier_ave, train_accuracy))
        loss_array[0, e_index] = train_expr_recon_ave
        loss_array[1, e_index] = train_kl_ave
        loss_array[2, e_index] = train_classifier_ave
        loss_array[3, e_index] = train_accuracy

        # TB
        train_writer.add_scalar('Total loss', train_total_loss_ave, e_index)
        train_writer.add_scalar('Expr recon loss', train_expr_recon_ave, e_index)
        train_writer.add_scalar('KL loss', train_kl_ave, e_index)
        train_writer.add_scalar('Classification loss', train_classifier_ave, e_index)
        train_writer.add_scalar('Accuracy', train_accuracy, e_index)

    if separate_testing:
        def val(e_index, get_metrics=False):
            vae_model.eval
            val_expr_recon = 0
            val_kl = 0
            val_classifier = 0
            val_correct_num = 0
            val_total_loss = 0
            y_store = torch.tensor([0])
            predicted_store = torch.tensor([0])

            with torch.no_grad():
                for batch_index, sample in enumerate(val_loader):
                    data = sample[0]
                    y = sample[1]
                    data = data.to(device)
                    y = y.to(device)
                    _, recon_data, mean, log_var, pred_y = vae_model(data)
                    if parallel:
                        recon_data = recon_data.to('cuda:0')
                        pred_y = pred_y.to('cuda:0')

                    expr_recon = expr_recon_loss(recon_data, data)
                    kl = kl_loss(mean, log_var)
                    class_loss = classifier_loss(pred_y, y)
                    loss = expr_recon + kl + class_loss

                    pred_y_softmax = F.softmax(pred_y, dim=1)
                    _, predicted = torch.max(pred_y_softmax, 1)
                    correct = (predicted == y).sum().item()

                    y_store = torch.cat((y_store, y.cpu()))
                    predicted_store = torch.cat((predicted_store, predicted.cpu()))

                    val_expr_recon += expr_recon.item()
                    val_kl += kl.item()
                    val_classifier += class_loss.item()
                    val_correct_num += correct
                    val_total_loss += loss.item()

            output_y = y_store[1:].numpy()
            output_pred_y = predicted_store[1:].numpy()

            if get_metrics:
                metrics_array[0] = metrics.accuracy_score(output_y, output_pred_y)
                metrics_array[1] = metrics.precision_score(output_y, output_pred_y, average='weighted')
                metrics_array[2] = metrics.recall_score(output_y, output_pred_y, average='weighted')
                metrics_array[3] = metrics.f1_score(output_y, output_pred_y, average='weighted')

            val_expr_recon_ave = val_expr_recon / len(val_dataset)
            val_kl_ave = val_kl / len(val_dataset)
            val_classifier_ave = val_classifier / len(val_dataset)
            val_accuracy = val_correct_num / len(val_dataset) * 100
            val_total_loss_ave = val_total_loss / len(val_dataset)

            print('Validation\n'
                  'Expr Recon Loss: {:.2f}   KL Loss: {:.2f}   Classification Loss: {:.2f}'
                  '\nACC: {:.2f}%\n'.
                  format(val_expr_recon_ave, val_kl_ave, val_classifier_ave, val_accuracy))
            loss_array[4, e_index] = val_expr_recon_ave
            loss_array[5, e_index] = val_kl_ave
            loss_array[6, e_index] = val_classifier_ave
            loss_array[7, e_index] = val_accuracy

            # TB
            val_writer.add_scalar('Total loss', val_total_loss_ave, e_index)
            val_writer.add_scalar('Expr recon loss', val_expr_recon_ave, e_index)
            val_writer.add_scalar('KL loss', val_kl_ave, e_index)
            val_writer.add_scalar('Classification loss', val_classifier_ave, e_index)
            val_writer.add_scalar('Accuracy', val_accuracy, e_index)

            return val_accuracy, output_pred_y

    # This is skipped if we already have a OmiVAE pickle model
    if Running:
        print('\nUNSUPERVISED PHASE\n')
        # unsupervised phase
        klstart = 3
        # weight of KL loss; 1 is normal vae ratio... this starts at 0 and increase s
        if kl_annealing:
            weight = 0
        if beta_vae:
            # if want beta vae then change weight here!
            weight=beta
        else:
            weight=1

        for epoch_index in range(p1_epoch_num):
            if kl_annealing:
                if epoch_index > klstart:
                    weight = 5
                if epoch_index > (klstart + 2):
                    weight = 0.7
            train(e_index=epoch_index, e_num=p1_epoch_num + p2_epoch_num, k_expr_recon=1, k_kl=1, k_c=0,K_weight=weight)
            if separate_testing:
                _, out_pred_y = val(epoch_index)

        print("finishing unsupervised phase")
        # find z (latent space) values after unsupervised phase
        vae_model.eval()
        expr_tensor = expr_df.values
        expr_tensor = expr_tensor.astype(dtype='float32')
        expr_tensor = torch.Tensor(expr_tensor)
        expr_tensor = torch.transpose(expr_tensor, 0, 1)
        with torch.no_grad():
            _, recon_data, mean, log_var, pred_y = vae_model(expr_tensor)
            np.savetxt('/content/drive/My Drive/data/mean_before_supervised.csv', mean)

        # Do not need to run this phase if performing disentangling experiments
        print('\nSUPERVISED PHASE\n')
        # supervised phase
        epoch_number = p1_epoch_num
        for epoch_index in range(p1_epoch_num, p1_epoch_num + p2_epoch_num):
            epoch_number += 1
            #Run like normal
            train(e_index=epoch_index, e_num=p1_epoch_num + p2_epoch_num, k_expr_recon=0, k_kl=0, k_c=1)
            if separate_testing:
                if epoch_index == p1_epoch_num + p2_epoch_num - 1:
                    val_classification_acc, out_pred_y = val(epoch_index, get_metrics=True)
                else:
                    val_classification_acc, out_pred_y = val(epoch_index)
                if early_stopping:
                    early_stop_ob(vae_model, val_classification_acc)
                    if early_stop_ob.stop_now:
                        print('Early stopping\n')
                        break

        if early_stopping:
            best_epoch = p1_epoch_num + early_stop_ob.best_epoch_num
            loss_array[8, 0] = best_epoch
            print('Load model of Epoch {:d}'.format(best_epoch))
            # Get ride of map location for colab
            vae_model.load_state_dict(
                torch.load('/content/drive/My Drive/data/checkpoint.pt', map_location=torch.device('cpu')))
            _, out_pred_y = val(epoch_number, get_metrics=True)
        print("finishing SUPERVISED phase")
        # find z values after unsupervised phase
        vae_model.eval()
        expr_tensor = expr_df.values
        expr_tensor = expr_tensor.astype(dtype='float32')
        expr_tensor = torch.Tensor(expr_tensor)
        expr_tensor = torch.transpose(expr_tensor, 0, 1)
        with torch.no_grad():
            _, recon_data, mean, log_var, pred_y = vae_model(expr_tensor)
            np.savetxt('/content/drive/My Drive/data/mean_after_supervised.csv', mean)

        print("y predictions are")
        print(out_pred_y)

        # Encode all of the data into the latent space
        print('Encoding all the data into latent space...')
        vae_model.eval()
        with torch.no_grad():
            d_z_store = torch.zeros(1, latent_dim).to(device)
            for i, sample in enumerate(full_loader):
                d = sample[0]
                d = d.to(device)
                _, _, d_z, _, _ = vae_model(d)
                d_z_store = torch.cat((d_z_store, d_z), 0)
        all_data_z = d_z_store[1:]
        all_data_z_np = all_data_z.cpu().numpy()

        latent_space_path = '/content/drive/My Drive/data/D_latent_space.tsv'

        all_data_z_df = pd.DataFrame(all_data_z_np, index=sample_id)
        all_data_z_df.to_csv(latent_space_path, sep='\t')

        if separate_testing:
            pred_y_path = '/content/drive/My Drive/data/' + '128_D_pred_y.tsv'
            np.savetxt(pred_y_path, out_pred_y, delimiter='\t')

            metrics_record_path = '/content/drive/My Drive/data/' + '128_D_metrics.tsv'
            np.savetxt(metrics_record_path, metrics_array, delimiter='\t')

        if output_loss_record:
            loss_record_path = '/content/drive/My Drive/data/' + 'D_loss_record.tsv'
            np.savetxt(loss_record_path, loss_array, delimiter='\t')


    return vae_model




