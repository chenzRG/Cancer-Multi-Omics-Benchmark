import argparse
import torch

import numpy as np

from src.tools.core_utils import train


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cohorts', help='list of cohorts to process', type=str)
parser.add_argument('-dv', '--device', help='torch device in which the computations will be done', type=str, default='cpu')
parser.add_argument('-dr', '--data_directory', help='folder in which the data are stored', type=str, default='../TCGA/')
parser.add_argument('-res', '--result_directory', help='folder in which the results should be stored', type=str, default='results/')
parser.add_argument('-t', '--task', help='task to perform', type=str, choices=['classification', 'survival'], default='classification')
parser.add_argument('-src', '--sources', help='list of sources to integrate', type=str, default='CNV,RNAseq,methyl')


parser.add_argument('-nc', '--num_classes', help='number of classes for the classification task', type=int, default=4)
parser.add_argument('-b', '--batch_size', help='batch size for the data loader', type=int, default=32)
parser.add_argument('-e', '--epochs', help='number of training epochs', type=int, default=20)
parser.add_argument('-p2', '--p2_switch', help='epoch to switch to phase 2', type=int, default=10)
parser.add_argument('-lr', '--lr', help='learning rate for the training optimizer', type=int, default=1e-3)
parser.add_argument('-bt', '--beta', help='value of the regulation coefficient for the beta-VAE', type=int, default=1)
parser.add_argument('-dp', '--dropout', help='dropout rate', type=float, default=0.2)


parser.add_argument('-hd', '--hidden_dim', help='list of neurones for the hidden layers of the intermediate autoencoders', type=str, default='1024,512,256')
parser.add_argument('-ct', '--central_dim', help='list of neurones for the hidden layers of the central autoencoder', type=str, default='2048,1024,512,256')
parser.add_argument('-lt', '--latent_dim', help='size of the latent vector', type=int, default=128)

parser.add_argument('-ch', '--classifier_dim', help='list of neurones for the classifier hidden layers', type=str, default='256,128')
parser.add_argument('-sh', '--survival_dim', help='list of neurones for the survival hidden layers', type=str, default='64,32')
parser.add_argument('-lc', '--lambda_classif', help='weight of the classification loss', type=float, default=5)
parser.add_argument('-ls', '--lambda_survival', help='weight of the survival loss', type=float, default=0)

parser.add_argument('-xp', '--explain', help='choose if you want to explain the results or not', type=bool, default=False)
parser.add_argument('-cxp', '--explained_class', help='class to explain', type=str, default='Her2')
parser.add_argument('-sxp', '--explained_source', help='source to explain', type=str, default='RNAseq')

args = parser.parse_args()

if args.device != 'cpu':
    cuda = torch.cuda.is_available()
    device = torch.device("cuda" if cuda else "cpu")
else:
    device = torch.device("cpu")

sources = args.sources.split(',')
hidden_dim = [int(element) for element in args.hidden_dim.split(',')]
central_dim = [int(element) for element in args.central_dim.split(',')]
classifier_dim = [int(element) for element in args.classifier_dim.split(',')]
survival_dim = [int(element) for element in args.survival_dim.split(',')]

if __name__ == "__main__":
    lt_metrics = []
    for split in range(1, 6):
        metric = train(args.task, args.cohorts, sources, split, device, num_classes=args.num_classes,
             batch_size=args.batch_size, n_epochs=args.epochs, beta=args.beta, lr=args.lr,
            hidden_dim=hidden_dim, central_dim=central_dim, latent_dim=args.latent_dim, dropout=args.dropout,
            classifier_dim=classifier_dim, survival_dim=survival_dim, lambda_classif=args.lambda_classif, 
            lambda_survival=args.lambda_survival, explain=args.explain, explained_source=args.explained_source, explained_class=args.explained_class)
        lt_metrics.append(metric)
        print(metric)
    print('C-index : {} +- {}'.format(np.mean(lt_metrics), np.std(lt_metrics)))