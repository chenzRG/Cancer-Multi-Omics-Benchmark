# -*- coding: utf-8 -*-
"""
Created on Wed 01 Sept 2021

@author: Hakim Benkirane

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Build the Multi-Layer Classifier netowrk.
"""
import torch
import torch.nn as nn
from collections import OrderedDict
from src.tools.net_utils import FullyConnectedLayer
from src.loss.classification_loss import classification_loss

from torch.optim import Adam


class MultiClassifier(nn.Module):
    """
    Defines a multi-layer fully-connected classifier
    """
    def __init__(self, n_class=2, latent_dim=256, norm_layer=nn.BatchNorm1d, leaky_slope=0.2, dropout=0,
                 class_dim = [128, 64]):
        """
        Construct a multi-layer fully-connected classifier
        Parameters:
            class_num (int)         -- the number of class
            latent_dim (int)        -- the dimensionality of the latent space and the input layer of the classifier
            norm_layer              -- normalization layer
            leaky_slope (float)     -- the negative slope of the Leaky ReLU activation function
            dropout_p (float)       -- probability of an element to be zeroed in a dropout layer
            layer_num (int)         -- the layer number of the classifier, >=3
        """
        super(MultiClassifier, self).__init__()

        self.dt_layers = OrderedDict()

        self.dt_layers['InputLayer'] = FullyConnectedLayer(latent_dim, class_dim[0], norm_layer=norm_layer, leaky_slope=leaky_slope, dropout=dropout,
                                activation=True)

        block_layer_num = len(class_dim)
        dropout_flag = True
        for num in range(1, block_layer_num):
            self.dt_layers['Layer{}'.format(num)] = FullyConnectedLayer(class_dim[num - 1], class_dim[num], norm_layer=norm_layer, leaky_slope=leaky_slope,
                                    dropout=dropout_flag*dropout, activation=True)
            # dropout for every other layer
            dropout_flag = not dropout_flag

        # the output fully-connected layer of the classifier
        self.dt_layers['OutputLayer'] = FullyConnectedLayer(class_dim[-1], n_class, norm_layer=norm_layer, leaky_slope=leaky_slope, dropout=0,
                                 activation=False, normalization=False)

        self.net = nn.Sequential(self.dt_layers)

    def forward(self, x):
        y = self.net(x)
        return y

    def predict(self, x):
        return torch.max(self.forward(x), dim=1).indices

    def compile(self, lr):
        self.optimizer = Adam(self.parameters(), lr=lr)

    def fit(self, x_train, y_train, n_epochs, verbose=False):
        self.train()
        for epoch in range(n_epochs):
            overall_loss = 0
            loss = 0
            self.optimizer.zero_grad()
            y_pred = self.forward(x_train)
            loss = classification_loss('CE', y_pred, y_train)
            overall_loss += loss.item()
        
            loss.backward()
            self.optimizer.step()
            if verbose:
                print("\tEpoch", epoch + 1, "complete!", "\tAverage Loss: ", overall_loss)