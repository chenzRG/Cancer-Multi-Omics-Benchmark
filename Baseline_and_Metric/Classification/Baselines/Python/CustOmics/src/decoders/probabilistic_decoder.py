# -*- coding: utf-8 -*-
"""
Created on Wed 01 Sept 2021

@author: Hakim Benkirane

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Build the Probabilistic Decoder module.
"""

import torch
import torch.nn as nn
from collections import OrderedDict
from src.tools.net_utils import FullyConnectedLayer


class ProbabilisticDecoder(nn.Module):
    """
    Module that performs the generation step for the variational autoencoder
    """
    def __init__(self, latent_dim, hidden_dim, output_dim, norm_layer=nn.BatchNorm1d, leaky_slope=0.2, dropout=0, debug = False):
        """
        Constructs the generative network for the VAE architecture
        Parameters
        ----------
        latent_dim: int
            Dimension of the latent representation.
        hidden_dim: list
            List of dimensions for the multiple intermediate layers.
        output_dim: int
            Dimension of the output layer.
        norm_layer: pytorch.nn
            Normalization Layer.
        leaky_slope: float
            Coefficient for the Leaky ReLU (must be between 0 and 1).
        dropout: float
            Dropout rate (must be between 0 and 1).
        debug: bool
            Debug parameter, prints the intermediate tensors during training.
        """
        super(ProbabilisticDecoder, self).__init__()

        self.dt_layers = OrderedDict()

        self.dt_layers['InputLayer'] = FullyConnectedLayer(latent_dim, hidden_dim[0], norm_layer=norm_layer, leaky_slope=leaky_slope, dropout=dropout,
                                activation=True)

        block_layer_num = len(hidden_dim)
        dropout_flag = True
        for num in range(1, block_layer_num):
            self.dt_layers['Layer{}'.format(num)] = FullyConnectedLayer(hidden_dim[num - 1], hidden_dim[num], norm_layer=norm_layer, leaky_slope=leaky_slope,
                                    dropout=dropout_flag*dropout, activation=True)
            # dropout for every other layer
            dropout_flag = not dropout_flag

        # the output fully-connected layer of the classifier
        self.dt_layers['OutputLayer'] = FullyConnectedLayer(hidden_dim[-1], output_dim, norm_layer=norm_layer, leaky_slope=leaky_slope, dropout=0,
                                 activation=False, normalization=False)

        self.net = nn.Sequential(self.dt_layers)
        
    def forward(self, x):
        
        x_hat = torch.sigmoid(self.net(x))
        return x_hat