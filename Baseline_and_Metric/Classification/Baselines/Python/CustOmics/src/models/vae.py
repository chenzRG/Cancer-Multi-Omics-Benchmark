# -*- coding: utf-8 -*-
"""
Created on Wed 01 Sept 2021

@author: Hakim Benkirane

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Build the Variational Autoencoder module.
"""

import torch
import torch.nn as nn

from src.loss.mmd_loss import compute_mmd


class VAE(nn.Module):
    """
    Module that performs variational inference on a single-source input (it serves as a base model for more complex architectures)
    """
    def __init__(self, encoder, decoder, device):
        """
        Constructs the variational autoencoder's (VAE) architecture
        Parameters:
            encoder (ProbabilisticEncoder) -- the encoder of the model
            decoder (ProbabilisticDecoder) -- the decoder of the model
            device (pytorch)               -- the device in which the computation is done
        """
        super(VAE, self).__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.device = device
        self._relocate()

    def _relocate(self):
        self.encoder.to(self.device)
        self.decoder.to(self.device)
        
    def reparameterization(self, mean, var):
        epsilon = torch.randn_like(var).to(self.device)        # sampling epsilon        
        z = mean + var*epsilon                          # reparameterization trick
        return z
        
                
    def forward(self, x):
        mean, log_var = self.encoder(x)
        z = self.reparameterization(mean, torch.exp(0.5 * log_var)) # takes exponential function (log var -> var)
        x_hat           = self.decoder(z)
        return x_hat, z

    def loss(self, x, beta):
        x_hat, z = self.forward(x)

        reconstruction_loss = nn.MSELoss()
        recon = reconstruction_loss(x, x_hat)

        true_samples=torch.randn(z.shape[0], z.shape[1]).to(self.device)
        MMD=torch.sum(compute_mmd(true_samples, z))

        return recon + beta * MMD


