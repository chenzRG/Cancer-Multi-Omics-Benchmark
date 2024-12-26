# -*- coding: utf-8 -*-
"""
Created on Wed 01 Sept 2021

@author: Hakim Benkirane

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Build the Standard Autoencoder module.
"""
import torch
import torch.nn as nn


class AutoEncoder(nn.Module):
    """
    A standard autoencoder on a signle-source input (it serves as a base model for more complex architectures)
    """
    def __init__(self, encoder, decoder, device):
        """
        Constructs the variational autoencoder's (VAE) architecture
        Parameters:
            encoder (Encoder) -- the encoder of the model
            decoder (Decoder) -- the decoder of the model
            device (pytorch)  -- the device in which the computation is done
        """
        super(AutoEncoder, self).__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.device = device
        self._relocate()

    def _relocate(self):
        self.encoder.to(self.device)
        self.decoder.to(self.device)
        
                
    def forward(self, x):
        z = self.encoder(x)  
        x_hat            = self.decoder(z)
        return x_hat, z

    def decode(self, z):
        x_hat = self.decoder(z)
        return x_hat
        
    def loss(self, x, beta):
        x_hat, z = self.forward(x)
        reconstruction_loss = nn.MSELoss()

        return reconstruction_loss(x, x_hat)




