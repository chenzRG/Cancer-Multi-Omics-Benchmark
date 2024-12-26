# -*- coding: utf-8 -*-
"""
Created on Wed 01 Sept 2021

@author: Hakim Benkirane

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Conception of the consensus integration loss.
"""
import torch.nn as nn


def consensus_loss(x, autoencoders):
    n_autoencoders = len(autoencoders)
    loss = 0
    reconstruction_loss = nn.MSELoss()
    for i in range(n_autoencoders):
        for j in range(n_autoencoders):
            loss += reconstruction_loss(autoencoders[j].decoder(autoencoders[i].encoder(x[i])), x[j])
    return loss