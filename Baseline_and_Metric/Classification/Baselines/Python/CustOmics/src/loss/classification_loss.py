# -*- coding: utf-8 -*-
"""
Created on Wed 01 Sept 2021

@author: Hakim Benkirane

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Conception of the classification loss.
"""

import torch.nn as nn



def classification_loss(loss_name, y_true, y_pred ,reduction='mean'):
    """
    Return the loss function.
    Parameters:
        loss_name (str)    -- the name of the loss function: BCE | MSE | L1 | CE
        reduction (str)    -- the reduction method applied to the loss function: sum | mean
    """
    if loss_name == 'BCE':
        return nn.BCEWithLogitsLoss(reduction=reduction)(y_true, y_pred)
    elif loss_name == 'CE':
        return nn.CrossEntropyLoss(reduction=reduction)(y_true, y_pred)
    else:
        raise NotImplementedError('Loss function %s is not found' % loss_name)