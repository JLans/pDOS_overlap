# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import numpy as np
from copy import deepcopy

def get_r2(y_true, y_pred):
    """R2 or the error.

    Parameters
    ----------
    y_true : numpy.ndarray or list
        Ground truth (correct) values

    y_pred : numpy.ndarray or list
        Predicted values, as returned by a regression estimator.

    Returns
    -------
    loss : float
        R2 value.
    """
    SStot = np.sum((y_true-y_true.mean())**2)
    SSres = np.sum((y_true-y_pred)**2)
    return 1 - SSres/SStot

def get_rmse(y_true, y_pred):
    """Compute maximum absolute error.

    Parameters
    ----------
    y_true : numpy.ndarray or list
        Ground truth (correct) values.

    y_pred : numpy.ndarray or list
        Predicted values, as returned by a regression estimator.

    Returns
    -------
    loss : float
        The maximum absolute error times the sign of the error.
    """
    SSres = np.mean((y_true-y_pred)**2)
    return SSres**0.5

def get_max_error(y_true, y_pred):
    """Compute maximum absolute error.

    Parameters
    ----------
    y_true : numpy.ndarray or list
        Ground truth (correct) values.

    y_pred : numpy.ndarray or list
        Predicted values, as returned by a regression estimator.

    Returns
    -------
    loss : float
        The maximum absolute error.
    """
    return np.array(y_pred-y_true)[np.argmax(np.abs(y_pred-y_true))]

def get_wasserstein_loss(y_true, y_pred,individual=False):
    """Compute the l2 wasserstein loss

    Parameters
    ----------
    y_true : numpy.ndarray or list
        Ground truth (correct) values.

    y_pred : numpy.ndarray or list
        Predicted values, as returned by a regression estimator.
        
    individual : bool
        If True, returns Wasserstein loss along first dimension.
        If False, returns average Wasserstein loss.

    Returns
    -------
    loss : float or numpy.ndarray
        The degree to which the samples are correctly predicted.
    """
    y_true = np.array(deepcopy(y_true))
    y_pred = np.array(deepcopy(y_pred))
    if len(y_true.shape) == 1:
        y_true = y_true.reshape(-1,1)
    if len(y_pred.shape) == 1:
        y_pred = y_pred.reshape(-1,1)
    Tcum = np.cumsum(y_true, axis=-1)
    Pcum = np.cumsum(y_pred, axis=-1)
    w_loss = (1/float(y_true.shape[1])*np.sum((Pcum-Tcum)**2, axis=1))**0.5
    if individual == True:
        return w_loss
    else:
        return w_loss.mean()
