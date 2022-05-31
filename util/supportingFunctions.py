__author__ = 'Haohan Wang'

import numpy as np

def constructLaplacian(network):
    L = np.zeros_like(network)
    for i in range(network.shape[0]):
        for j in range(network.shape[1]):
            if i == j and np.sum(network[i])!=0:
                L[i, j] = 1
            elif network[i, j] == 1:
                L[i, j] = -1.0/np.sqrt(np.sum(network[i])*np.sum(network[j]))

def approximateInverseLaplacian(L):
    I = np.eye(L.shape[0])

    H = L - I
    T = np.dot(H, H)
    J = I + H + T + np.dot(T, H)

    return J

def calculateCorrelationMatrix():
    pass

def covariateRegression():
    pass

