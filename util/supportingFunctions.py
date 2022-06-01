__author__ = 'Haohan Wang'

import numpy as np
from scipy.stats import spearmanr

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

def calculateCorrelationMatrix(X):
    C = np.zeros([X.shape[1], X.shape[1]])
    for i in range(X.shape[1]):
        C[i, i] = 1
        for j in range(i + 1, X.shape[1]):
            C[i, j] = spearmanr(X[:, i], X[:, j])[0]
            C[j, i] = C[i, j]

    return C

def covariateRegression(X, c):
    X_new = np.zeros_like(X)

    tmp = np.dot(np.linalg.inv(np.dot(c.T, c)), c.T)
    for i in range(X.shape[1]):
        X_new[:, i] = X[:, i] - np.dot(c, np.dot(tmp, X[:, i]))

    return X_new

