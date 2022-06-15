__author__ = 'Haohan Wang'

import scipy.optimize as opt
import scipy
import numpy as np
import numpy.linalg as linalg

import sys
sys.path.append('../')

import time

from helpingMethods import *

class KernelMixedModel:
    def __init__(self, numintervals=100, ldeltamin=-5, ldeltamax=5, scale=0, alpha=0.05, fdr=False, covT=0.5, quiet=True):
        self.numintervals = numintervals
        self.ldeltamin = ldeltamin
        self.ldeltamax = ldeltamax
        self.scale = scale
        self.alpha = alpha
        self.fdr = fdr
        self.ldelta0 = None
        self.covT = covT
        self.quiet = quiet

    def estimateDelta(self, y):
        if y.ndim == 1:
            y = y.reshape((y.shape[0], 1))
        self.y = y
        self.S, self.U, self.ldelta0= self.train_nullmodel(y, self.K)
        print(self.ldelta0)

    def setK(self, K):
        self.K = K

    def setLDelta(self, delta0):
        self.ldelta0=delta0

    def setCovarianceThreshold(self, covT):
        self.covT = covT

    def setJ(self, J):
        self.J = J

    def setC(self, C):
        self.C = C

    def setN(self, N):
        self.N = N

    def fit(self, X, y):
        [n_s, n_f] = X.shape
        K_all = np.copy(self.K)

        self.neg_log_p = np.zeros(n_f)
        self.beta = np.zeros(n_f)

        stepSize = n_f*0.05
        count = 0
        countStep = 1

        start = time.time()

        for i in range(n_f):

            idx1 = np.where(np.abs(self.C[i]) > self.covT)[0]
            idx2 = np.where(self.N[i] == 0)
            idx = np.intersect1d(idx1, idx2)

            K = np.dot(np.dot(X[:, idx].reshape([n_s, idx.shape[0]]), self.J[idx,:].reshape([idx.shape[0], n_f])), X.T)
            self.setK(K_all - K)

            if not self.quiet:
                count += 1
                if count > countStep*stepSize:
                    print ('Finished ', countStep*5, '% of the progress, ', end = '')
                    print ('with', time.time() - start, 'seconds')
                    countStep += 1

            try:
                self.neg_log_p[i], self.beta[i] = self.fit_helper(X[:, i].reshape([n_s, 1]), y)
            except:
                self.neg_log_p[i] = 0
                self.beta[i] = 0

    def fit_helper(self, X, y):
        [n_s, n_f] = X.shape
        if y.ndim == 1:
            y = y.reshape((y.shape[0], 1))
        self.y = y
        self.S, self.U, ldelta0= self.train_nullmodel(y, self.K)

        # if self.ldelta0 is None:
        #     self.ldelta0 = ldelta0
        self.ldelta0 = ldelta0

        X0 = np.ones(len(self.y)).reshape(len(self.y), 1)

        delta0 = np.exp(self.ldelta0)
        Sdi = 1. / (self.S + delta0)
        Sdi_sqrt = np.sqrt(Sdi)
        SUX = np.dot(self.U.T, X)
        SUX = SUX * np.tile(Sdi_sqrt, (n_f, 1)).T
        SUy = np.dot(self.U.T, self.y)
        SUy = SUy * Sdi_sqrt.reshape((n_s, 1))

        self.SUX = SUX
        self.SUy = SUy
        SUX0 = np.dot(self.U.T, X0)
        self.SUX0 = SUX0 * np.tile(Sdi_sqrt, (1, 1)).T
        self.X0 = X0
        neg_log_p, beta = self.hypothesisTest(self.SUX, self.SUy, X, self.SUX0, self.X0)
        return neg_log_p, beta

    def fdrControl(self):
        tmp = np.exp(-self.neg_log_p)
        tmp = sorted(tmp)
        threshold = 1e-8
        n = len(tmp)
        for i in range(n):
            if tmp[i] < (i+1)*self.alpha/n:
                threshold = tmp[i]
        self.neg_log_p[self.neg_log_p<-np.log(threshold)] = 0

    def getNegLogP(self):
        if not self.fdr:
            return self.neg_log_p
        else:
            self.fdrControl()
            return self.neg_log_p

    def hypothesisTest(self, UX, Uy, X, UX0, X0):
        [m, n] = X.shape
        p = []
        for i in range(n):
            if UX0 is not None:
                UXi = np.hstack([UX0, UX[:, i].reshape(m, 1)])
                XX = matrixMult(UXi.T, UXi)
                XX_i = linalg.pinv(XX)
                beta = matrixMult(matrixMult(XX_i, UXi.T), Uy)
                Uyr = Uy - matrixMult(UXi, beta)
                Q = np.dot(Uyr.T, Uyr)
                sigma = Q * 1.0 / m
            else:
                Xi = np.hstack([X0, UX[:, i].reshape(m, 1)])
                XX = matrixMult(Xi.T, Xi)
                XX_i = linalg.pinv(XX)
                beta = matrixMult(matrixMult(XX_i, Xi.T), Uy)
                Uyr = Uy - matrixMult(Xi, beta)
                Q = np.dot(Uyr.T, Uyr)
                sigma = Q * 1.0 / m
            ts, ps = tstat(beta[1], XX_i[1, 1], sigma, 1, m)
            if -1e100 < ts < 1e100:
                p.append(ps)
            else:
                p.append(1)
        p = np.array(p)
        return -np.log(p), beta[1]

    def train_nullmodel(self, y, K, S=None, U=None, numintervals=500):
        self.ldeltamin += self.scale
        self.ldeltamax += self.scale

        if S is None or U is None:
            S, U = linalg.eigh(K)

        Uy = np.dot(U.T, y)

        # grid search
        nllgrid = np.ones(numintervals + 1) * np.inf
        ldeltagrid = np.arange(numintervals + 1) / (numintervals * 1.0) * (self.ldeltamax - self.ldeltamin) + self.ldeltamin
        for i in np.arange(numintervals + 1):
            nllgrid[i] = nLLeval(ldeltagrid[i], Uy, S) # the method is in helpingMethods

        nllmin = nllgrid.min()
        ldeltaopt_glob = ldeltagrid[nllgrid.argmin()]

        for i in np.arange(numintervals - 1) + 1:
            if (nllgrid[i] < nllgrid[i - 1] and nllgrid[i] < nllgrid[i + 1]):
                ldeltaopt, nllopt, iter, funcalls = opt.brent(nLLeval, (Uy, S),
                                                              (ldeltagrid[i - 1], ldeltagrid[i], ldeltagrid[i + 1]),
                                                              full_output=True)
                if nllopt < nllmin:
                    nllmin = nllopt
                    ldeltaopt_glob = ldeltaopt
        return S, U, ldeltaopt_glob

    def getBeta(self):
        return self.beta