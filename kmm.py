# Main file for usage of CMM (Couple Mixed Model)
# If you find this work useful, please cite:
# Wang H., Pei F., Vanyukov MM., Bahar I., Wu W. and Xing EP.
# Coupled Mixed Model for Joint GWAS of Pairs of Complex Disorders with Independently Collected Datasets
# Contact: {haohanw,weiwu2,epxing}@cs.cmu.edu

import sys
from model.KMM import KernelMixedModel
from util.dataLoader import *
from util.supportingFunctions import *

from optparse import OptionParser, OptionGroup

usage = """usage: %prog [options] --file1 fileName1 --file2 fileName2
This program provides the basic usage to CMM, e.g:
python cmm.py --file1 data/mice1.plink --file2 data/mice2.plink
	    """
parser = OptionParser(usage=usage)

dataGroup = OptionGroup(parser, "Data Options")
modelGroup = OptionGroup(parser, "Model Options")

## data options
dataGroup.add_option("--gene", dest='geneFile', help="name of the gene expression file")
dataGroup.add_option("--pheno", dest='phenoFile', help="name of the phenotype file")
dataGroup.add_option("--cov", dest='covariateFile', default=None, help="name of the covariate file")
dataGroup.add_option("--net", dest='networkFile', help="name of the network file")
dataGroup.add_option("--out", dest='outputFile', default='output', help="name of the output file")

## model options
modelGroup.add_option("--threshold", dest="threshold", default=0.5,
                      help="the threshold of the correlation network")
modelGroup.add_option('-f', action='store_true', dest='fdr', default=False, help='to use FDR')
modelGroup.add_option('-q', action='store_true', dest='quiet', default=False, help='Run in quiet mode')

## advanced options
parser.add_option_group(dataGroup)
parser.add_option_group(modelGroup)

(options, args) = parser.parse_args()


def outputFiles():
    pass

fileType = 0
IN = None

if len(args) != 0:
    parser.print_help()
    sys.exit()

outFile = options.outputFile + '.txt'

print ('Running ... ')

X, genes = readFiles(options.geneFile)
y, _ = readFiles(options.phenoFile)
if options.covariateFile is not None:
    c, _ = readFiles(options.covariateFile)
    X = covariateRegression(X, c)
N = readNetworkFiles(options.networkFile, genes)


model = KernelMixedModel(quiet=options.quiet, fdr=options.fdr)

print ('Computation starts ... ')

L = constructLaplacian(N)
J = approximateInverseLaplacian(L)
C = calculateCorrelationMatrix(X)

K = np.dot(np.dot(X, J), X.T)
model.setJ(J)
model.setK(K)
model.setC(C)
model.setN(N)
model.setCovarianceThreshold(options.threshold)

model.fit(X, y)
neg_log_p = model.getNegLogP()
pvalue = np.exp(-neg_log_p)
beta = model.getBeta()

idx = np.argsort(pvalue)

f = open(options.outputFile, 'w')
f.writelines('gene\tp-value\teffect size')
for i in range(len(idx)):
    if pvalue[idx[i]] < 0.05:
        f.writelines(genes[idx[i]] + '\t' + str(pvalue[idx[i]]) + '\t' + str(beta[idx[i]]) + '\n')

f.close()