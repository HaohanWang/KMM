# Main file for usage of CMM (Couple Mixed Model)
# If you find this work useful, please cite:
# Wang H., Pei F., Vanyukov MM., Bahar I., Wu W. and Xing EP.
# Coupled Mixed Model for Joint GWAS of Pairs of Complex Disorders with Independently Collected Datasets
# Contact: {haohanw,weiwu2,epxing}@cs.cmu.edu

import sys
from model.KMM import KernelMixedModel
from util.dataLoader import *
from util.supportingFunctions import *

def printOutHead(out): out.write("\t".join(["RANK", "SNP_ID", "EFFECT_SIZE_ABS"]) + "\n")

def outputResult(out, rank, id, beta):
    out.write("\t".join([str(x) for x in [rank, id, beta]]) + "\n")


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
dataGroup.add_option("--cov", dest='covariateFile', help="name of the covariate file")
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

reader = FileReader(fileName=options.fileName1, fileType=options.fileType, imputation=(not options.missing))
X1, Y1, Xname1 = reader.readFiles()

reader = FileReader(fileName=options.fileName2, fileType=options.fileType, imputation=(not options.missing))
X2, Y2, Xname2 = reader.readFiles()


model = KernelMixedModel(quiet=options.quiet, fdr=options.fdr)

print ('Computation starts ... ')

if options.snum is not None:  # select for a fix number of variable
    snum = float(options.snum)
    Beta1, Beta2 = binarySearch_Rho(model, X1, Y1, X2, Y2, snum, quiet=options.quiet)
elif options.lmbd is not None:
    lmbd = float(options.lmbd)
    model.setLambda(lmbd)
    model.setRho(float(options.rho))
    model.setLearningRate(learningRate)  # learning rate must be set again every time we run it.
    model.fit(X1, Y1, X2, Y2)
    Beta1 = model.getBeta1()
    Beta2 = model.getBeta2()
else:
    Beta1, Beta2 = crossValidation(model, X1, Y1, X2, Y2, learningRate)

ind1 = np.where(Beta1 != 0)[0]
bs1 = Beta1[ind1].tolist()
xname1 = []
for i in ind1:
    xname1.append(i)

beta_name1 = zip(bs1, xname1)
bn = sorted(beta_name1)
bn.reverse()

out1 = open(outFile1, 'w')
printOutHead(out1)

for i in range(len(bn)):
    outputResult(out1, i + 1, bn[i][1], bn[i][0])

out1.close()

ind2 = np.where(Beta2 != 0)[0]
bs2 = Beta2[ind2].tolist()
xname2 = []
for i in ind2:
    xname1.append(i)

beta_name2 = zip(bs2, xname2)
bn = sorted(beta_name2)
bn.reverse()

out2 = open(outFile2, 'w')
printOutHead(out2)

for i in range(len(bn)):
    outputResult(out2, i + 1, bn[i][1], bn[i][0])

out1.close()

print ('\nComputation ends normally, check the output file at ', outFile)