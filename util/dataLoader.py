__author__ = 'Haohan Wang'

import numpy as np

def readFiles(fileName, skipFirstColumn=True, splitter=','):
    '''
    :param fileName: the fileName for phenotype, gene expression, and covariates
    the file is structured as the first row is the name of variables
    ......................... the first column is the name of the id of the samples
    ......................... the remaining are the values of the variables
    :param skipFirstColumn: whether the firstColumn holds the sample IDs,
    ......................... set False if the input has sample IDs
    :return: npy matrix, variable names (list)
    ......... note the order of the matrix and the variable list must be consistent
    '''

    text = [line.strip() for line in open(fileName)]

    if skipFirstColumn:
        startingCol = 0
    else:
        startingCol = 1

    ColNames = text[0][startingCol:]

    data = []

    for i in range(1, len(text)):
        row = [float(t) for t in text[i].split(splitter)[startingCol:]]
        data.append(row)

    data = np.array(data)

    return data, ColNames

def readNetworkFiles(fileName, names, splitter=','):
    '''
    :param fileName: the fileName of the network file
    the file is structured as each line is starts with a gene and follows with the gene(s) it interacts
    :param names: the geneNames returned from the fileReading
    :return: the matrix sparse of network
    '''
    network = np.zeros([len(names), len(names)])
    n2i = {}
    for i in range(len(names)):
        n2i[names[i].lower()] = i

    text = [line.strip() for line in open(fileName)]
    for line in text:
        items = line.split(splitter)
        n1 = items[0].lower()
        for item in items[1:]:
            n2 = item.lower()
            if n1 in n2i and n2 in n2i:
                network[n2i[n1], n2i[n2]] = 1
                network[n2i[n2], n2i[n1]] = 1

    return network

