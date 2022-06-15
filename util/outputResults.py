__author__ = 'Haohan Wang'

import numpy as np

def outputResults(pvalue, genes, betas, fileName):

    idx = np.argsort(pvalue)

    f = open(fileName + '.html', 'w')

    header = [line.strip() for line in open('resultHead.txt')]
    for line in header:
        f.writelines(line + '\n')

    for i in range(len(idx)):
        if pvalue[idx[i]] < 0.05:
            f.writelines('<tr>\n')
            f.writelines('<th scope="row">' + genes[idx[i]] +'</th>\n')
            f.writelines('<td>'+ str(round(pvalue[idx[i]], 6)) + '</td>\n')
            f.writelines('<td>' + str(betas[idx[i]]) + '</td>\n')
            f.writelines('</tr>\n')


    tail = [line.strip() for line in open('resultTail.txt')]
    for line in tail:
        f.writelines(line + '\n')

    f.close()

