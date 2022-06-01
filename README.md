![KMM](pl.PNG "KMM")

# Kernel Mixed Model

Implementation of the Kernel Mixed Model in the following papers:

   Wang, H., Lopez, O. L., Wu, W., & Xing, E. P. (2022, April). Gene Set Priorization Guided by Regulatory Networks with p-values through Kernel Mixed Model. In Research in Computational Molecular Biology: 26th Annual International Conference, RECOMB 2022, San Diego, CA, USA, May 22–25, 2022, Proceedings (pp. 107-125). Cham: Springer International Publishing.
   
   and others under review

## Introduction

The Kernel Mixed Model incorporates the network structure to improve the study of transcriptome association study

**Replication:** This repository serves for the purpose to guide others to use our tool, if you are interested in the scripts to replicate our results, please contact us and we will share the repository for replication. Contact information is at the bottom of this page. 

## File Structure:

* model/ main method for the Kernel Mixed Model
* util/ other helper files
* kmm.py main entry point of using the Kernel Mixed Model to work with your own data

## An Example Command:

```
python kmm.py -gene <expression values> -pheno <phenotype values> -net <network structure> -cov <covariate values to regress out> -out <output file> 
```
#### Data Support
* The implementation currently supports CSV files. 
* Extensions to other data format can be easily implemented through `FileReader` in `util/dataLoadear`. Feel free to contact us for the support of other data format. 

## Installation (Not Required)
You will need to have numpy and scipy installed on your current system.
You can install precision lasso using pip by doing the following

```
   pip install git+https://github.com/HaohanWang/KMM
```

You can also clone the repository and do a manual install.
```
   git clone https://github.com/HaohanWang/KMM
   python setup.py install
```


## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
&middot;
[@HaohanWang](https://twitter.com/HaohanWang)
