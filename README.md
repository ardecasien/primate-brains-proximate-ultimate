# primate-brains-proximate-ultimate

#### Repository for the joint analysis of proximate and ultimate factors on primate brain size

This repository contains scripts used in the joint analysis of proximate and ultimate factors on primate brain size.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

* A .csv file including phenotype data ```data/phenotype_data.csv```
* A .csv file including dN/dS data ```data/dnds.csv```
* A .nex file including the phylogenetic tree data ```data/consensusTree_10kTrees_Primates_Version3.nex```
  
# Pipeline

* **Key libraries:** mvMORPH, phylopath, phytools

```
scripts/teja-et-al-code.R
scripts/teja-et-al-functions.R

```
