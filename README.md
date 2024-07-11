# primate-brains-proximate-ultimate

#### Repository for the joint analysis of proximate and ultimate factors on primate brain size

This repository contains scripts used in the joint analysis of proximate and ultimate factors on primate brain size.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

* A .csv file including phenotype data ```data/brain_data.csv```
  
# Pipeline
  
### Read dN/dS data

* **Key libraries:** phylobase

```
# Import metadata and combine
scripts/read_dnds.R
```
