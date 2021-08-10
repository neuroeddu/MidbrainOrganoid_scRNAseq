# install and update packages

install.packages("gglot2","Seurat","patchwork","dplyr","Matrix","limma")

# install devtools and install scCATCH
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
# install library

# install clustifyr

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("clustifyr")

# import all the libraries

library("ggplot2")
library("Seurat")
library("cowplot")
library("clustree")
library(patchwork)
library(dplyr)
library("Matrix")
library("scCATCH")

