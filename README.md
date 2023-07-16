# AgingPred
![](https://img.shields.io/badge/source%20code-support-blue) ![](https://img.shields.io/badge/R-package-green) ![](https://img.shields.io/badge/Version-0.1.2-yellow)<br>
## Introduction
This is an R package developed using Shiny that allows users to upload human peripheral blood RNA-seq data and predict the age and aging rate of the uploaded samples using our pre-built random forest age prediction model.<br>

## Installation
Note: Before installing this App, you will need to install some dependent R packages on your R.<br>
From Bioconductor:
```R
if (!require("BiocManager"))
  install.packages("BiocManager")
if (!require("devtools"))
  install.packages("devtools")
library(BiocManager)
if (!require("limma"))
  install.packages("limma")
if (!require("sva"))
  install.packages("sva")
  ```
  From CRAN (Comprehensive R Archive Network): 
  ```R
  install.packages("shiny")
  install.packages("shinydashboard")
  install.packages("ggplot2")
  install.packages("ggpubr")
  install.packages("shinycssloaders")
  install.packages("randomForest")
  install.packages("dplyr")

```
How to install AgingPred:
```R
devtools::install_github("xingao0612/AgingPred")
```

## AgingPred usage:<br>
```R
library(AgingPred)
AgingPred()
```
## APP Interface

Start interface

![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface11.png)

Result interface
![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/results inter.jpg)  

## Contact us
E-mail: gaoxin_0612@163.com

