# AgingPred
![](https://img.shields.io/badge/source%20code-support-blue) ![](https://img.shields.io/badge/R-package-green) <br><br>
## Introduction
This is an R package developed based on Shiny app, which allows users to upload peripheral blood leukocyte expression profile data. It uses a pre-built random forest age prediction model to predict the age of the samples and identify the aging rate.<br><br>

## Installation
Note: Before installing this App, you will need to install some dependent R packages on your R.
```R
if (!require("BiocManager"))
  install.packages("BiocManager")
if (!require("devtools"))
  install.packages("devtools")
library(BiocManager)
if (!require("limma"))
  BiocManager::install("sva")
  
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

![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface.jpg)

Result interface
![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface2.jpg)  

## Contact us
E-mail: gaoxin_0612@163.com

