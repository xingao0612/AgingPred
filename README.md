# AgingPred package

This is an R package developed based on Shiny app, which allows users to upload peripheral blood leukocyte expression profile data. It uses a pre-built random forest age prediction model to predict the age of the samples and identify the aging rate.<br><br>

## Installation
Note: Before installing this App, you will need to install some dependent R packages on your R.
```R
if (!require("BiocManager"))
  install.packages("BiocManager")
library(BiocManager)
if (!require("limma"))
  BiocManager::install("sva")
if (!require("devtools"))
  install.packages("devtools")
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
Start interface
![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface.jpg)  
Result interface
![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface2.jpg)  
