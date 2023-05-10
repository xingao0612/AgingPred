# AgingPred

This is an R package developed based on Shiny app, which allows users to upload peripheral blood leukocyte expression profile data. It uses a pre-built random forest age prediction model to predict the age of the samples and identify the aging rate.<br><br>

How to install:
```R
library(devtools)
install_github("xingao0612/AgingPred")
```
Note: before installing "AgingPred", please make sure that you have installed the "sva" package. If not, install the "sva" package with the following code:
```R
BiocManager::install("sva")
```
AgingPred usage:<br>
```R
library(AgingPred)
AgingPred()
```
Start interface
![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface.jpg)  
Result interface
![](https://github.com/xingao0612/AgingPred/blob/master/inst/www/interface2.jpg)  
