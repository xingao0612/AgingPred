# AgingPred

This is an R package developed based on Shiny app, which allows users to upload peripheral blood leukocyte expression profile data. It uses a pre-built random forest age prediction model to predict the age of the samples and identify the aging rate.<br><br>

How to install:
```R
library(devtools)
install_github("xingao0612/AgingPred")
```
注意：在安装“AgingPred”之前，请确认您已经安装了"sva"包，如没有安装，请先使用以下代码安装"sva"包： 
```R
BiocManager::install("sva")
```
usage:<br>
```R
library(AgingPred)
AgingPred()
```
![](https://github.com/xingao0612/AgingPred/blob/master/interface.jpg)  
