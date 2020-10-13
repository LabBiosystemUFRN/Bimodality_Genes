# <center> Bimodal Gene Expression Detection</center>
## <center> A new method to identify genes with bimodal gene expression</center>

####  <span style="color:red">The code located here is provided "AS-IS", with no warranties, express or implied, and hereby disclaims all implied warranties, including any warranty of merchantability and warranty of fitness for a particular purpose. It's released under GPL v2.</span> </center>
<br>

This repository contains the files necessary to reproduce the results reported in our analysis of bimodality in gene expression. The *bin* folder contains the code used to generate analysis and figures. It's based on R software (version 3.6.3 2020-02-29) and uses the libraries bnlearn, signal, ggplot2, mclust, gridExtra, grid,lattice, cowplot, reshape2 and doParallel.

The *data_files*  folder contains data from 522 patients with primary tumor of clear cell renal carcinoma (KIRC).

The original data is splited in two files: KIRC_FPKM_TP.tsv.tar.bz.part.aa and KIRC_FPKM_TP.tsv.tar.bz.part.ab

The script unpack.sh reassmeble and stract the data

<a name="any"></a>
###  Our pipeline  
To run our pipeline, download our files from [github](https://github.com/LabBiosystemUFRN/Bimodality_Genes) 
> mkdir yourFolder
> cd yourFolder
> git clone https://github.com/LabBiosystemUFRN/Bimodality_Genes

Then open the file yourFolber/bin/bimodalityDetection.R, locate the line
> dirBase<-"Place here the correct name of yourFolder"

and replace it with your corret folder name, like

> dirBase<-"/home/smith/yourFolder"

Then you can run :
> Rscript bin/bimodalityDetection.R

### R dependencies

> bnlearn, signal, ggplot2, mclust, gridExtra, grid,lattice, cowplot, reshape2 and doParallel

