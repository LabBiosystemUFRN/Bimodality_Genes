# <center> Bimodal Gene Expression Detection</center>
## <center> A new method to identify genes with bimodal gene expression</center>

####  <span style="color:red">The code located here is provided "AS-IS", with no warranties, express or implied, and hereby disclaims all implied warranties, including any warranty of merchantability and warranty of fitness for a particular purpose. It's released under GPL v2.</span> </center>
<br>

This repository contains the files necessary to reproduce the results reported in our bimodality analysis in gene expression. The bin folder contains the code used to generate analysis and figures. It's based on R software (version 3.6.3 2020-02-29) and uses the libraries' bnlearn, signal, ggplot2, mclust, gridExtra, grid,lattice, cowplot, reshape2 and doParallel.

The *data_files*  folder contains data from 522 patients with a primary tumor of clear cell renal carcinoma (KIRC).

The original data is splited in two files: KIRC_FPKM_TP.tsv.tar.bz.part.aa and KIRC_FPKM_TP.tsv.tar.bz.part.ab

The script unpack.sh reassemble and extract the data. 


If your need to process another dataset without compaction, first comment the line

> \# unpack(dirBase)

Your file must be placed inside the folder named data, and the file must have the following structure, separated by tabs:

* sample - with the sample identifier. The final output files will use this identifier to list what samples belong to each cluster;

* symbol - with the gene identifier, usually the gene symbol gene. It will be used to report bimodal genes;

* FPKM - with the values ​​of expression.

Below you can see an example of the expression file format accepted by the algorithm:

"sample"	|	"symbol"	|	"FPKM"
--------	|	--------	|	--------
"TCGA-66-2783"	|	"LINC02082"	|	0
"TCGA-66-2783"	|	"RAB4B"	|	2.44648524683
"TCGA-66-2783"	|	"TIGAR"	|	3.01211011845
"TCGA-66-2783"	|	"RNF44"	|	24.7499170861
"TCGA-66-2783"	|	"DNAH3"	|	0.809002751355
"TCGA-66-2783"	|	"RPL23A"	|	80.0582471777
"TCGA-66-2783"	|	"ARL8B"	|	44.5540484135
"TCGA-66-2783"	|	"CALB2"	|	10.2637161003
"TCGA-66-2783"	|	"MFSD3"	|	4.88189970348
"TCGA-66-2783"	|	"LINC00636"	|	0
"TCGA-66-2783"	|	"PIGV"	|	8.14476960589




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

By default, the script will run in parallel mode. If you have memory issues using this mode, please, change the following line in the script:

> parallelMode = F

Then you can run :
> Rscript Bimodality_Genes/bin/bimodalityDetection.R

### R dependencies

> bnlearn, signal, ggplot2, mclust, gridExtra, grid,lattice, cowplot, reshape2 and doParallel

