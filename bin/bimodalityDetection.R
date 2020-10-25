#*************************************************
#main ----
#
#   ATTENTION!!! This code is provided "AS-IS", 
#   with no warranties, express or implied, and 
#   hereby disclaims all implied warranties, 
#   including any warranty of merchantability and 
#   warranty of fitness for a particular purpose. 
#   It's released under GPL v2.
#*************************************************

#*************************************************
# Updated in 16/10/2020
#   Parallel mode
#   New graphics
#   Filter for low samples clusters
#   cluster classification by N
#*************************************************


#Clean all variables ----
rm(list=ls(all=TRUE))

#Are you having trouble using parallel processing?
#Set parallelMode = F
parallelMode = T

#File location ----
#Did you change it to your base location?
dirBase<-"Place here the correct name of your work folder"

#dirBase<-"/home/clovis/Doutorado/Projetos/Bimodalidade/TD0/"

#figures
dirFig<-file.path(dirBase,"figures/")
#bin dir
binDir<-file.path(dirBase,"bin/")


source(file.path(binDir,"allFunctions.R"))

loadDependencies()
unpack(dirBase)

#Threshold variables and parameters----
#only samples with expression above this value will be consider
minExpression = 0.02 
#Discard all genes with a sample count below this value 
minSampleSize = 50
#Ignore all cluster with size below this percentage of total samples
minClusterSize = 0.1
#minimum difference between peaks and valleys
limiarUp = 0.1 #threshold Up
#ignore peaks with density values below this value
limiarDw = 0.2 #threshold Down
#smoothing factor
atenuacao = 0.05
#Transform all expression values using log2
useLog = F

#Sample identification
#You can use more than one sample if necessary
#It will be processed in the appearance order
vtipo<-c("KIRKNew")

#Sample file name
#Must be in the same order above
vfileName<-c("KIRC_FPKM_TP.tsv")

i=1
#Main loop ----
#Loop to process all samples
for(i in 1:length(vfileName)){
  tipo<-vtipo[i]
  fileName<-vfileName[i]
  cat("Processing",tipo,fileName,"\n")

  dirFigAtu = paste0(dirFig,tipo,"/")
  
  if(parallelMode){
    processaPar(dirBase = dirBase, 
             dirFig = dirFigAtu, 
             atenuacao = atenuacao,
             minExpression = minExpression, 
             minSampleSize =minSampleSize,
             minClusterSize = minClusterSize,
             tipo = tipo,
             fileName=fileName,
             useLog = useLog)
    
  }else{
    processa(dirBase = dirBase, 
           dirFig = dirFigAtu, 
           atenuacao = atenuacao,
           minExpression = minExpression, 
           minSampleSize =minSampleSize,
           minClusterSize = minClusterSize,
           tipo = tipo,
           fileName=fileName,
           useLog = useLog)
  }
  
}

#dirFig = dirFigAtu
