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


#Clean all variables ----
rm(list=ls(all=TRUE))

#File location ----
#Did you change it to your base location?
dirBase<-"Place here the correct name of your work folder"

#figures
dirFig<-paste0(dirBase,"figures/")
#bin dir
binDir<-paste0(dirBase,"bin/")


source(paste0(binDir,"allFunctions.R"))

loadDependencies()
unpack(dirBase)

#Threshold variables ----
#minimum difference between peaks and valleys
limiarUp = 0.05 #threshold Up
#ignore peaks with density values below this value
limiarDw = 0.1 #threshold Down
#smoothing factor
atenuacao = 0.05
#only samples with expression above this value will be consider
minExpression = 0.02 

#Sample identification
#You can use more than one sample if necessary
#It will be processed in the appearance order
vtipo<-c("KIRC")

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
  
  
  processa(dirBase = dirBase, 
           dirFig = dirFigAtu, 
           atenuacao = atenuacao,
           minExpression = minExpression, 
           tipo = tipo,
           fileName=fileName)
  
}

