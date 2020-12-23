#######################################################################
#################### SURVIVAL ANALYSIS #############################
#######################################################################

library(data.table)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(stringr)
library(clipr)

##########################################################
###   HNSC DATABASE ORGANIZATION AND FILTERING  ####
##########################################################
##########################################################
setwd("E:/ArtigoJosivan/samples/hnscNon")

## Creates a spreadsheet with gene data per peak

# lista todos os arquivos
files1.1 <- list.files(pattern = "*_0.02_peak1_cluster1.lst")
files1.2 <- list.files(pattern = "*_0.02_peak1_cluster2.lst")
files1.3 <- list.files(pattern = "*_0.02_peak1_cluster3.lst")

files1 <- c(files1.1,files1.2,files1.3)

files2.1 <- list.files(pattern = "*_0.02_peak2_cluster1.lst")
files2.2 <- list.files(pattern = "*_0.02_peak2_cluster2.lst")
files2.3 <- list.files(pattern = "*_0.02_peak2_cluster3.lst")

files2 <- c(files2.1,files2.2,files2.3)



lista1 <- data.frame() 
for(i in files1) 
{
  temp <- fread(i, header = F)
  temp$gene <- str_sub(i, end = -25)
  temp$pico <- "P1"
  temp$tumor <- "hnsc"
  lista1 <- rbind(lista1, temp) 
}

#  Do the same for list 2

lista2 <- data.frame()
for(i in files2)
{
  temp <- fread(i, header = F)
  temp$gene <- str_sub(i, end = -25)
  temp$pico <- "P2"
  temp$tumor <- "hnsc"
  lista2 <- rbind(lista2, temp)
}

### Organizes data for the merge with data clinical

dados <- rbind(lista1, lista2)
setnames(dados,c("CASE_ID","Gene","Pico","Tumor"))
dados$CASE_ID <- str_sub(dados$CASE_ID, end = 12)
dados$CASE_ID <- gsub("\\.", "\\-", dados$CASE_ID)
head(dados)


##### TCGA clinical data - cbioportal

clinical <- read.csv("hnsc_clinical_patient.txt", sep = '\t', header = F, stringsAsFactors = F)
clinical <- clinical[-c(1:4),]
colnames(clinical) <- clinical[1,]
clinical <- clinical[-1,]
clinical <- clinical %>% dplyr::select("PATIENT_ID","OS_STATUS","OS_MONTHS")
setnames(clinical,c("CASE_ID","OS_STATUS","OS_MONTHS"))
clinical$OS_MONTHS <- as.numeric(as.character(clinical$OS_MONTHS))

clinical <- clinical[(complete.cases(clinical)),]

# https://stat.ethz.ch/R-manual/R-devel/library/survival/html/Surv.html
clinical$OS_STATUS <- as.numeric(ifelse(clinical$OS_STATUS == "0:LIVING", 0, 1))
head(clinical)



### Function Survival

ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, cens.size = 1, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col,
                      size = cens.size)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ind[[i]]]),
        cens = c(0, s$n.censor[ind[[i]]]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if (length(surv.col == 1)) {
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if (surv.col[1] != 'gg.def') {
      pl + col
    } else {
      pl + scale_colour_discrete(name = gr.name)
    }
    
    line <- if (length(lty.est) == 1) {
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {
      scale_linetype_manual(name = gr.name, values = lty.est)
    }
    
    pl <- pl + line
    
    pl <- if (CI == T) {
      if (length(surv.col) > 1 && length(lty.est) > 1) {
        stop(
          'Either surv.col or lty.est should be of length 1 in order
          to plot 95% CI with multiple strata'
        )
      } else if ((length(surv.col) > 1 | surv.col == 'gg.def')[1]) {
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{
        pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low, lty = group), col = surv.col)
      }
    } else {
      pl
    }
    
    
    pl <- if (plot.cens == T & length(dat.cens) > 0) {
      pl + geom_point(data = dat.cens,
                      aes(y = surv, col = status),
                      shape = cens.shape,
                      col = cens.col,
                      size = cens.size)
    } else if (plot.cens == T & length(dat.cens) == 0) {
      stop ('There are no censored observations')
    } else
      (pl)
    
    pl <- if (back.white == T) {
      pl + theme_bw()
    } else
      (pl)
    pl
  }
  pl <- if (strata == 1) {
    ggsurv.s(
      s,
      CI ,
      plot.cens,
      surv.col ,
      cens.col,
      lty.est,
      lty.ci,
      cens.shape,
      back.white,
      xlab,
      ylab,
      main
    )
  } else {
    ggsurv.m(
      s,
      CI,
      plot.cens,
      surv.col ,
      cens.col,
      lty.est,
      lty.ci,
      cens.shape,
      back.white,
      xlab,
      ylab,
      main
    )
  }
  pl
}

setwd("E:/ArtigoJosivan/survival/hnscNon")


# Function to generate graphs and summary of survival tests

tumor="hnsc"
aux = dados[dados$Tumor == "hnsc",]
for(gene in unique(aux$Gene)){
  aux2 = dados[dados$Gene == gene,]
  colnames(aux2) <- c("CASE_ID", "Gene", "Pico","Tumor")
  #aux2$CASE_ID <- str_sub(aux2$CASE_ID, end = 12)
  BIMODAL <- merge(aux2, clinical, by = "CASE_ID")
  fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Pico, data = BIMODAL)
  dif <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ Pico, data = BIMODAL)
  print("pdf")
  pdf(paste0(tumor, "_", gene, ".pdf"), width = 10, height = 10)
  print(ggsurv(fit) + theme_minimal())
  dev.off()
  
  sink(paste0(tumor, "_", gene, ".txt"))
  print(dif)
  sink()
}
