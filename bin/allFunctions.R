#install required packages 
loadDependencies <- function(){
  if (!require("bnlearn")) {
    install.packages("bnlearn")
  }
  # if (!require("modes")) {
  #   install.packages("modes")
  # }
  if (!require("signal")) {
    install.packages("signal")
  }
  if (!require("ggplot2")) {
    install.packages("ggplot2")
  }
  if (!require("mclust")) {
    install.packages("mclust")
  }
  if (!require("gridExtra")) {
    install.packages("gridExtra")
  }
  if (!require("grid")) {
    install.packages("grid")
  }
  if (!require("lattice")) {
    install.packages("lattice")
  }
  if (!require("cowplot")) {
    install.packages("cowplot")
  }
  if (!require("reshape2")) {
    install.packages("reshape2")
  }
  
  if (!require("doParallel")) {
    install.packages("doParallel")
  }
  
}


#unpack data
unpack<- function(dirBase){
  
  exec<-paste0(dirBase,"bin/unpack.sh")
  command<-paste0(exec," ",dirBase)
  system(command)
  
}
#calculates the first derivative of the expression density
fderivada = function(densidade){
  derivada <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(derivada)<-c("x","y")
  for(i in seq(1,length(densidade$x)-1,1)){
    #janela <- densidade[(i-2):(i+2),]
    derivada[i,2] <- (densidade[i,2]-densidade[i+1,2])/
      (densidade[i,1]-densidade[i+1,1])
    derivada[i,1] <- densidade[i,1]
  }
  return(derivada)
}

#cutoff of peaks and valleys using provided threshold

aplicaTH= function (densidade, picos, limiarUp, limiarDw, arqLog, gene){
  
  #Cancel peaks and valleys lower than limiarDw
  #There must be N peaks and N + 1 valleys. Add at the beginning and end if necessary
  #As values are normalized between 0 and 1, place first value at -0.1 and last at 1
  #There must be alternation between peaks and valleys too. If not, generates an error
  #Test if the difference between peaks and valleys exceedes the value provided in limiarUp
  #otherwise, cut peak and posterior valley
  
  
  #cria variavel de retorno
  resultado<-data.frame(matrix(ncol = 3, nrow = 0))
  #df temporario juntando tabela original e picos
  tmp <- merge(densidade,picos, by= 1)
  colnames(tmp)<- c("x","y","tipo")
  #elimina valores menores que limiarDw do sinal
  tmp <- tmp[tmp$y > limiarDw,]
  #erro se tmp ficar vazio
  if(nrow(tmp)<=1){
    cat(paste("Nenhum pico e vale a ser processado... ",gene,"\n"),file=arqLog, append = T)
    #print("Erro: ")
    colnames(resultado)<-c("x","y","tipo")
    return(resultado)
  }
  
  #acrescenta o primeiro vale, caso este não exista
  if(tmp$tipo[1] == 1){
    tmp<-rbind(c(-0.1,0,-1),tmp)
  }
  #acrescenta o ultimo vale, caso este não exista
  if(tmp$tipo[nrow(tmp)] == 1){
    tmp<-rbind(tmp,c(1,0,-1))
  }
  #verifica se existe alternancia entre picos e vales
  for(i in seq(1,nrow(tmp)-1,2)){
    #primeiro vale depois pico
    if(tmp$tipo[i]== 1 | tmp$tipo[i+1] == -1 ){
      cat(paste("Inconsistência no número de picos e vales... ",gene,"\n"),file=arqLog, append = T)
      #      print("Erro: ")
      colnames(resultado)<-c("x","y","tipo")
      return(resultado)
    }
  }
  #retorna sempre o primerio vale
  resultado <- data.frame(tmp[1,])
  
  #realiza verificação onde vale/pico/vales > limiarUp
  while(nrow(tmp) > 1){
    valeAnt<-tmp$y[1]
    pico   <-tmp$y[2]
    valePos<-tmp$y[3]
    if(limiarUp < pico - valeAnt & limiarUp < pico - valePos){
      resultado <- rbind(resultado, tmp[2,],tmp[3,])
      tmp <- tmp[-c(1,2),]
    }else{
      tmp <- tmp[-c(2,3),]
    }
  }
  colnames(resultado)<-c("x","y","tipo")
  return(resultado)
}

#locates all peaks
achaPico = function (derivada){
  #compara os sinais da derivada. Havendo diferença houve inflexão no gráfico
  sinal<- sign(derivada$y)
  #variavel de resultado
  resultado<-data.frame(matrix(ncol = 2, nrow = 0))
  tmp<-data.frame(matrix(ncol = 2, nrow = 1))
  for(i in seq(1,length(sinal)-1,1)){
    if(sinal[i] != sinal[i+1]){
      #em tipo picos recebem 1, vales recebem -1
      if(sinal[i] > sinal[i+1]){
        tmp[1,2]<-1
      }else{
        tmp[1,2]<- -1
      }
      tmp[1,1]<-derivada$x[i]
      
      resultado<- rbind(resultado,tmp)
    }
  }
  colnames(resultado)<-c("x","tipo")
  return(resultado)
}

normaliza <- function(x){
  M<-max(x)
  m<-min(x)
  return(c((x-m)/(M-m)) )
}


#Uses GMM to cluster samples groups
intervalo <- function(densidade, amostras, 
                      maxExpr, minExpr, 
                      nomeGene, dirFig, minExpression, 
                      atenuacao, tipo, qtdPicos){
  
  #parameter matching:
  #densidade,  
  #amostras, maxExpr,minExpr,
  #gene,dirFig,minExpression, 
  #atenuacao, tipo, qtdPicos
  
  #variavel de resultado onde cada elemento da lista conterá um vetor com o nome das amostras
  result <- list()
  library("mclust")
  #clusteriza
  mCluster<-Mclust(data = amostras$V3,G = qtdPicos+1)
  #vetor de break points
  breaks<-(seq(0,maxExpr,(maxExpr-minExpr)/nrow(densidade)))
  breaks[length(breaks)]<-maxExpr
  
  valHist <- data.frame(x =      c(rep(NA, length(breaks)-1)), 
                        dens =   rep(0, length(breaks)-1),
                        counts = rep(0, length(breaks)-1),
                        cor =    rep(qtdPicos+1, length(breaks)-1),
                        stringsAsFactors = F)
  valHist$x<-breaks[1:length(breaks)-1]
  #df de médias e sd
  mdsd <- data.frame(cluster = c(rep(NA, qtdPicos)), 
                     mean = rep(0, qtdPicos),
                     sd   = rep(0, qtdPicos),
                     stringsAsFactors = F)
  
  #número de clusters a serem encontrados será qtdPicos +1
  
  #cores<-seq(1:qtdPicos)
  cores<-c("Distr 1"="#FF8C69","Distr 2"="#93CCEA","No Distr"="#C4C3D0")
  
  hist2<-list()
  i=1
  for(i in seq(1:qtdPicos)){
    #separa dados dos clusters
    dt <- data.frame(x = mCluster$data[mCluster$classification == i & mCluster$uncertainty < 0.45])
    if(nrow(dt)==0){
      result<-0
      return(result)
    }
    #dt<-dt[order(dt$x),]
    #calcula média e sd do subset
    mdsd$cluster[i] <- i
    mdsd$mean[i] <- mean(dt$x)
    mdsd$sd[i] <- sd(dt$x)
    #calcula histograma
    den<-round(nrow(dt)/20)
    #não deixa dar 0 no denominador
    den<-ifelse(den,den,1)
    
    histograma <- hist(dt$x, 
                       breaks = (round(nrow(dt)/den,0)+1) ,
                       plot = F)
    hist2[[i]]<-data.frame(x=histograma$breaks[2:length(histograma$breaks)], y=histograma$counts)
    histograma <- hist(dt$x, breaks = breaks ,plot = F)
    tmp<-data.frame(x=histograma$breaks[1:length(histograma$breaks)-1], 
                    dens=histograma$density,
                    counts=histograma$counts,
                    cor=i,
                    stringsAsFactors = F)
    tmp<-merge(valHist,tmp,by = "x" )
    valHist$x<-tmp$x
    valHist$dens[tmp$dens.y != 0]<- tmp$dens.y[tmp$dens.y != 0]
    valHist$cor[tmp$dens.y != 0]<- tmp$cor.y[tmp$dens.y != 0]
    valHist$counts<- tmp$counts.y+tmp$counts.x
  }
  
  #integra os demais dados
  dt <- data.frame(x = mCluster$data)
  #calcula histograma
  histograma <- hist(dt$x, breaks = breaks ,plot = F)
  tmp<-data.frame(x=histograma$breaks[1:length(histograma$breaks)-1], 
                  dens=histograma$density,
                  counts=histograma$counts,
                  #                        cor=cores[i],
                  cor=qtdPicos+1,
                  stringsAsFactors = F)
  
  #copia dados ainda não computados
  valHist$dens[valHist$cor == qtdPicos+1]<- tmp$dens[valHist$cor == qtdPicos+1]
  valHist$counts[valHist$cor == qtdPicos+1]<- tmp$counts[valHist$cor == qtdPicos+1]
  
  
  
  
  #Normaliza os valores de expressão apenas para fins do gráfico
  amostras$norm <- normaliza(amostras$V3)
  
  if(minExpression == 0){
    xtext<-"Expression Values - all values"
  }else{
    xtext<-paste0("Expression Values - over ", minExpression,"% of maximum expression")
  }
  
  xtext<-"Expression Values"
  library(ggplot2)
  subG<-list()
  g<-ggplot()+
    ggtitle(paste(tipo," ",nomeGene," Threshold Y = ",limiarUp," Threshold X = ",minExpression))
  i=2
  for(i in 1:qtdPicos){
    if(nrow(valHist[valHist$cor == i,])==0){
      result<-0
      return(result)
    }
    
    media <- mdsd$mean[i]
    sd <- mdsd$sd[i]
    
    x=c(media-sd,media+sd,media+sd,media-sd)
    y=c(max(densidade$yDens),max(densidade$yDens),min(densidade$yDens),min(densidade$yDens))
    area<-data.frame(x=x,y=y)
    area$x[area$x<0]=0
    
    g<-g+geom_polygon(data = area, aes(x,y),col='Gray91', fill='Gray91')
    g<-g+geom_vline(xintercept = media,lty = 2, col="Gray51")
    g<-g+geom_point(data= valHist[valHist$cor == i,], aes(x,counts*max(densidade$yDens)/ max(valHist$counts),shape="ponto"),col=i+1)
    
    #separa as amostras resultado
    result[[i]] <- data.frame(sample=amostras$V1[mCluster$classification == i & mCluster$uncertainty < 0.45],
                              expression=amostras$V3[mCluster$classification == i & mCluster$uncertainty < 0.45],
                              stringsAsFactors = F)    
    
    subG[[i]]<-ggplot()+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90,
                                       size = 8))+
      xlab(paste0(xtext, " - Group ",i," (bin = ", hist2[[i]][2,1]-hist2[[i]][1,1],")"))+
      ylab("Density") +
      #xlim(c(0,max(densidade$xDens)))+
      scale_y_continuous(name = "Counts")+
      geom_col(data= hist2[[i]], aes(x=as.factor(x),y=y),
               colour="gray", fill=i+1)
    
  }
  if(nrow(valHist[valHist$cor == qtdPicos+1 & valHist$counts > 1 ,])>0){
    g<-g+geom_point(data= valHist[valHist$cor == qtdPicos+1 & 
                                    valHist$counts > 1 ,], 
                    aes(x,counts*max(densidade$yDens)/ max(valHist$counts),
                        shape="ponto"), 
                    col = 4)
  }
  if(nrow(valHist[valHist$cor == qtdPicos+1 & valHist$counts == 1 ,])>0){
    g<-g+geom_point(data= valHist[valHist$cor == qtdPicos+1 & 
                                    valHist$counts == 1 ,], 
                    aes(x,counts*max(densidade$yDens)/ max(valHist$counts),
                        shape="ponto"),
                    col = 4,
                    alpha=1/2)
  }
  g<-g+
    theme_bw()+
    xlab(xtext)+
    ylab("Density") +
    scale_y_continuous(name = expression("Density"), 
                       #limits = c(0, max(densidade$yDens)),
                       sec.axis = sec_axis(~ ./max(densidade$yDens)* max(histograma$counts) , 
                                           name = "Counts"))+
    geom_line(data=densidade, aes(xDens,yDens,linetype="before"),col="blue")
  
  
  g<-g+scale_shape_manual("",guide = F,values = c("ponto"=21,"none"=NA))+
    scale_linetype_manual("Density",guide = F,values = c("before"=1,"after"=2),
                          labels = c("Before","After"))
  
  result[[qtdPicos+1]] <- data.frame(sample=amostras$V1[mCluster$classification == qtdPicos+1 |
                                                          (mCluster$classification != qtdPicos+1 & mCluster$uncertainty >= 0.45)],
                                     expression=amostras$V3[mCluster$classification == qtdPicos+1 |
                                                              (mCluster$classification != qtdPicos+1 & mCluster$uncertainty >= 0.45)],
                                     stringsAsFactors = F)    
  return(list(result,g,subG,max(densidade$yDens)))
}

#Print and show the detection chart
printGraf<-function(dirFig,
                    g,
                    subG,
                    gene,
                    minExpression){

  library("gridExtra")
  library("grid")
  
  #cria o grid do gráfico
  layout<-rbind(c(1,1),
                c(2,3))
  
  g1<-grid.arrange(g,subG[[1]],subG[[2]],
                   layout_matrix=layout)
  
  ggsave(filename = paste0(dirFig,"/",gene,"_",minExpression,".pdf"),
         plot = g1,
         device = "pdf",
         width = 11,height = 8)
  
}


#Data processing and peak detection
processaPicos<- function(amostras = amostras,
                         arqLog = arqLog,
                         gene){
  
  #calcula densidade ----
  g <- density(amostras$V3)
  #g <- density(amostras$V3,bw="SJ") #outro tipo de cálculo de densidade
  densidade<-data.frame(cbind(g$x,g$y))
  colnames(densidade) <- c("xDens","yDens")
  
  #Valores maximos e mínimos de expressão
  medExpe <- mean(amostras$V3)
  sdExpe <- sd(amostras$V3)
  maxExpr <- max(amostras$V3)
  minExpr <- min(amostras$V3)
  
  #Normaliza x da densidade entre 0 e 1
  densidade$xNorm<- (densidade$xDens)/(max(densidade$xDens))
  maxExprN<-(maxExpr)/(max(densidade$xDens))
  minExprN<-(minExpr)/(max(densidade$xDens))
  medDens <- (medExpe-min(densidade$xDens))/(max(densidade$xDens)-min(densidade$xDens))
  sdDens <- (sdExpe-min(densidade$xDens))/(max(densidade$xDens)-min(densidade$xDens))
  
  #normaliza y
  M<-max(densidade$yDens)
  m<-min(densidade$yDens)
  densidade$yDNorm<-sapply(densidade$yDens,function(x){
    tmp<-c((x-m)/(M-m))  
  })
  
  #Calcula derivadas
  derivadaV <- fderivada(densidade[c("xNorm","yDNorm")])
  #filtro passa baixa
  derivadaS<-smooth.spline(derivadaV, spar = atenuacao, df = 5)
  derivadaS<- data.frame(cbind(derivadaS$x,derivadaS$y))
  colnames(derivadaS)<- c("x","y")
  #sem filtro
  derivadaS<-derivadaV
  
  
  #Normaliza derivadas
  M<-max(derivadaS$y)
  m<-min(derivadaS$y)
  derivadaS$y<-sapply(derivadaS$y,function(x){
    tmp<-c((x)/(M))
  })
  
  M<-max(derivadaV$y)
  derivadaV$y<-sapply(derivadaV$y,function(x){
    tmp<-c((x)/(M))
  })
  
  picos<-achaPico(derivadaS)
  valPicos<-aplicaTH(densidade=densidade[c("xNorm","yDNorm")], 
                     picos=picos, limiarUp = limiarUp, 
                     limiarDw = limiarDw, arqLog = arqLog,gene = gene )
  d<-densidade
  densidade$picos <- merge(densidade[c("xNorm")], valPicos, by = 1, all.x = T)[,3]
  
  qtdPicos <-nrow(valPicos[valPicos$tipo==1,])
  
  return(list(densidade,
              qtdPicos,
              maxExpr,
              minExpr,
              medExpe,
              sdExpe ))
  
}

#Read samples file and start processing
processa = function(dirBase, 
                    dirFig, 
                    atenuacao,
                    minExpression,
                    tipo,
                    fileName){
  #Create folder for figures
  dirFigFake<-paste0(dirFig,"fake/")
  if (!dir.exists(dirFig)){
    dir.create(dirFig)
    dir.create(dirFigFake)
  }
  #Create folder for samples
  dirSamples<-paste0(dirBase,"/samples/",tipo,"/")
  dirSamplesFake<-paste0(dirSamples,"fake/")
  if (!dir.exists(dirSamples)){
    dir.create(dirSamples)
    dir.create(dirSamplesFake)
  }
  
  
  #realiza o processamento propriamente dito
  
  #Open the log file
  arqLog<-paste0(dirBase,
                 "log/",
                 tipo,
                 format(Sys.time(), "%X_%Y_%m_%d"),
                 "log.txt")
  
  cont=0
  lista <- data.frame(matrix(ncol = 2, nrow = 0))
  erros <- data.frame(matrix(ncol = 1, nrow = 0))
  dfTmp <- data.frame(matrix(ncol = 2, nrow = 0))
  
  dirDados<-paste0(dirBase,"data/")
  
  
  #Read the data file
  setwd(dirDados)
  
  cat("Reading File",fileName,"\n")
  tryCatch(allAmostras<-read.table(paste0(dirDados,
                                          fileName), 
                                   sep = "\t",
                                   header = T,
                                   stringsAsFactors = F),
           error = function(e) {print(paste("Error opening file ",fileName))
             writeLines(paste("Error opening file ",fileName))
             cat(paste("Error opening file ",fileName,"\n"),file=arqLog, append = T)
             return(1)})
  cat("Filtering file",fileName,"\n")
  
  allAmostras<-allAmostras[!allAmostras$symbol == '-',]
  genes <- unique(allAmostras$symbol)
  genes<-genes[order(genes)]
  #genes="AP5Z1"
  
  #Prossesing each genes in sample
  for(gene in genes){
    #le amostras ----
    amostras <- na.omit(allAmostras[allAmostras$symbol == gene,
                                    c("sample","symbol","FPKM")])
    colnames(amostras)<-c("V1","V2","V3")
    if(nrow(amostras) == 0){
      writeLines(paste("Erro processing ",gene,": no expressions values found."))
      cat(paste("Erro processing",gene,":  no expressions values found.","\n"),file=arqLog,append = T)
      next()}
    #nome do gene ----
    #genes com alias tem nomes separados por "//" e devem ser eleiminados
    if(grepl("/",gene)){
      gene<-strsplit(gene,"/")[[1]][1]
    }
    
    writeLines(paste("Processing ", gene))
    cat(paste("Processing ", gene,"\n"),file=arqLog,append = T)
    # print(nrow(amostras))
    
    
    #realiza o corte de tudo que ficar abaixo do percentual da expressão máxima informado
    if(minExpression != 0){
      #amostras<-amostras[amostras$V3>=maxExpr*minExpression/100,]
      amAbaixo<-amostras[amostras$V3<minExpression,] #usa valor fixo
      amostras<-amostras[amostras$V3>=minExpression,] #usa valor fixo
    }
    if(nrow(amostras)<50){
      writeLines(paste("Not enough data - just ", nrow(amostras), " samples for ", gene))
      cat(paste("Not enough data - just ", nrow(amostras), " samples for ", gene,"\n"),file=arqLog,append = T)
      next()
    }
    pPicos<-processaPicos(amostras,arqLog,gene)
    densidade<-pPicos[[1]]
    qtdPicos<-pPicos[[2]]
    maxExpr<-pPicos[[3]]
    minExpr<-pPicos[[4]]
    medExpe<-pPicos[[5]]
    sdExpe<-pPicos[[6]]
    nomeGene = gene
    if(qtdPicos >= 2){
      # relevant results ----
      resultado <- intervalo(densidade = densidade,  
                             amostras = amostras, 
                             maxExpr = maxExpr,
                             minExpr = minExpr,
                             nomeGene = gene,
                             dirFig = dirFig,
                             minExpression = minExpression, 
                             atenuacao = atenuacao, 
                             tipo = tipo, 
                             qtdPicos = qtdPicos)
      
      if(class(resultado) == "numeric"){
        writeLines(paste("Error processing ",gene,": GMM did not find two peaks."))
        cat(paste("Erro processing",gene,": GMM did not find two peaks.","\n"),file=arqLog,append = T)
        next()}
      #extrai conteudo do retorno da funcao
      g<-resultado[[2]]
      subG<-resultado[[3]]
      maxYDens<-resultado[[4]]
      resultado<-resultado[[1]]
      
      #confere se há mesmo dois picos 
      P1e2<-c(resultado[[1]]$sample,
              resultado[[2]]$sample)
      #filtra amostras dos picos
      amostrasPicos<-amostras[amostras$V1 %in% P1e2,]
      #teste picos novamente 
      pPicos2<-processaPicos(amostrasPicos,arqLog,gene)
      densidade2<-pPicos2[[1]]
      qtdPicos2<-pPicos2[[2]]
      #adiciona linha dos picos sem amostras escluidas pela GMM
      maxD2<-maxYDens/max(densidade2$yDens)
      g<-g+geom_line(data=densidade2,
                     aes(xDens,yDens*maxD2,linetype="after"),
                     col="orange")
      if(qtdPicos==qtdPicos2){
        dirFigDest<-dirFig
        dirSamplesDest<-dirSamples
      }else{
        dirFigDest<-dirFigFake
        dirSamplesDest<-dirSamplesFake
      }
      
      printGraf(dirFig = dirFigDest,
                g = g,
                subG = subG,
                gene = gene,
                minExpression = minExpression)
      
      #salva identificador das amostras
      for(i in 1:length(resultado)){
        write.table(as.character(resultado[[i]]$sample),
                    paste0(dirSamplesDest,"/",gene,"_",minExpression,"_",i,".lst"),
                    row.names = F, 
                    col.names = F)
        #totaliza as contagens
        
      }
      sumario<-data.frame(th=nrow(amAbaixo),
                          low=nrow(resultado[[3]][resultado[[3]]$expression<min(resultado[[1]]$expression),]),
                          p1=nrow(resultado[[1]]),
                          btw=nrow(resultado[[3]][(resultado[[3]]$expression>max(resultado[[1]]$expression)&
                                                     resultado[[3]]$expression<min(resultado[[2]]$expression)),]),
                          p2=nrow(resultado[[2]]),
                          high=nrow(resultado[[3]][resultado[[3]]$expression>max(resultado[[2]]$expression),]))
      sumario$tot<-sum(sumario)
      write.table(sumario,paste0(dirSamplesDest,"/",gene,"_",minExpression,".summary"),
                  row.names = F, 
                  col.names = T,
                  sep = "\t" )
      
      cont <- cont+1
      dfTmp[1,1]<-cont
      dfTmp[1,2]<-as.character(gene)
      
      lista <- rbind(lista,dfTmp)
      
      cat(c("*******************************************************\n",
            paste("Bimodality for ", gene,"\n"),
            "*******************************************************\n"),
          file = arqLog, append = T)
      writeLines(c("*******************************************************",
                   paste("Bimodality for ", gene),
                   "*******************************************************"))

    }
    
    
  }
  
  
  #dev.off()
  
  colnames(lista)<-c("Nr","Gene")
  colnames(erros)<-c("Arquivo")
  
  #write.csv(lista,paste0(dirFig,"lista.csv"),row.names = F)
  #write.csv(erros,paste0(dirFig,"erros.csv"),row.names = F)
  #close(arqLog)
  
}

