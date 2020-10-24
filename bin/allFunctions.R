#install required packages #
loadDependencies <- function(){
  if (!require("bnlearn")) {
    install.packages("bnlearn", repos = "http://cran.us.r-project.org")
  }
  if (!require("signal")) {
    install.packages("signal", repos = "http://cran.us.r-project.org")
  }
  if (!require("ggplot2")) {
    install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  }
  if (!require("mclust")) {
    install.packages("mclust", repos = "http://cran.us.r-project.org")
  }
  if (!require("gridExtra")) {
    install.packages("gridExtra", repos = "http://cran.us.r-project.org")
  }
  if (!require("grid")) {
    install.packages("grid", repos = "http://cran.us.r-project.org")
  }
  if (!require("lattice")) {
    install.packages("lattice", repos = "http://cran.us.r-project.org")
  }
  if (!require("cowplot")) {
    install.packages("cowplot", repos = "http://cran.us.r-project.org")
  }
  if (!require("reshape2")) {
    install.packages("reshape2", repos = "http://cran.us.r-project.org")
  }
  
  if (!require("doParallel")) {
    install.packages("doParallel", repos = "http://cran.us.r-project.org")
  }
  
}


#unpack data
unpack<- function(dirBase){
  
  cat("Unpacking: ")
  
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
  
  #acrescenta o primeiro vale, caso este não exista
  if(tmp$tipo[1] == 1){
    tmp<-rbind(c(-0.1,0,-1),tmp)
  }
  #acrescenta o ultimo vale, caso este não exista
  if(tmp$tipo[nrow(tmp)] == 1){
    tmp<-rbind(tmp,c(1,0,-1))
  }
  #verifica se existe alternancia entre picos e vales
  i=5
  for(i in seq(1,nrow(tmp)-1,2)){
    #primeiro vale depois pico
    if(tmp$tipo[i]== 1 | tmp$tipo[i+1] == -1 ){
      cat(paste("Number of peaks and valleys doesn't match... ",gene,"\n"),file=arqLog, append = T)
      #      print("Erro: ")
      colnames(resultado)<-c("x","y","tipo")
     return(resultado)
    }
  }
  
  #elimina picos menores que limiarDw do sinal
  # eliminando o vale subsequente tb
  tmp2<-tmp[1,]
  i=2
  #testa somente picos
  for(i in seq(from=2,to=nrow(tmp),by=2)){
    if(tmp$y[i] > limiarDw){
      tmp2<-rbind(tmp2,tmp[i:(i+1),])  
    }    
  }
  #tmp <- tmp[tmp$y > limiarDw,]
  #erro se tmp ficar vazio
  if(nrow(tmp)<=1){
    cat(paste("No peak and valley to be processed... ",gene,"\n"),file=arqLog, append = T)
    #print("Erro: ")
    colnames(resultado)<-c("x","y","tipo")
    return(resultado)
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
achaPico = function (derivadaS){
  #compara os sinais da derivada. Havendo diferença houve inflexão no gráfico
  sinal<- sign(derivadaS$y)
  sinal<-sinal[sinal!=0]
  #variavel de resultado
  resultado<-data.frame(matrix(ncol = 2, nrow = 0))
  tmp<-data.frame(X1=0.0,X2=0.0)
  i=142
  for(i in seq(1,length(sinal)-1,1)){
    if(sinal[i] != sinal[i+1]){
      #em tipo picos recebem 1, vales recebem -1
      if(sinal[i] > sinal[i+1]){
        tmp[1,2]<-1
      }else{
        tmp[1,2]<- -1
      }
      tmp[1,1]<-derivadaS$x[i]
      
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
                      nomeGene, dirFig, minExpression, minClusterSize,
                      atenuacao, tipo, qtdPicos,coordPicos){
  library(ggplot2)
  
  #parameter matching:
  #densidade,  
  #amostras, maxExpr,minExpr,
  #gene,dirFig,minExpression, 
  #atenuacao, tipo, qtdPicos
  
  #variavel de resultado onde cada elemento da lista conterá um vetor com o nome das amostras
  result <- list()
  #g contém todos os 3 graficos de fases de detecção gerados
  g<-list()
  
  #vetor de break points
  breaks<-(seq(0,maxExpr,(maxExpr-minExpr)/nrow(densidade)))
  breaks[length(breaks)]<-maxExpr
  
  
  #gráfico 1 ----
  #contém a densidade e os bins de todas as amostras
  
  #calcula histograma completo
  tmp <- hist(amostras$V3, 
                     breaks = breaks ,
                     plot = F)
  histTodos<-data.frame(breaks=breaks[2:length(breaks)], 
                       counts=tmp$counts,
                       stringsAsFactors = F)
  histTodos<-histTodos[histTodos$counts!=0,]
  
  xtext<-"Expression Values"
  
  g1<-ggplot()+
    ggtitle(paste0("Gene ",nomeGene," - Peak detection"))+
    theme_bw()+
    xlab(xtext)+
    ylab("Density") +
    scale_y_continuous(name = expression("Density"), 
                       #limits = c(0, max(densidade$yDens)),
                       sec.axis = sec_axis(~ ./max(densidade$yDens)* max(histTodos$counts) , 
                                           name = "Counts"))+
    geom_line(data=densidade, aes(xDens,
                                  yDens,
                                  linetype="Density",
                                  color="Density"))+
    geom_point(data=histTodos, aes(x=breaks, 
                                   y=(counts*max(densidade$yDens)/ max(counts)),
                                   color= "Samples count",
                                   shape="Samples count"))+
    geom_point(data = coordPicos,aes(x,
                                     y,
                                     shape="Peaks", 
                                     color= "Peaks",
                                     fill = "Peaks"))+
    scale_fill_manual("",guide=F,values = c("Peaks"="magenta2"))+
    scale_shape_manual("",values = c("Samples count"=21,
                                     "Peaks"=25))+
    scale_color_manual("",guide=F,values = c("Samples count"=alpha("black",0.5),
                                     "Peaks"="magenta2", 
                                     "Density"= "orange"))+
    scale_linetype_manual("",values = c("Density"=1))+
    guides(shape = guide_legend(override.aes = list(shape = c(25,21), 
                                                   fill=c("magenta2",NA), 
                                                   color= c("magenta2", "black"))),
           linetype = guide_legend(override.aes = list(color = "orange")))
  
  g1
  
  
  library("mclust")
  #clusterização ----
  mCluster<-Mclust(data = amostras$V3,G = qtdPicos+1,modelNames = "V")
  
  if(is.null(mCluster$data)){
    result<-1
    return(result)
  }

  clusteres<-data.frame(amostras$V1,
                       mCluster$data,
                       mCluster$classification,
                       mCluster$uncertainty)  
  colnames(clusteres)<-c("sample","expression","cluster","uncertainty")
  clusteres<-clusteres[order(clusteres$expression),]
  
  #separa e ordena clusteres pelo tamanho
  #cria df com o valor dos histogramas por cluster
  valHist <- data.frame(x =      numeric(), 
                        counts = numeric(),
                        cor =    character(),
                        stringsAsFactors = F)

  #quantidades por clusters
  #teste de continuidade
  #verifica se cada cluster é uma unidade única e contínua
  #caso não seja, o segmento de menor quantidade deixa de ser computado
  
  fragm<-data.frame(cluster = character(),
                    inicio = numeric(),
                    fim = numeric(),
                    total = numeric(),
                    stringsAsFactors = F)
  clOld="0"
  i=1
  for (i in 1:nrow(clusteres)) {
    clNew<-clusteres$cluster[i]
    if(clNew != clOld){
      if(clOld !=0){
        fragm$fim[nrow(fragm)]<- (i-1)
      }
      fragm[(nrow(fragm)+1),]<-NA
      fragm$cluster[nrow(fragm)]<-clNew
      fragm$inicio[nrow(fragm)]<-i
      clOld<-clNew
    }
  }
  fragm$fim[nrow(fragm)]<- i
  fragm$total<-fragm$fim - fragm$inicio + 1
  
  #controle de continuidade
  continuid<- as.data.frame(t(table(fragm$cluster)),stringsAsFactors = F)
  continuid$Var1<-NULL
  continuid<-continuid[continuid$Freq >1,]
  if(nrow(continuid)>0){
    i="4"
    for(i in unique(continuid$Var2)){
      tmp<-fragm[fragm$cluster == i,]
      #Apenas a ultima linha contendo o maior numero de amostras
      lastLin<-which(tmp$total == max(tmp$total))
      lastLin<-lastLin[length(lastLin)]
      tmp<-tmp[-lastLin,]
      #transfere para cluster picos +2 (inutil) todos os descontínuos
      for (j  in 1:nrow(tmp)) {
        clusteres$cluster[tmp$inicio[j]:tmp$fim[j]]<-qtdPicos+2
        
      }
    }
  }

  #transfere amostras com baixa confiabilidade para cluster picos +2 (inutil)
  clusteres$cluster[clusteres$uncertainty >= 0.45]<- qtdPicos+2
  
  #regitra contagem de amostras por cluster
  tmp<-as.data.frame(t(table(clusteres$cluster)),stringsAsFactors = F)
  #renumera clusteres caso algum tenha sido apagado, exceto qtdPicos+2
  maxCl<-nrow(tmp[tmp$Var2 != qtdPicos+2,])
  tmp$Var2[tmp$Var2 != qtdPicos+2]<-c(1:maxCl)
  #df de ordenamento
  ordem<-data.frame(prev=c(tmp$Var2),
                    new = NA,
                    count=NA)
  
  #renumera clusteres para 90 e alguma coisa para depois 
  #garantir que o cluster 1 é sempre o maior
  # e o 2 o segundo maior
  ordem$prev<-as.numeric(tmp$Var2)
  ordem$count<-as.numeric(tmp$Freq)
  ordem<-ordem[ordem$prev != qtdPicos+2,]
  minOrdem<-min(ordem$prev)
  maxOrdem<-max(ordem$prev)
  ordem$prev<-ordem$prev+90
  ordem<-ordem[order(ordem$count,decreasing = T),]
  ordem$new<-c(minOrdem:maxOrdem)

  #Garante que cada pico possui pelo menos 10% das amostras totais
  threshosdAmostras<- ceiling(nrow(clusteres)*minClusterSize)
  i=3
  for(i in 1:qtdPicos){
    #cat(ordem$count[i])
    if(ordem$count[i] < threshosdAmostras){
      result<-0
      return(result)
    }
  }

  #ordem Clusteres ----
  #trocando a ordem dos clusteres
  clusteres$cluster[clusteres$cluster <= qtdPicos +1]<-clusteres$cluster[clusteres$cluster <= qtdPicos +1]+90
  i=1
  for(i in 1:nrow(ordem)) {
    clusteres$cluster[clusteres$cluster == ordem$prev[i]]<-ordem$new[i]
    
  }

  #totalização ----
  #totaliza resultados e histogramas
  i=4
  for(i in 1:(qtdPicos+2)){
    result[[i]]<- clusteres[clusteres$cluster == i,]    
    #calcula histograma do cluster
    tmp <- hist(result[[i]]$expression, 
                breaks = breaks ,
                plot = F)
    tmp<-data.frame(x=breaks[2:length(breaks)], 
                    counts=tmp$counts,
                    stringsAsFactors = F)
    tmp<-tmp[tmp$counts!=0,]
    if(nrow(tmp)==0){
      next()
    }
    tmp$cor<-as.character(i)
    valHist<-rbind(valHist,tmp)
  }  
  
  #grafico 2 ----
  #define cores e labels
  cores<-c("1"=alpha("red",0.75),
           "2"=alpha("blue",0.75),
           "3"= alpha("green4",0.75),
           "4"= alpha("green3",0.75),
           "5"= alpha("lightblue4",0.75),
           "6"= alpha("lightsalmon2",0.75),
           "7"=alpha("gray50",1))
  labels= c("1st Cluster", 
            "2nd Cluster",
            "3rd Cluster",
            "4th Cluster",
            "5th Cluster",
            "Other Clusters",
            "Uncertain")
  #ajusta cores e labels para incerto
  cores[qtdPicos+2]<-c(alpha("gray50",1))
  labels[qtdPicos+2]<-c("Uncertain")
  
  g2<-ggplot()+
    ggtitle("GMM clusterization")+
    theme_bw()+
    xlab(xtext)+
    ylab("Counts") +
    scale_y_continuous(name = expression("Counts"),
                       limits = c(0, max(valHist$counts)))+
    geom_point(data=valHist, aes(x=x, 
                                   y=counts,
                                   color= cor), pch = 21)+
    scale_color_manual("Order by N",values = cores,
                       labels= labels)
    
  g[[1]]<-g1
  g[[2]]<-g2
  subG<-NA
  return(list(result,g,subG,max(densidade$yDens)))
}



#Print and show the detection chart
printGraf<-function(dirFig,
                    g,
                    subG,
                    gene,
                    minExpression){
  

  library(cowplot)
  gFinal <- plot_grid(g[[1]], g[[2]],g[[3]], ncol=1,  align = "v")

  ggsave(filename = paste0(dirFig,"/",gene,".pdf"),
         plot = gFinal,
         device = "pdf",
         width = 11,height = 8,
         dpi=600)
  
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
  coordPicos<-merge(valPicos, densidade, by.x = "x", by.y = "xNorm")
  coordPicos<-coordPicos[coordPicos$tipo == 1,c(4,5)]
  colnames(coordPicos)<-c("x","y")
  d<-densidade
  densidade$picos <- merge(densidade[c("xNorm")], valPicos, by = 1, all.x = T)[,3]
  
  qtdPicos <-nrow(valPicos[valPicos$tipo==1,])
  
  return(list(densidade,
              qtdPicos,
              maxExpr,
              minExpr,
              medExpe,
              sdExpe,
              coordPicos))
  
}

#prepare and plot graphic 3
graph3<- function(densidade2,
                  qtdPicos,
                  coordPicos,
                  maxExpr,
                  minExpr, 
                  resultado){
  
  # preparação de dados para o gráfico 3
  #cria df com o valor dos histogramas por cluster
  valHist2 <- data.frame(x =      numeric(), 
                         counts = numeric(),
                         cor =    character(),
                         stringsAsFactors = F)
  
  breaks<-(seq(0,maxExpr,(maxExpr-minExpr)/nrow(densidade2)))
  breaks[length(breaks)]<-maxExpr
  #calcula histograma completo
  #totaliza resultados e histogramas
  i=1
  for(i in 1:(qtdPicos)){
    #calcula histograma do cluster
    tmp <- hist(resultado[[i]]$expression, 
                breaks = breaks ,
                plot = F)
    tmp<-data.frame(x=breaks[2:length(breaks)], 
                    counts=tmp$counts,
                    stringsAsFactors = F)
    tmp<-tmp[tmp$counts!=0,]
    if(nrow(tmp) == 0){
      next
    }else{
      tmp$cor<-as.character(i)
      valHist2<-rbind(valHist2,tmp)
    }
  }  
  
  #grafico 3 ----
  xtext<-"Expression Values"
  #define cores e labels
  cores<-c("1"=alpha("red",0.75),
           "2"=alpha("blue",0.75),
           "3"= alpha("green4",0.75),
           "4"= alpha("green3",0.75),
           "5"= alpha("lightblue4",0.75),
           "6"= alpha("lightsalmon2",0.75),
           "7"=alpha("gray50",1))
  labels= c("1st Cluster", 
            "2nd Cluster",
            "3rd Cluster",
            "4th Cluster",
            "5th Cluster",
            "Other Clusters",
            "Uncertain")
  #ajusta cores e labels para incerto
  cores[qtdPicos+2]<-c(alpha("gray50",1))
  labels[qtdPicos+2]<-c("Uncertain")
  
  if(nrow(valHist2)==0){
    text = paste("No graphic to show")
    g3<-ggplot() + 
      annotate("text", x = 4, y = 25, size=8, label = text) + 
      theme_void()
  }else{
    g3<-ggplot()+
      ggtitle("Peaks after filtration")+
      theme_bw()+
      xlab(xtext)+
      ylab("Density") +
      scale_y_continuous(name = expression("Density"), 
                         #limits = c(0, max(densidade$yDens)),
                         sec.axis = sec_axis(~ ./max(densidade2$yDens)* max(valHist2$counts) , 
                                             name = "Counts"))+
      geom_line(data=densidade2, aes(xDens,
                                     yDens,
                                     linetype="New Density",
                                     color="New Density"))+
      geom_point(data=valHist2, aes(x=x, 
                                    y=(counts*max(densidade2$yDens)/ max(counts)),
                                    color= cor,
                                    shape= "Samples count"))
    if(qtdPicos!=0){
      g3<-g3+geom_point(data = coordPicos,aes(x,
                                              y,
                                              shape="Peaks", 
                                              color= "Peaks",
                                              fill = "Peaks"))+
        
        scale_fill_manual("",guide=F,values = c("Peaks"="magenta2"))+
        scale_shape_manual("",values = c("Samples count"=21,
                                         "Peaks"=25))+
        scale_color_manual("",guide=F,values = c(cores,
                                                 "Peaks"="magenta2", 
                                                 "New Density"= "orange"))+
        scale_linetype_manual("",values = c("New Density"=1))+
        guides(shape = guide_legend(override.aes = list(shape = c(25,21), 
                                                        fill=c("magenta2",NA), 
                                                        color= c("magenta2", "black"))),
               linetype = guide_legend(override.aes = list(color = "orange")))
    }
  }
  
}

#save the samples from clusters and peaks
saveSamples<-function(dirSamplesDest,
                      resultado,
                      coordPicos, 
                      gene){
  #correlation between peaks and clusters
  #ensure peak coord order
  coordPicos<- coordPicos[order(coordPicos$x),]
  coordPicos$peak<- seq(1:nrow(coordPicos))
  coordPicos$cluster<- NA
  
  i=4
  for(i in 1:(length(resultado[[1]])-1)){
    #Find correlation peak / cluster
    minVal <- min(resultado[[i]]$expression)
    maxVal <- max(resultado[[i]]$expression)
    j=2
    for(j in 1:nrow(coordPicos)){
      if(minVal <= coordPicos$x[j] & 
         coordPicos$x[j] <= maxVal){
        coordPicos$cluster[j]<-i
      }  
    }
  }
  
  

  if(nrow(coordPicos[is.na(coordPicos$cluster),])>0){
    if(1 %in% coordPicos$cluster){
      coordPicos$cluster[is.na(coordPicos$cluster)]<-2
    }else{
      coordPicos$cluster[is.na(coordPicos$cluster)]<-1
    }
  }

  coordPicos<-coordPicos[!duplicated(coordPicos$cluster),]
  
  #salva identificador das amostras
  i=2
  for(i in 1:length(resultado)){
    if(i %in% coordPicos$cluster){
      peak<-coordPicos$peak[coordPicos$cluster == i]
      idt<- paste0("_peak",peak,"_cluster",i)
    }else{
      idt<- paste0("_cluster",i)
      
    }
    
    write.table(as.character(resultado[[i]]$sample),
                paste0(dirSamplesDest,"/",gene,"_",minExpression,idt,".lst"),
                row.names = F, 
                col.names = F)
    #totaliza as contagens
    
  }
  return(coordPicos) 
}

#write the summary
writeSummary <- function(dirSamplesDest,
                         resultado,
                         coordPicos,
                         amAbaixo, 
                         gene,
                         qtdPicos2){
  soma<-0
  lastCl<- length(resultado)
  nCl<-seq(1,(length(resultado)-1))
  summFile<-paste0(dirSamplesDest,"/",gene,"_",minExpression,".summary")
  
  cat(paste("Below min Expression: ",nrow(amAbaixo),"\n"),file=summFile, append = F)
  soma<-soma + nrow(amAbaixo)
  
  if(qtdPicos2 > nrow(coordPicos)){
    qtdPicos2 <- nrow(coordPicos)
  }
  #write clusters that correspond to peaks
  i=1
  for(i in 1:qtdPicos2){
    clust<-coordPicos$cluster[coordPicos$peak == i]
    cat(paste0("Peak ",i,": ",
               nrow(resultado[[clust]])," (cluster ",clust,")\n"),
        file=summFile, append = T)
    nCl<-nCl[nCl != clust]
    soma<-soma + nrow(resultado[[clust]])
  }
  #other clusters
  for(i in nCl){
    cat(paste0("Cluster ",i,": ",
               nrow(resultado[[i]]),"\n"),
        file=summFile, append = T)
    soma<-soma + nrow(resultado[[i]])
  }
  
  #Uncertain samples
  cat(paste0("Uncertain: ",
             nrow(resultado[[lastCl]]),"\n"),
      file=summFile, append = T)
  soma<-soma + nrow(resultado[[lastCl]])
  
  cat(paste0("\nTotal samples: ",
             soma
             ,"\n"),
      file=summFile, append = T)
}


#Read samples file and start processing
processa = function(dirBase, 
                    dirFig, 
                    atenuacao,
                    minExpression,
                    minSampleSize,
                    minClusterSize,
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
  #genes<-genes[which(genes == "MIS18A"):length(genes)]
  #gene="KIRREL3-AS2"
  
  #loop ----
  #Prossesing each genes in sample
  for(gene in genes){
    #le amostras ----
    amostras <- na.omit(allAmostras[allAmostras$symbol == gene,
                                    c("sample","symbol","FPKM")])
    colnames(amostras)<-c("V1","V2","V3")
    if(nrow(amostras) == 0){
      writeLines(paste("Error processing ",gene,": no expressions values found."))
      cat(paste("Error processing",gene,":  no expressions values found.","\n"),file=arqLog,append = T)
      next()}
    #nome do gene ----
    #genes com alias tem nomes separados por "//" e devem ser eliminados
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
    if(nrow(amostras)<minSampleSize){
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
    coordPicos<-pPicos[[7]]
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
                             minClusterSize = minClusterSize,
                             atenuacao = atenuacao, 
                             tipo = tipo, 
                             qtdPicos = qtdPicos,
                             coordPicos = coordPicos)
      
      if(class(resultado) == "numeric"){
        if(resultado == 0){
          writeLines(paste("Clusters with insufficient samples ",gene," (less then 10% of total)"))
          cat(paste("Clusters with insufficient samples ",gene," (less then 10% of total)","\n"),file=arqLog,append = T)
        }else if (resultado == 1){
          writeLines(paste("GMM fail to clusterize ",gene))
          cat(paste("GMM fail to clusterize ",gene,"\n"),file=arqLog,append = T)
        }
        next()}
      #extrai conteudo do retorno da funcao
      g<-resultado[[2]]
      subG<-resultado[[3]]
      maxYDens<-resultado[[4]]
      resultado<-resultado[[1]]
      
      
      #confere se há mesmo dois picos 
      #cria novo dataset
      amostrasFiltrado<-resultado[[1]][0,]
      for(idx in 1:qtdPicos){
        amostrasFiltrado<-rbind(amostrasFiltrado,resultado[[idx]])
      }
      amostrasTmp<-amostrasFiltrado[,1:2]
      #ajusta df
      amostrasTmp$V2<-gene
      colnames(amostrasTmp)<-c("V1","V3","V2")
      amostrasTmp<-amostrasTmp[,c("V1","V2","V3")]
      
      #Test if have more than 2 peaks after filtering
      #Novo teste de picos 
      pPicos2<-processaPicos(amostrasTmp,arqLog,gene)
      densidade2<-pPicos2[[1]]
      qtdPicos2<-pPicos2[[2]]
      coordPicos<-pPicos2[[7]]
      #graph 3 ----
      g[[3]]<-graph3(densidade2 = densidade2,
                    qtdPicos = qtdPicos2,
                     coordPicos = coordPicos,
                     maxExpr = maxExpr,
                     minExpr = minExpr,
                    resultado = resultado)
      
      if(qtdPicos2>1){
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
      
      #save samples ----
      coordPicos<-saveSamples(dirSamplesDest = dirSamplesDest,
                  resultado = resultado ,
                  coordPicos = coordPicos, 
                  gene = gene)
      
      
      #write summary ----
      writeSummary(dirSamplesDest = dirSamplesDest,
                   resultado = resultado,
                   coordPicos = coordPicos,
                   amAbaixo = amAbaixo,
                   gene = gene,
                   qtdPicos2 = qtdPicos2 )
      

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

#Read samples file and start processing in parallel
processaPar = function(dirBase, 
                       dirFig, 
                       atenuacao,
                       minExpression,
                       minSampleSize,
                       minClusterSize,
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
  #genes<-genes[which(genes == "MIS18A"):length(genes)]
  #genes="BRDT"#
  
  library(doParallel)
  nucleos<-detectCores()-2
  cl <- makeCluster(nucleos)
  registerDoParallel(cl)
  
  cat(paste("Starting the parallel processing for ", tipo,"\n"),file=arqLog,append = T)
  system(paste0('tail -f ',arqLog), wait=F)
  # foreach ----
  foreach(geneIdx = seq(1,length(genes)),
          .combine=rbind,
          .export = c("genes",
                      "arqLog",
                      "allAmostras",
                      "dirBase", 
                      "dirFig", 
                      "atenuacao",
                      "minExpression",
                      "minSampleSize",
                      "minClusterSize",
                      "tipo",
                      "fileName",
                      "dirDados",
                      "limiarDw",
                      "limiarUp",
                      #functions
                      "achaPico",
                      "processaPicos",
                      "aplicaTH",
                      "fderivada",
                      "intervalo",
                      "normaliza",
                      "printGraf",
                      "processaPicos",
                      "graph3",
                      "saveSamples",
                      "writeSummary")) %dopar% {
                        gene<-genes[geneIdx]
                        
                        #le amostras ----
                        amostras <- na.omit(allAmostras[allAmostras$symbol == gene,
                                                        c("sample","symbol","FPKM")])
                        colnames(amostras)<-c("V1","V2","V3")
                        if(nrow(amostras) == 0){
                          writeLines(paste("Error processing ",gene,": no expressions values found."))
                          cat(paste("Error processing",gene,":  no expressions values found.","\n"),file=arqLog,append = T)
                          return()}
                        #nome do gene ----
                        #genes com alias tem nomes separados por "//" e devem ser eliminados
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
                        if(nrow(amostras)<minSampleSize){
                          writeLines(paste("Not enough data - just ", nrow(amostras), " samples for ", gene))
                          cat(paste("Not enough data - just ", nrow(amostras), " samples for ", gene,"\n"),file=arqLog,append = T)
                          return()
                        }
                        pPicos<-processaPicos(amostras,arqLog,gene)
                        densidade<-pPicos[[1]]
                        qtdPicos<-pPicos[[2]]
                        maxExpr<-pPicos[[3]]
                        minExpr<-pPicos[[4]]
                        medExpe<-pPicos[[5]]
                        sdExpe<-pPicos[[6]]
                        coordPicos<-pPicos[[7]]
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
                                                 minClusterSize = minClusterSize,
                                                 atenuacao = atenuacao, 
                                                 tipo = tipo, 
                                                 qtdPicos = qtdPicos,
                                                 coordPicos = coordPicos)
                          
                          if(class(resultado) == "numeric"){
                            if(resultado == 0){
                              writeLines(paste("Clusters with insufficient samples ",gene," (less then 10% of total)"))
                              cat(paste("Clusters with insufficient samples ",gene," (less then 10% of total)","\n"),file=arqLog,append = T)
                            }else if (resultado == 1){
                              writeLines(paste("GMM fail to clusterize ",gene))
                              cat(paste("GMM fail to clusterize ",gene,"\n"),file=arqLog,append = T)
                            }
                            return()}
                          #extrai conteudo do retorno da funcao
                          g<-resultado[[2]]
                          subG<-resultado[[3]]
                          maxYDens<-resultado[[4]]
                          resultado<-resultado[[1]]
                          
                          
                          #confere se há mesmo dois picos 
                          #cria novo dataset
                          amostrasFiltrado<-resultado[[1]][0,]
                          for(idx in 1:qtdPicos){
                            amostrasFiltrado<-rbind(amostrasFiltrado,resultado[[idx]])
                          }
                          amostrasTmp<-amostrasFiltrado[,1:2]
                          #ajusta df
                          amostrasTmp$V2<-gene
                          colnames(amostrasTmp)<-c("V1","V3","V2")
                          amostrasTmp<-amostrasTmp[,c("V1","V2","V3")]
                          
                          #Test if have more than 2 peaks after filtering
                          #Novo teste de picos 
                          pPicos2<-processaPicos(amostrasTmp,arqLog,gene)
                          densidade2<-pPicos2[[1]]
                          qtdPicos2<-pPicos2[[2]]
                          coordPicos<-pPicos2[[7]]
                          #graph 3 ----
                          g[[3]]<-graph3(densidade2 = densidade2,
                                         qtdPicos = qtdPicos2,
                                         coordPicos = coordPicos,
                                         maxExpr = maxExpr,
                                         minExpr = minExpr,
                                         resultado = resultado)
                          
                          if(qtdPicos2>1){
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
                          
                          #save samples ----
                          coordPicos<-saveSamples(dirSamplesDest = dirSamplesDest,
                                                  resultado = resultado ,
                                                  coordPicos = coordPicos, 
                                                  gene = gene)
                          
                          
                          #write summary ----
                          writeSummary(dirSamplesDest = dirSamplesDest,
                                       resultado = resultado,
                                       coordPicos = coordPicos,
                                       amAbaixo = amAbaixo,
                                       gene = gene,
                                       qtdPicos2 = qtdPicos2 )
                          
                          
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
                        
                        
                        #fim foreach ----
                        
                      }
  
  stopCluster(cl)
  #dev.off()
  
  colnames(lista)<-c("Nr","Gene")
  colnames(erros)<-c("Arquivo")
  
  #write.csv(lista,paste0(dirFig,"lista.csv"),row.names = F)
  #write.csv(erros,paste0(dirFig,"erros.csv"),row.names = F)
  #close(arqLog)
  
}
