rm(list = ls())
setwd("/Users/Zozo/Documents/Project/Network entropy code/mouse")
load(sprintf('CoreMouseSTD.rda'))
mouseSTDCorex[1:2,1:3]
#REARRANGE FOR date
RawDate=colnames(mouseSTDCorex)
Genes=rownames(mouseSTDCorex)
head(RawDate)
Date1=gsub(pattern = "1.*", replacement = "1", x = RawDate)
Date2=gsub(pattern = "2.*", replacement = "2", x = Date1)
Date=gsub(pattern = "3.*", replacement = "3", x = Date2)
head(Date)
ID <-replace(Date,Date== "X1","1")
ID <-replace(ID, ID=="X2", "2")
ID <-replace(ID, ID== "X3", "3")
CellType=as.numeric(unlist(ID))
CellType
#if using a sub matrix
#bloodCorex=submatrix
#CellType=subcell

#correlation matrix for genes (with p value for each pair)
library("Hmisc")
corex=mouseSTDCorex
typeof(corex)
#normalize expression matrix
Ncorex=t(apply(corex, 1, function(x)(x-min(x))/(max(x)-min(x))))
#correlation
res2 <- rcorr(as.matrix(t(Ncorex)))
#res2<-as.matrix(res2)
str(res2)
#res2[1:2]
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P
#res2[is.na(res2)] <- 0

#A simple function to format the correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
dim(res2$r)
sum(is.na(res2$r))
#res2$r[is.na(res2$r)] <- 0
#res2$r=abs((res2$r))
#res2$p[is.na(res2$p)] <- 0
res2$p=abs(res2$p)
min(res2$r)
max(res2$r)
outputCor=flattenCorrMatrix(res2$r, res2$P)
outputCor[1:5,1:4]
sum(is.na(res2$r))
#res2$r[is.na(res2$r)] <- 0
#res2$p[is.na(res2$p)] <- 0
res2$p=abs(res2$p)
min(res2$r)
max(res2$r)
outputCor=flattenCorrMatrix(res2$r, res2$P)
#outputCor[is.na(outputCor)] <- 0
outputCor[1:5,1:4]
head(outputCor[,3])
#find the pencentile of the correlation
#hist(abs(outputCor[,3]))
#indexNA=which(is.na(res2$r))
#res2$r[indexNA[1]]

#check the quantile for correlation
correlationDistribution=abs(outputCor[,3])
quan=quantile(correlationDistribution, c(.9, .95, .999)) 
quan=as.numeric(quan)
max(correlationDistribution)
min(correlationDistribution)
mean(correlationDistribution)
hist(correlationDistribution)
str(correlationDistribution)

#Make the correlation network binary 
#make network binary by quantile
networks=symnum(abs(res2$r), cutpoints = c(0,0.5,1),
                symbols = c("0", "1"),
                abbr.colnames = FALSE)
dim(networks)
#check total number of connection
sum(as.numeric(networks))

#change the network to numeric
networks[1:3,1:3]
typeof(networks)
network=matrix(as.numeric(unlist(networks)),ncol=ncol(networks),byrow=FALSE)
is.numeric(network)
dim(network)
#network=abs(network)
colnames(network) = Genes
rownames(network) = Genes
network[1:3,1:3]

#degree of each node
GeneConnectivity=rowSums(network)-1
hist(GeneConnectivity)
NetworkThreshold=mean(GeneConnectivity)
#find the total connection each node(gene) 
listB=rowSums(network)
names(listB)=rownames(network)
head(listB)
max(listB)
#find the genes with high connectivity
maxA=sum(listB>=NetworkThreshold)
maxA


#find top genes
#get index of the top 2000 genes
df=as.data.frame(listB)
head=df[order(df, decreasing=TRUE)[1:2000],,drop=FALSE]
head(head)
topGenes=rownames(head)
topIndex=matrix(1,ncol=length(topGenes),1)
for (i in 1:length(topGenes)){
  topIndex[i]=which(rownames(network)==topGenes[i])
}
topIndex=as.numeric(unlist(topIndex))
#find the matching gene names
Genes=rownames(corex)
Genes
corex[1:3,1:4]
cells=CellType
topGenes=Genes[topIndex]
head(topGenes)
#write.csv(topGenes,"TopBlood.csv")
#topgene connectivity matrix
str(topIndex)
topConnectivity=network[topIndex,topIndex]
topCorex=corex[topIndex,]
topCorex[1:3,1:4]
dim(topCorex)
#change the expression matrix to binary
#make the expression matrix binary
expression <- function(geneMatrix) {
  #creat a empty matrix,dim(geneMatrix)[2] is #cells
  nCell=dim(geneMatrix)[2]
  nGene=dim(geneMatrix)[1]
  a <-matrix(0, nrow = nGene, ncol = nCell)
  #calculate entropy for each gene in the cell
  for(i in 1:nCell) {
    for(j in 1:nGene) {
      #      p=sum(geneMatrix[,i]) # let mean be the threshold
      #p=mean(geneMatrix[,i]) # let mean be the threshold
      #p=quantile(geneMatrix[,i], 0.75) 
      p=0.1
      if(geneMatrix[j,i] == 0){
        a[j,i] <- 0 
      } 
      if(geneMatrix[j,i] <= p){
        a[j,i] <- 0 
      }
      
      else{
        a[j,i] <- 1
      }
      
      
    }
  }
  data.frame(a)
}
#topCorexB=expression(topCorex)
corexB=expression(Ncorex)
dim(corex)
sum(corexB)

corexB[1:3,1:4]
#save(topCorexB,file="topbinaryExpressionMatrixNorm01.rda")
#save(corexB,file="allbinaryExpressionMatrixNorm01mouseSTD.rda")
load(sprintf('allbinaryExpressionMatrixNorm01mouseSTD.rda'))

#find functional modules from correlation matrix
library("igraph")
#if using the whole 5000 genes
network[1:3,1:4]
topConnectivity=network
dim(topConnectivity)
#if using the top 100 genes
topConnectivity=matrix(unlist(topConnectivity),ncol=ncol(topConnectivity),byrow=FALSE)
rownames(topConnectivity) <- Genes
colnames(topConnectivity) <- Genes
#rownames(topConnectivity) <- Genes[topIndex]
#colnames(topConnectivity) <- Genes[topIndex]
g <- graph_from_adjacency_matrix(topConnectivity, weighted=NULL,mode="undirected")
wtc <- cluster_walktrap(g)
#degree matrix
degree_distribution(g)
#only works for undirected graphs, which this example is fine since symetric
fc <- fastgreedy.community(as.undirected(g))
#csg <- cluster_spinglass(g,vertex=5)
ml <- cluster_louvain(g, weights = NULL)
#save(fc,file="fcM1.14.rda")


comm=membership(ml)
max(comm)
unique(comm)

#make colors for different communities
V(g)$color <- ifelse(membership(fc)==1,"red",
                      ifelse(membership(fc)==2,"orange",
                             ifelse(membership(fc)==3,"yellow","green")))


#plot(g2)
#plot(g2)
topGenes=Genes
index1=which(comm==1)
comm1=topGenes[index1]
#comm1
index2=which(comm==2)
comm2=topGenes[index2]
index3=which(comm==3)
comm3=topGenes[index3]
index4=which(comm==4)
comm4=topGenes[index4]
max(comm)
commcsv=sort(comm)
index5=which(comm==5)
comm5=topGenes[index5]
index6=which(comm==6)
index7=which(comm==7)
index8=which(comm==8)
index9=which(comm==9)
typeof(index9)
commall=c(comm1,comm2,comm3,comm4,comm5)
commidx=c(rep(1,length(comm1)),rep(2,length(comm2)),rep(3,length(comm3)),rep(4,length(comm4)),rep(5,length(comm5)))
commcsv=rbind(commall,commidx)
#write.csv(commcsv,"TopCardioCommConnect.csv")

#calculation connectivity of each module
moduleConnecitivity <- function(comm,topCorex){
  numModule=max(comm)
  nGene=length(comm)
  nCell=dim(topCorex)[2]
  mcMatrix=matrix(0,ncol=nCell,nrow=numModule)
  ifactivate=matrix(0,ncol=nCell,nrow=numModule)
  #ts is the mean expression
  t1=mean(rowSums(mouseSTDCorex))/ncol(mouseSTDCorex)
  for (j in 1:nCell){
    for (i in 1:numModule){
      NumNodes=sum(comm==i)
      N=which(comm==i)
      #t1=as.numeric(quantile(unlist(topCorex[N,]),0.5))
      #t2=as.numeric(quantile(unlist(topCorex[N,]),0.05))
      MC=0 #sum of gene expressed
      for (n in 1:NumNodes){
        iN=N[n]
        MC=MC+topCorex[iN,j]
      }
      #mcMatrix[i,j]=MC/NumNodes # normalized pathway expression 【0-1】
      mcMatrix[i,j]=MC # total pathway expression 【0-1】
      mcvalue=as.numeric(MC)
      if (mcvalue>=t1){
        ifactivate[i,j]=1
      }
      else{
        ifactivate[i,j]=0
      }
      
      #mcMatrix[i,j]=MC
    }
  }
  data.frame(rbind(mcMatrix,ifactivate))
  
}
#calculation connectivity of each module using normalized expression
NmoduleConnecitivity <- function(comm,topCorex){
  numModule=max(comm)
  nGene=length(comm)
  nCell=dim(topCorex)[2]
  mcMatrix=matrix(0,ncol=nCell,nrow=numModule)
  ifactivate=matrix(0,ncol=nCell,nrow=numModule)
  for (j in 1:nCell){
    for (i in 1:numModule){
      NumNodes=sum(comm==i)
      N=which(comm==i)
      t1=0.1
      t2=0.01
      #t1=as.numeric(quantile(unlist(topCorex[N,]),0.5))
      #t2=as.numeric(quantile(unlist(topCorex[N,]),0.05))
      MC=0 #sum of gene expressed
      for (n in 1:NumNodes){
        iN=N[n]
        MC=MC+topCorex[iN,j]
      }
      mcMatrix[i,j]=MC/NumNodes # normalized pathway expression 【0-1】
      #mcMatrix[i,j]=MC # total pathway expression 【0-1】
      mcvalue=as.numeric(MC)
      if (mcvalue>=t1){
        ifactivate[i,j]=1
      }
      else if (mcvalue<=t1 && mcvalue>=t2 && NumNodes>=10){
        ifactivate[i,j]=1
      }
      else{
        ifactivate[i,j]=0
      }
      
      #mcMatrix[i,j]=MC
    }
  }
  data.frame(rbind(mcMatrix,ifactivate))
  
}

#mcMatrix=moduleConnecitivity(comm,Ncorex)
mcMatrix=moduleConnecitivity(comm,corex)
NmcMatrix=moduleConnecitivity(comm,Ncorex)

mcMatrix[3:6,1:2]
#save(mcMatrix,file='M1_7Connectivity.rda')
#save(mcMatrix,file='M1_9Connectivity.rda')
#save(mcMatrix,file='M1_12Connectivity.rda')
save(mcMatrix,file='M2_1mlConnectivitymouseSTD.rda')
#save(mcMatrix,file='M1_15Connectivity.rda')
#load(sprintf('M1_14Connectivity.rda'))
#load(sprintf('M1_14Connectivity.rda'))

mcMatrix[1:4,1:3]
mcMatrix[-1:-2,1:3]

corexB[1:2,1:3]
dim(mcMatrix)
max(comm)
mcM=mcMatrix[1:6,]
ifa=mcMatrix[7:12,]
#mcMatrix=rbind(mcM,ifa)
typeof(mcMatrix)
dim(ifa)
mcM[1:13,1:5]
ifa[1:13,1:5]
#ifa[is.na(ifa)] <- 0
#mcMatrix=mcMatrix[1:10,]
MConnectivity=matrix(unlist(ifa),ncol=ncol(ifa),byrow = FALSE)
dim(MConnectivity)
MC=colSums(MConnectivity)
ct=unlist(CellType)
boxplot(MC ~ ct,  xlab = "Celltype",
        ylab = "ModularConnectivity",main = "mouseSTD",varwidth=TRUE )

McConnectivity=matrix(unlist(mcM),ncol=ncol(mcM),byrow = FALSE)
dim(McConnectivity)
McC=colSums(McConnectivity)
boxplot(McC ~ ct,  xlab = "Celltype",
        ylab = "Connectivity",main = "Cardiomyocyte",varwidth=TRUE )


load(sprintf('cardioCounts.rda'))
CountGene=colSums(bloodCorex)
boxplot(CountGene ~ ct,  xlab = "Celltype",
        ylab = "CellConnectivity",main = "Celltype",varwidth=TRUE )


meanGene <- function(mcMatrix,CellType){
  ihsc=which(CellType=='1')
  impp=which(CellType=='2')
  imep=which(CellType=='3')
  icmp=which(CellType=='4')
  igmp=which(CellType=='5')
  
  D0=mean(mcMatrix[ihsc])
  D2=mean(mcMatrix[impp])
  D5=mean(mcMatrix[imep])
  D15=mean(mcMatrix[icmp])
  D30=mean(mcMatrix[igmp])
  
  data.frame(cbind(D0,D2,D5,D15,D30))
}
meanGene(MC,ct)
