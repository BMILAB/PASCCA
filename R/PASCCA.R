#' @title Data transformation
#'
#' @description  To transform the count data to the log2 scale.
#'
#' @usage PAlog(data)
#'
#' @param data Either a numeric matrix or a data frame.
#'     The expressions of RNA-seq data at poly(A) levels,
#'     exon levels or position levels.
#' @return An object of expression matrix or data frame
#'     the same as the input data.
PAlog <- function(data){
  Te <- data[which(data!=0)]
  Te <- log2(Te)
  data[ which(data==0)] =0
  data[ which(data!=0)] =Te
  return(data)
}


#' @title Significance test
#'
#' @description Shrinkage correlation coefficient analysis
#' of the significance test function.
#'
#' @usage Sccoef_test(r, n, p, q)
Sccoef_test<-function(r, n, p, q)
{
  m <- length(r)
  Q <- rep(0, m)
  P = rep(0, m)
  lambda <- 1
  for (k in m:1) {
    lambda <- lambda * (1 - r[k]^2)
    if(lambda<=0){
      Q[k]<-0
    }else{
      Q[k] <- -log(lambda)
    }
  }
  s <- 0
  for (k in 1:m) {
    Q[k] <- (n - k + 1 - 1/2 * (p + q + 3) + s) * Q[k]
    if((p - k + 1) * (q - k + 1)>0){
      P[k] <- 1 - pchisq(Q[k], (p - k + 1) * (q - k + 1))
    }else{
      P[k] <-0
    }
  }
  Sccoef_test<-t(cbind(r, Q, P))
}

#' @title Compute shrinkage correlation coefficient distance
#'
#' @description To compute shrinkage correlation coefficient distance
#' matrix by parallel.
#'
#' @usage SCCAdist(data)
#'
#' @param data Either a numeric matrix or a data frame.
#'     The expressions of RNA-seq data at poly(A) levels,
#'     exon levels or position levels.
#' @return A shrinkage correlation coefficient distance
#' matrix.
SCCAdist<- function(data){
  chose <<-data
  gene1 <- mldata[[chose[1]]]
  gene2 <- mldata[[chose[2]]]
  g1_g2<-cbind(gene1,gene2)
  g1_g2<-data.frame(g1_g2)
  gene1<-data.frame(gene1)
  Num_gene1<-ncol(gene1)
  Num_gene2<-ncol(gene2)
  Num_gene3<-ncol(g1_g2)
  N=Num_gene3;
  G1=matrix(0,N,J);S1=matrix(0,N,J);D1=matrix(0,N,J)
  S2=matrix(0,N,J);S3=matrix(0,N,J);VARS1=matrix(0,N,J)
  lamda1=matrix(0,N,J);T1=matrix(0,N,J);FRY=matrix(0,N,J)
  G2=matrix(0,N,J);st=matrix(0,J,1);ed=matrix(0,J,1)
  for(i in 1:N){
    for(j in 1:J){
      if(j<=tissue1){
        st[j,]=(j-1)*R[j]+1
        ed[j,]=st[j,]+R[j]-1
      }else{
        st[j,]=ed[j-1,]+1
        ed[j,]=st[j,]+R[j]-1
      }
      G1[i,j]=mean(g1_g2[st[j,]:ed[j,],i])
      S1[i,j]=sum((g1_g2[st[j,]:ed[j,],i]-G1[i,j])^2)/(R[j]-1)
      D1[i,j]=sum(S1[i,j]^2)
      S2[i,j]=(1-((J-2)/D1[i,j]))*S1[i,j]
      S3[i,j]=sum((R[j]-1)/(sum(R)-J)*S1[i,j])
      VARS1[i,j]=sum(((g1_g2[st[j,]:ed[j,],i]-G1[i,j])^2-S1[i,j])^2)*R[j]/((R[j]-1)^3)
      lamda1[i,j]= sum((sum(R)-J-R[j]+1)*VARS1[i,j]/(sum(R)-J))/sum((S1[i,j]-S3[i,j])^2)
      lamda1[which(is.nan(lamda1))]<-0
      if(lamda1[i,j]<1&&lamda1[i,j]>0){
        lamda1[i,j]<-lamda1[i,j]
      }else{
        lamda1[i,j]<-0
      }
      T1[i,j]=(1-lamda1[i,j])*S1[i,j]+lamda1[i,j]*S3[i,j]
      FRY[i,j]=sqrt(T1[i,j]/R[j])
      FRY[FRY==0]<-0.0001
      G2[i,j]=sum(G1[i,j]/FRY[i,j]^2)/sum(1/FRY[i,j]^2)
      G2[which(is.nan(G2))]<-0.0001
      G2[G2==0]<-0.0001
    }
  }
  NUO1=as.matrix(diag(1,N,N))
  for(x in 1:N){
    for(y in 1:N){
      NUO1[x,y]=sum(((G1[x,]- G2[x,])/FRY[x,])*((G1[y,]- G2[y,])/FRY[y,]))/sqrt(sum(((G1[x,]- G2[x,])/FRY[x,])^2)*sum(((G1[y,]- G2[y,])/FRY[y,])^2))
      NUO1[which(is.nan(NUO1))]<-0
    }
  }
  R11=NUO1[1:Num_gene1,1:Num_gene1];
  R22=NUO1[(Num_gene1+1):Num_gene3,(Num_gene1+1):Num_gene3];
  R12=NUO1[1:Num_gene1,(Num_gene1+1):Num_gene3];
  R21=NUO1[(Num_gene1+1):Num_gene3,1:Num_gene1];
  if(qr(R22)$rank>=nrow(R22)&&qr(R11)$rank>=nrow(R11)&&qr(R12)$rank>=nrow(R12)&&qr(R21)$rank>=nrow(R21)){
    A=solve(R11)%*%R12%*%solve(R22)%*%R21
    ev=eigen(A)$values
    cca<-rep(0,length(ev))
    for(s in 1:length(ev) )
    {
      if(Im(ev[s])==0&&Re(ev[s])>0){
        cca[s]<-sqrt(Re(ev[s]))
      }else{
        cca[s]<-0
      }
    }
  }else{
    cca<-0
  }
  Sccoef_test <- Sccoef_test(r=cca,n=nsample,p=ncol(gene1),q=ncol(gene2))
  Pnum <- Sccoef_test[3,]
  Pnum=ifelse(Pnum<=alpha & Pnum!=0  , -(log(Pnum)), 0)
  numerator <- sum(Pnum[])
  if (numerator==0) {
    dmat[chose[1],chose[2]]<-0
    dmat[chose[2],chose[1]]<-0
  } else {
    denomiator <- sum(Sccoef_test[1,]*Pnum)
    dmat[chose[1],chose[2]]<-1-denomiator /numerator
    dmat[chose[2],chose[1]]<-1-denomiator /numerator
  }
}


#' @title Produce poly(A) list
#'
#' @description This function can produce a poly(A) list.
#'
#' @usage PAlist(data)
#'
#' @param data Either a numeric matrix or a data frame.
#'     The expressions of RNA-seq data at poly(A) levels,
#'     exon levels or position levels.
#' @return PAlist returns an object of poly(A) list.
PAlist<-function (datafile)
{
  library("plyr")
  mm <- datafile
  mm$gene<- as.character(mm$gene)
  lldata <- function(data) {
    dl <- (data[,3:ncol(mm)])
    rownames(dl) = data[,1]
    return(dl)
  }
  mldata <<- dlply(.data = mm, .variables = "gene", .fun = lldata)
  mldata <<- sapply(mldata, t)
  return(mldata)
}




#' @title  APA-related gene expression data preprocessing
#'
#' @description  It鈥檚 to do preliminary processing of APA-related gene
#'     expression data.User can choose whether to perform log-transformation
#'     by parameter log (default: TRUE).
#'
#' @details This function provides two steps taken to pre-process data. 1) Data cleaning:
#' To filter out genes with one poly(A) site from gene expression data. 2) Data
#' transformation: To transform the count data to the log2 scale by setting the
#' parameter log=TRUE (default: TRUE).
#'
#' @usage PAprocess(data, log = TRUE)
#'
#' @note If the set of samples have repeated measurements, the order of the
#'     samples in data must be arranged in from big to small according to
#'     the number of replicates.
#'
#' @param data Either a numeric matrix or a data frame.The expressions
#'     of RNA-seq data  at exon levels or position levels. The first column
#'     is poly(A) or exon names, the second column is gene names, and the
#'     remaining column is sample names under different biological conditions
#'     such as different tissues, cell types and developmental stages.
#'     If the set of samples have repeated measurements, the order of the
#'     samples in data must be arranged in from big to small according to
#'     the number of replicates.
#' @param log Choose Whether or not to perform log2 transformation for data.
#' @return An object of expression matrix or data frame the same as the input
#'     data.
#' @author Yuqi Long, Wenbin Ye
#' @references Thomas, P. E., X. Wu, M. Liu, B. Gaffney, G. Ji et al.,
#'     2012 Genome-Wide Control of Polyadenylation Site Choice by CPSF30
#'     in Arabidopsis. Plant Cell 24: 4376-4388.
#' @examples
#'  ##Loading example data
#'  data(polyA_example_data2)
#'  dim(data2)
#'  str(data2)
#'  data2[1:5,1:5]
#'
#'  ##Data preprocessing
#'  pre_data <- PAprocess(data2,log="T")
#'  dim(pre_data)
#'  pre_data[1:5,1:5]
PAprocess <- function(data,log=TRUE){
  pgene=data[,2]
  geneRepeat <- table(pgene)
  geneChoose <- names(geneRepeat[(geneRepeat>1)])
  geneJu   <- is.element(pgene,geneChoose)
  data <- data[which(geneJu==TRUE),]
  if (log){
    log_data=apply( data[,3:ncol(data)],2, PAlog)
    log_data=cbind( data[,1:2],log_data)
  }else{
    log_data= data
  }
  return(log_data)
}


#' @title  Distance matrix computation
#'
#' @description  This function computes and returns the distance matrix by
#' using the PASCCA distance measure to computer the distances
#' between the genes (rows) of a APA-related gene experssion
#' data with replicates.
#'
#' @details The function PASCCA uses the shrinkage canonical correlation
#' analysis to calculate the relationship between genes. The input
#' of this function is a matrix or a data frame, for example,
#' it can be the expressions of RNA-seq data  at exon levels or
#' position levels.The output of this function is a distance matrix.
#'
#' @usage PASCCA(data, alpha = 0.05, repli, tissues, tiss)
#'
#' @param data The APA-related gene expression, the first column
#'    is poly(A) or exon names, the second column is gene names,
#'    and the remaining column is sample names under different
#'    biological conditions such as different tissues, cell types
#'    and developmental stages. If the set of samples have repeated
#'    measurements, the order of the samples in data must be arranged
#'    in from big to small according to the number of replicates.
#' @param alpha The cut-off value of the significance level.
#'    We accept the null hypothesis if the significance level
#'    is above the cut-off value. It means the confidence interval
#'    is 0.95 when the alpha is 0.05. The default value of alpha is 0.05.
#' @param repli The numbers of replicates per biological condition
#'    such as different tissues, cell types and developmental stages.
#'    Note that it needs to be in the same order as the input.
#' @param tissues The total number of biological conditions.
#'    If the input data consists of root with three biological
#'    replicates, seed with three biological replicates and flower
#'    with two biological replicates, the tissues will be three
#'    because there are three conditions (root,seed and flower).
#' @param tiss  The frequency of the first type of repetition.
#'    If the input data consists of root with three biological
#'    replicates, seed with three biological replicates and
#'    flower with two biological replicates, the tiss will be
#'    two since both root and seed have three biological replicates.
#' @return PASCCA returns an object of distance matrix.
#' @author Yuqi Long, Wenbin Ye
#' @references Yao J, Chang C, Salmi M L, et al. Genome-scale cluster
#'    analysis of replicated microarrays using shrinkage correlation
#'    coefficient[J]. BMC bioinformatics, 2008, 9(1): 288.
#'
#'    Hong S, Chen X, Jin L, et al. Canonical correlation analysis
#'    for RNA-seq co-expression networks[J]. Nucleic acids research,
#'    2013, 41(8): e95-e95.
#' @examples
#'   ##example1---------------------------------------------
#'   ##Loading example data
#'   data(polyA_example_data2)
#'   dim(data2)
#'
#'   ##Data preprocessing
#'   pre_data <- PAprocess(data2,log=TRUE)
#'   dim(pre_data)
#'
#'   ##Getting information of the samples
#'   sample_name <- colnames(pre_data)[3:ncol(pre_data)]
#'   sample_name <- strsplit(sample_name,"\\d$")
#'   sample_name <- paste("",lapply(sample_name,"[[",1),sep="");
#'   table(sample_name)
#'
#'   ##Calculationg PASCCA distance matrix
#'   gene_dist <- PASCCA(pre_data, alpha = 0.05,
#'                    repli=as.numeric(table(sample_name)),
#'                    tissues=length(unique(sample_name)),
#'                    tiss=as.numeric(table(table(sample_name))))
#'   str(gene_dist)
#'   gene_dist[1:3,1:3]
#'   #or
#'   gene_dist <- PASCCA(pre_data, alpha = 0.05,repli=c(rep(3,14)), tissues=14, tiss=14)
#'
#'   ##Example2---------------------------------------------
#'   data(polyA_example_data1)
#'   pre_data <- PAprocess(data1,log=TRUE)
#'   gene_dist <- PASCCA(pre_data, alpha = 0.05,repli=c(rep(3,13),rep(2,1)), tissues=14, tiss=13)
#'   ##Example3---------------------------------------------
#'   data(polyA_example_data3)
#'   pre_data <- PAprocess(data3,log=TRUE)
#'   gene_dist <- PASCCA(pre_data, alpha = 0.05,repli=c(rep(4,1),rep(3,8),rep(2,5)), tissues=14, tiss=1)
PASCCA<-function (data,alpha=0.05,repli,tissues,tiss)
{
  alpha <<- alpha
  R<<-repli
  J<<-tissues
  tissue1<<-tiss
  mldata <- PAlist(data)
  Listnumber <- length(mldata)
  N <- choose(Listnumber, 2)
  combnumber = matrix(0, N, 2)
  i = 0
  for (j in 1:(Listnumber-1)) {
    for (k in (j+1):Listnumber) {
      i= i+1
      combnumber[i,1] = j
      combnumber[i,2] = k
    }
  }
  dmat <<- matrix(0, nrow = Listnumber, ncol = Listnumber)
  rownames(dmat) = names(mldata)
  colnames(dmat) = names(mldata)
  nsample <<- nrow(mldata[[1]])
  cl.cores <- detectCores()-2
  cl<- makeCluster(cl.cores)
  clusterExport(cl = cl, envir = .GlobalEnv, varlist = c("Sccoef_test","nsample", "alpha", "SCCAdist", "mldata", "dmat","R","J","tissue1"))
  results = parApply(cl,X=as.array(combnumber), MARGIN = 1,SCCAdist)
  stopCluster(cl)
  i = 0
  for (j in 1:(Listnumber - 1)) {
    for (k in (j + 1):Listnumber) {
      i = i + 1
      dmat[j, k] = dmat[k, j] = results[i]
    }
  }
  return(dmat)
}



#' @title Cluster analysis
#'
#' @description Distances of all gene pairs obtained by the function PASCCA,
#' then the distance matrix is further used for clustering by
#' the function PASCCluster. We adopted the widely-used
#' clustering method, hierarchical clustering, to cluster genes,
#' which was implemented by the R function using hclust default
#' parameters.
#'
#' @details This function PASCCluster has three parameters, dist nc and plot.
#' Based on hierarchical clustering, dist is a distance matrix,
#' nc is the number of clusters and plot is used to plot a cluster
#' denfrogram.
#'
#' @usage PASCCluster(dist,nc,plot)
#'
#' @param dist a dissimilarity matrix as produced by the function PASCCA.
#' @param nc umeric scalar (OR a vector) with the number of clusters
#'     the tree should be cut into.
#' @param plot plot clustering tree of a hierarchical clustering
#'     if the value is TRUE (default: FALSE)
#' @return PASCCluster returns a list, including an object of
#'     class [stats::hclust] which describes the tree produced by the
#'     clustering process and a vector with group memberships
#'     by [stats::cutree]. Besides, when the parameter plot is TRUE,
#'     it will generate the a dendrogram.
#' @author Yuqi Long, Wenbin Ye
#' @examples
#'  ##Loading example data
#'  data(polyA_example_data2)
#'
#'  ##Data preprocessing
#'  pre_data <- PAprocess(data2,log=TRUE)
#'
#' ##Getting information of the samples
#' sample_name <- colnames(pre_data)[3:ncol(pre_data)]
#' sample_name <- strsplit(sample_name,"\\d$")
#' sample_name <- paste("",lapply(sample_name,"[[",1),sep="");
#' table(sample_name)
#' ##Getting the number of repetitions per sample
#' sample_replicates <- as.numeric(table(sample_name))
#' sample_replicates <- sample_replicates[order(sample_replicates,decreasing = TRUE)]
#'
#' ##Calculationg PASCCA distance matrix
#' gene_dist <- PASCCA(pre_data, alpha = 0.05,
#'                    repli=sample_replicates,
#'                    tissues=length(unique(sample_name)),
#'                    tiss=sum(sample_replicates==sample_replicates[1]))
#'
#' ##Hierarchical clustering
#' gene_cluster <- PASCCluster(gene_dist,nc=5,plot = TRUE)

PASCCluster <- function(dist,nc=6,plot=FALSE){
  genedist=as.dist(dist)
  genehc = hclust(genedist, method="complete")
  if(plot){
    plot(genehc,hang = -1, cex = 0.6)
  }
  cluster <- cutree(genehc,nc)
  PAhclu=list(hclust= genehc,cutree = cluster)
  return(PAhclu)
}
