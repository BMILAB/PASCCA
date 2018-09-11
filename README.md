PASCCA R package
====================

Cluster analysis of replicated alternative polyadenylation data using shrinkage canonical correlation analysis

About
====================
The package PASCCA is an easy-to-use R package for analyses of APA related gene expression, including the characterization of poly(A) sites, quantification of association between genes with/without repeated measurements, clustering of APA-related genes to infer significant APA specific gene modules, and the evaluation of clustering performance with a variety of indexes. By providing a better treatment of the noise inherent in repeated measurements and taking into account multiple layers of poly(A) site data, PASCCA could be a general tool for clustering and analyzing APA-specific gene expression data. For the convenience of users, we provide synthetic [poly(A) site data sets](https://github.com/BMILAB/Synthetic-PolyA-Site-Data). 

Installing PASCCA
=============
Mandatory 
---------

* R (>3.1). [R 3.3.3](https://www.r-project.org/) is recommended.

Required R Packages
---------
* [plyr](https://CRAN.R-project.org/package=plyr), [parallel](https://www.rdocumentation.org/packages/parallel)

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/PASCCA")
library(PASCCA)
```

Using PASCCA
=============
In order to facilitate user understanding, we use the provided example dataset to illustrate the standard analysis work-flow of PASCCA. Please refer to the [User Guide](https://github.com/BMILAB/PASCCA/tree/master/doc) for full details.

Section 1 Data preprocessing
---------
Data preprocessing is an important step in the data mining process. First, we use the function "PAprocess" to do preliminary processing of APA-related gene expression data.
```
##Loading example data
data(polyA_example_data2)
##Data preprocessing
pre_data <- PAprocess(data2,log=TRUE)
```

Section 2 Distance matrix computation
---------
Second, based on processed expression data, to calculate the distance between genes with multiple poly(A) sites by the function "PASCCA" with seven parameters (data, alpha, repli, tissues, tiss).
* data: the APA-related gene expression.
* alpha: the cut-off value of the significance level.
* repli: the numbers of replicates per biological condition such as different tissues, cell types and developmental stages.
* tissues: the total number of biological conditions.
* tiss: the frequency of the first type of repetition.
```
##Getting information of the samples
sample_name <- colnames(pre_data)[3:ncol(pre_data)]
sample_name <- strsplit(sample_name,"\\d$")
sample_name <- paste("",lapply(sample_name,"[[",1),sep="");
table(sample_name)
##Getting the number of repetitions per sample
sample_replicates <- as.numeric(table(sample_name))
sample_replicates <- sample_replicates[order(sample_replicates,decreasing = TRUE)]

##Calculationg PASCCA distance matrix
gene_dist <- PASCCA(pre_data, alpha = 0.05,
                    repli=sample_replicates,
                    tissues=length(unique(sample_name)),
                    tiss=sum(sample_replicates==sample_replicates[1]))
#or
gene_dist <- PASCCA(pre_data, alpha = 0.05,repli=c(rep(3,14)), tissues=14, tiss=14)
```

Section 3 Clustering analysis
---------
Distances of all gene pairs obtained by the function "PASCCA", then the distance matrix is further used for clustering by the function "PASCCluster". We adopted the widely-used clustering method, hierarchical clustering, to cluster genes, which was implemented by the R function using "hclust" default parameters. PASCCluster returns a list, including an object of class hclust which describes the tree produced by the clustering process and a vector with group memberships by "cutree".
```
##Hierarchical clustering
gene_cluster <- PASCCluster(gene_dist,nc=5,plot = TRUE)
```
