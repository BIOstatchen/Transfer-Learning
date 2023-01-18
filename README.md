# GAMM：trans-ethnic GRS-informed gene-based association mixed model
# Background
Due to ethnic heterogeneity in genetic structure, genetic risk scores (GRS) constructed in European populations typically possess poor portability in underrepresented non-European populations. We aim here to explore the predictive performance of the GRS for traits by exploiting existing knowledge of genetic similarity obtained from Europeans. In this work, we refer to the primary population (e.g. EUR) as the auxiliary population and the underrepresented population (e.g. AFR/EAS) as the target population


**(https://github.com/BIOstatchen/Transfer-Learning)** is implemented in R statistical environment.

# Example
For the parameter estimation in GAMM
```ruby
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("lmm_pxem2_ZPcisSNP_variance_NotationPXEM_EM.cpp")
source("LRTSim_lmm_PXEM_Rcpp.R")
source("estimate_beta.R")
BETA=read.table("BETA.txt",head=T) ## BETA is the effect size of SNPs in base population from summary statistics and SE is their standard error.
SE=read.table("SE.txt",head=T)
R=read.table("R.txt",head=T) ## R is the linkage disequilibrium (LD) matrix which can be calculated with genotypes of population-matched individuals from external reference panels such as the 1000 Genomes Project in the base population. We calculate R in a shrinkage fashion R = 0.95 * as.matrix(cor(G_base_population)) + diag(1-0.95, nrow=nSNPs, ncol=nSNPs)
G=read.table("G.txt",head=T)
y=read.table("y.txt",head=T)
X=read.table("X.txt",head=T)
BETA=as.matrix(BETA)
SE=as.matrix(SE)
R=as.matrix(R)
G=as.matrix(G)
y=as.matrix(y)
X=as.matrix(X)
weight=estimate_beta(BETA,SE,R)
g = as.matrix(apply(G%*%weight,2,scale)[,1])
fit = lmm_pxem2_ZPcisSNP(y, X=cbind(1, X, g), G=G, PXEM=TRUE, maxIter=1000)

v1 = var(g%*%fit$alpha[4,1])
v2 = var(G%*%fit$mub)
v3 = fit$sigma2e
v4 = var(cbind(1,X)%*%fit$alpha[1:3,1])
pve = (v1+v2)/(v1+v2+v3+v4)
pge = v1/(v1+v2)

//' @param y  response variable for GWAS data
//' @param X  covariates for GWAS data
//' @param g  GRS = G*weight is the trans-enthic GRS information
//' @param G  genotype (cis-SNPs) matrix for GWAS
//' @param maxIter  maximum iteration (default is 1000)

```
For joint effect test using LRT in GAMM
```ruby
library(glmnet)
library(Rcpp)
library(RcppArmadillo)
source("LRTSim_lmm_PXEM_Rcpp.R")
sourceCpp("lmm_pxem2_ZPcisSNP_variance_NotationPXEM_EM.cpp")
sourceCpp("LRTSim_lmm_PXEM_Rcpp.cpp")
source("estimate_beta.R")
BETA=read.table("BETA.txt",head=T) ## BETA is the effect size of SNPs in base population from summary statistics and SE is their standard error.
SE=read.table("SE.txt",head=T)
R=read.table("R.txt",head=T) ## R is the linkage disequilibrium (LD) matrix which can be calculated with genotypes of population-matched individuals from external reference panels such as the 1000 Genomes Project in the base population. We calculate R in a shrinkage fashion R = 0.95 * as.matrix(cor(G_base_population)) + diag(1-0.95, nrow=nSNPs, ncol=nSNPs)
G=read.table("G.txt",head=T)
y=read.table("y.txt",head=T)
X=read.table("X.txt",head=T)
BETA=as.matrix(BETA)
SE=as.matrix(SE)
R=as.matrix(R)
G=as.matrix(G)
y=as.matrix(y)
X=as.matrix(X)
weight=estimate_beta(BETA,SE,R)
g = as.matrix(apply(G%*%weight,2,scale)[,1])
fit = lmm_pxem2_ZPcisSNP(y, X=cbind(1, X, g), G=G, PXEM=TRUE, maxIter=1000)
simLike <- LRTSimZP(Z = X, E = g, G = G, nsim=1e5, parallel=c("multicore"), ncpus = 4L) ## exact LRT

fitx=lm(y~X)
obsLike = c((fit$loglik - logLik(fitx))*2)
plrt = mean(simLike >= obsLike)
fitmix = pmixLRT(simLike, c(1e3,1e4,1e5)) ## approximate LRT

//' @param weight  true effects for SNPs available for the same trait in another base population
//' @param y  response variable for GWAS data
//' @param X  covariates for GWAS data, not include the vector of ones
//' @param g  GRS = G*weight is the trans-enthic GRS information
//' @param G  genotype (cis-SNPs) matrix for GWAS

```

# Cite
Haojie Lu<sup>$</sup>, Shuo Zhang<sup>$</sup>, Zhou Jiang<sup>$</sup> and Ping Zeng<sup>#</sup> (2022). Leveraging trans-ethnic genetic risk scores to improve association power for complex traits in underrepresented populations.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.

# Update
2022-07-22 GAMM version 1.0
