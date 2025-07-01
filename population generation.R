# Population generation
library(tidyverse)
library(MASS)
library(psych)
library(rlist)

#####################################################################################################################################
# This first section was used to generate population correlation matrices that serve as the "population truth" for each condition.
# Each matrix is a item-level correlation matrix generated based on the conditions specified (item factor loadings, intercorrelations)
# I ran this once, saved all of the matrices into a RDS object (Population Correlation Matrices.rds), and then the actual simulation parts calls upon the respective matrix
# That's why this section is commented out
#####################################################################################################################################
# Population Correlation Matrix based on Factor Matrices and Factor Correlations (Hong, 1999; Gnambs & Staufenbiel, 2016)
# P=L %*% C %*% t(L) + D
# L = k x r factor loading matrix for r factors
# C = r x r correlation matrix for common factors
# D = k x k diagonal matrix of unique variances = I - diag(L %*% C %*% t(L))

# Generate population parameters (correlation matrices)
set.seed(2020)
L <- list(matrix(c(.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,10),
                   rep(0,10),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30),
                 nrow=20,ncol=2), # two factor
          matrix(c(.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,20),
                   rep(0,10),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,10),
                   rep(0,20),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30),
                 nrow=30,ncol=3), # three factor
          matrix(c(.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,30),
                   rep(0,10),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,20),
                   rep(0,20),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,10),
                   rep(0,30),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30),
                 nrow=40,ncol=4), # four factor
          matrix(c(.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,40),
                   rep(0,10),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,30),
                   rep(0,20),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,20),
                   rep(0,30),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30,rep(0,10),
                   rep(0,40),.75,.70,.65,.60,.55,.50,.45,.40,.35,.30),
                 nrow=50,ncol=5)) # five factor

intercor_m <- seq(0,.7,.1)
intercor_sd <- seq(.05,.30,.05)
## Function to generate correlation matrix from a vector
vec2symmat <- function(invec, diag = 1, byrow = TRUE) {
  Nrow <- ceiling(sqrt(2*length(invec)))
  if (!sqrt(length(invec)*2 + Nrow) %% 1 == 0) {
    stop("invec is wrong length to create a square symmetrical matrix")
  }
  mempty <- matrix(0, nrow = Nrow, ncol = Nrow)
  mindex <- matrix(sequence(Nrow^2), nrow = Nrow, ncol = Nrow, byrow = byrow)
  if (isTRUE(byrow)) {
    mempty[mindex[lower.tri(mindex)]] <- invec
    mempty[lower.tri(mempty)] <- t(mempty)[lower.tri(t(mempty))]
  } else {
    mempty[mindex[upper.tri(mindex)]] <- invec
    mempty[lower.tri(mempty)] <- t(mempty)[lower.tri(t(mempty))]
  }
  diag(mempty) <- diag
  mempty
}
##  Function to generate correlation matrices
factor_cormat <- function(m,sd,n){
  corvec <- MASS::mvrnorm(n=n,mu=m,Sigma=sd*sd,empirical=T)
  corvec <- ifelse(corvec>=1,.99,ifelse(corvec<=-1,-.99,corvec))
  vec2symmat(corvec)
}
C <- c(lapply(intercor_m,vec2symmat), # two factor (uniform)
       lapply(intercor_m,function(x) vec2symmat(rep(x,3))), # three factor (uniform)
       list.flatten(lapply(intercor_m,function(m) lapply(intercor_sd,function(sd) factor_cormat(m,sd,3)))), # three factor (varied)
       lapply(intercor_m,function(x) vec2symmat(rep(x,6))), # four factor (uniform)
       list.flatten(lapply(intercor_m,function(m) lapply(intercor_sd,function(sd) factor_cormat(m,sd,6)))), # four factor (varied)
       lapply(intercor_m,function(x) vec2symmat(rep(x,10))), # five factor (uniform)
       list.flatten(lapply(intercor_m,function(m) lapply(intercor_sd,function(sd) factor_cormat(m,sd,10))))) # five factor (varied)
I <- c(rep(list(diag(20)),8), # two factor
       rep(list(diag(30)),56), # three factor
       rep(list(diag(40)),56), # four factor
       rep(list(diag(50)),56)) # five factor
R <- c(lapply(C[1:8],function(x) as.matrix(data.frame(L[1]))%*%as.matrix(x)%*%t(as.matrix(data.frame(L[1])))), # two factor
       lapply(C[9:64],function(x) as.matrix(data.frame(L[2]))%*%as.matrix(x)%*%t(as.matrix(data.frame(L[2])))), # three factor
       lapply(C[65:120],function(x) as.matrix(data.frame(L[3]))%*%as.matrix(x)%*%t(as.matrix(data.frame(L[3])))), # four factor
       lapply(C[121:176],function(x) as.matrix(data.frame(L[4]))%*%as.matrix(x)%*%t(as.matrix(data.frame(L[4]))))) # five factor
P <- vector(mode="list",length=176)
D <- I
for(i in 1:176){
  diag(D[i][[1]]) <- diag(R[i][[1]])
  D[i][[1]] <- I[i][[1]]-D[i][[1]]
  P[i][[1]] <- R[i][[1]]+D[i][[1]]
}
## Manually adjust the machine zero negative eigenvalues to ensure all matrices are positive definite (https://www.mathworks.com/matlabcentral/answers/320134-make-sample-covariance-correlation-matrix-positive-definite)
eigen_counter <- 0 # counter for the number of negative eigenvalues I'm adjusting
matrix_counter <- 0 # counter for the number of matrices I'm adjusting
for(i in 1:176){
  V <- eigen(P[i][[1]])$vectors
  eigen <- eigen(P[i][[1]])$values
  eigen_counter <- length(eigen[eigen<0]) + eigen_counter
  if(length(eigen[eigen<0])>0){
    matrix_counter <- matrix_counter+1
  } else{
    matrix_counter <- matrix_counter
  }
  eigen[eigen<0] <- .0000001
  P[i][[1]] <- V%*%diag(eigen,dim(P[i][[1]])[1],dim(P[i][[1]])[2])%*%t(V)
}
names(P) <- paste0("P",1:176)
saveRDS(P,"Population Correlation Matrices.rds")

