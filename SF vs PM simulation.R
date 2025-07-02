#SF vs PM simulation
library(tidyverse)
library(MASS)
library(psych)
library(mice)
library(miceadds)
library(doParallel)
library(parallel)
library(doRNG)
library(rlist)
library(here)

#####################################################################################################################################
# Actual simulation begins here
# It is broken into four different sections: two-, three-, four-, and five-factors (factor = construct). 
# E.g., two-factor data means it has two constructs, with each construct being measured by 10 items, so it's a 20-item dataset; three-factor = 30-item dataset, etc.
# The base structure of the code is very similar for each of the four sections, each subsequent one is slightly adapted to accommodate the number of constructs/items
# I focused on annotating the two-factor section in a little more detail, and the same comments would apply to the other sections
# This is also how I broke up the jobs for MSI--I had separate scripts for two-, three-, four-, five-construct runs (because the scripts were slightly different and also as kind of a manual parallelization)
#####################################################################################################################################

P <- readRDS("Population Correlation Matrices.rds") # read in the population matrices object
set.seed(2020)

# two-factor
# first create table into which results are going to be written
spreadsheet_2 <- data.frame("Pop_Condition"=character(), # identifier for population condition (i.e., P1 has a different true population correlation matrix than P2)
                            "Sample_Condition"=character(), # identifier for sample condition (i.e., S1 has a different NFac-SampleSize-NMissing combination than S2)
                            "NFac"=integer(), # number of factors/constructs
                            "SampleSize"=integer(), # sample size (below I use "validation sample" to differentiate from other draws from the population, but it's really just the sample of the condition)
                            "NMissing"=integer()) # number of items missing per factor (total number of items per factor is 10 so 1 missing = 9-item factor)
spreadsheet_2[1:640,"Pop_Condition"] <- rep(paste0("P",1:8),each=80)
spreadsheet_2[1:640,"Sample_Condition"] <- rep(paste0("S",1:80),8)
spreadsheet_2[1:640,"NFac"] <- rep(2,640)
spreadsheet_2[1:640,"SampleSize"] <- rep(c(100,200,300,400,500,600,700,800,900,1000),64)
spreadsheet_2[1:640,"NMissing"] <- rep(rep(1:8,each=10),8)

# this is the function that actually performs the simulation
SFvsPM_2 <- function(ss){
  pop_param <- P[ss$Pop_Condition][[1]] # get the population correlation matrix for this specific condition
  
  # Population truth--source of truth against which Short Form A, Short Form B, Short Form C, and PM are going to be compared against
  cor_truth <- psych::cluster.cor(keys=matrix(c(rep(1,10),rep(0,10),
                                         rep(0,10),rep(1,10)),
                                       nrow=20,ncol=2),r.mat=P[ss$Pop_Condition][[1]])$cor 
  cor_pop <- cor_truth[1,2] # based on the item-level correlation matrix, get the construct-level correlation

  # Short Form A (a empirically developed Short Form already exists so the entirety of the sample size in each condition can be used to compute results)
  dat_dev_A <- data.frame(MASS::mvrnorm(n=500,mu=rep(0,20),Sigma=pop_param,empirical=F)) # randomly draw a sufficiently large developmental sample to choose SF items based on EFA loadings
  dat_val <- data.frame(MASS::mvrnorm(n=ss$SampleSize,mu=rep(0,20),Sigma=pop_param,empirical=F)) # randomly draw validation sample for analysis based on sample size in each condition
  efa1_A <- psych::fa(r=dat_dev_A[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_A <- rownames(efa1_A)[order(efa1_A,decreasing=T)][1:(10-ss$NMissing)] # short form items for factor 1
  efa2_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10))],nFactors=1,rotate="oblimin")$loadings
  var_sf2_A <- rownames(efa2_A)[order(efa2_A,decreasing=T)][1:(10-ss$NMissing)] # short form items for factor 2
  
  dat_val_sf_A <- dat_val # the entire val dataset
  dat_val_sf_A$c1 <- rowMeans(dat_val_sf_A[,var_sf1_A])
  dat_val_sf_A$c2 <- rowMeans(dat_val_sf_A[,var_sf2_A])
  cor_sf_A <- cor(dat_val_sf_A[,c("c1","c2")])[1,2] # correlation based on short form A
  
  # Short Form B
  dat_dev_B <- dat_val[sample(c(1:nrow(dat_val)),size=ss$SampleSize/2,replace=F),] # half of the validation sample needs to be used to develop SF
  # I'm essentially extending the same time constraint implied by NMissing for the validation sample to the development sample.
  # e.g., assume NMissing is 5. That means for each of the two factors, participants only have time to complete 5 out of 10 items--in total only have time for 10 items. 
  # this constraint should also apply to the short form development sample such that each participant can only complete 10 items, which means only half of the dev sample can be used for each factor
  # for other values of NMissing, it's not as clean, but the two lines below calculates what that sample size would be 
  dat_dev1_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),] 
  dat_dev2_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  efa1_B <- psych::fa(r=dat_dev1_B[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_B <- rownames(efa1_B)[order(efa1_B,decreasing=T)][1:(10-ss$NMissing)] # short form items for factor 1
  efa2_B <- psych::fa(r=dat_dev2_B[,paste0("X",1:10+10)],nFactors=1,rotate="oblimin")$loadings
  var_sf2_B <- rownames(efa2_B)[order(efa2_B,decreasing=T)][1:(10-ss$NMissing)] # short form items for factor 2
  
  dat_val_sf_B <- dat_val[!row.names(dat_val) %in% row.names(dat_dev_B),] # the half of val dataset that wasn't used for developing SF
  dat_val_sf_B$c1 <- rowMeans(dat_val_sf_B[,var_sf1_B])
  dat_val_sf_B$c2 <- rowMeans(dat_val_sf_B[,var_sf2_B])
  cor_sf_B <- cor(dat_val_sf_B[,c("c1","c2")])[1,2] # correlation based on short form B
  
  # Short Form C
  var_sf1_C <- paste("X",(1:10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)] ## randomly choose indicators for short form
  var_sf2_C <- paste("X",(1:10+10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  
  dat_val_sf_C <- dat_val # the entire val dataset
  dat_val_sf_C$c1 <- rowMeans(dat_val_sf_C[,var_sf1_C])
  dat_val_sf_C$c2 <- rowMeans(dat_val_sf_C[,var_sf2_C])
  cor_sf_C <- cor(dat_val_sf_C[,c("c1","c2")])[1,2] # correlation based on short form C
  
  ## Full
  dat_val_full <- dat_val
  dat_val_full$c1 <- rowMeans(dat_val_full[,paste("X",1:10,sep="")])
  dat_val_full$c2 <- rowMeans(dat_val_full[,paste("X",(1:10)+10,sep="")])
  cor_full <- cor(dat_val_full[,c("c1","c2")])[1,2] # correlation based on the entire validation sample with no missing data (represents sampling error only)
  
  ## PM
  dat_val_pm <- dat_val
  dat_val_pm <- apply(dat_val_pm,1,function(x){
    index <- c(sample(c(1:10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10),size=ss$NMissing,replace=F)) # Identify the missing variables for each row
    x[as.numeric(index)] <- NA
    data.frame(t(x))
  }
  )
  dat_val_pm <- do.call(rbind,dat_val_pm)
  imp <- mice::mice(dat_val_pm,method="pmm",m=5,maxit=5,seed=iteration) 
  imp1 <- imp
  imp1_long <- mice::complete(imp1,action="long",include=TRUE) ## convert imp object into long format
  imp1_long$c1 <- rowMeans(imp1_long[,paste("X",1:10,sep="")]) # create composite scores
  imp1_long$c2 <- rowMeans(imp1_long[,paste("X",(1:10+10),sep="")])
  imp1_pm <- mice::as.mids(imp1_long) ## convert back into mice object
  cor_mi <- miceadds::micombine.cor(mi.res=imp1_pm,variables=c("c1","c2"))[1,"r"] # correlation based on PM
  
  unlist(c(ss,cor_pop,cor_full,cor_sf_A,cor_sf_B,cor_sf_C,cor_mi))
}

# Run SFvsPM_2 Function 
# this is where the function created above is actually run
local_cluster <- makeCluster(10)
registerDoParallel(local_cluster)
for(iteration in 1:100){
  res_Pop_2 <- foreach(i=1:640,.export=c("SFvsPM_2","spreadsheet_2"),.combine=rbind,.options.RNG=iteration) %dorng%
    SFvsPM_2(ss=spreadsheet_2[i,])
  colnames(res_Pop_2)[6:11] <- c("Pop12","Full12","SFA12","SFB12","SFC12","MI12")
  write_csv(data.frame(res_Pop_2),paste0(here("data","Two-Factor Raw Results","i"),iteration,".csv"))
}
stopCluster(local_cluster)

#####################################################################################################################################
# three-factor
spreadsheet_3 <- data.frame("Pop_Condition"=character(),
                            "Sample_Condition"=character(),
                            "NFac"=integer(),
                            "SampleSize"=integer(),
                            "NMissing"=integer())
spreadsheet_3[1:4480,"Pop_Condition"] <- rep(paste0("P",9:64),each=80)
spreadsheet_3[1:4480,"Sample_Condition"] <- rep(paste0("S",1:80),56)
spreadsheet_3[1:4480,"NFac"] <- rep(3,640)
spreadsheet_3[1:4480,"SampleSize"] <- rep(c(100,200,300,400,500,600,700,800,900,1000),448)
spreadsheet_3[1:4480,"NMissing"] <- rep(rep(1:8,each=10),56)

SFvsPM_3 <- function(ss){
  pop_param <- P[ss$Pop_Condition][[1]] # population correlation matrix

  # Population truth
  cor_truth <- psych::cluster.cor(keys=matrix(c(rep(1,10),rep(0,20),
                                                rep(0,10),rep(1,10),rep(0,10),
                                                rep(0,20),rep(1,10)),
                                              nrow=30,ncol=3),r.mat=P[ss$Pop_Condition][[1]])$cor 
  cor_pop <- c(cor_truth[1,2],
               cor_truth[1,3],
               cor_truth[2,3]) # true population intercorrelations
  
  # Short Form A
  dat_dev_A <- data.frame(MASS::mvrnorm(n=500,mu=rep(0,30),Sigma=pop_param,empirical=F)) # randomly draw a sufficiently large developmental sample to choose SF items based on EFA loadings
  dat_val <- data.frame(MASS::mvrnorm(n=ss$SampleSize,mu=rep(0,30),Sigma=pop_param,empirical=F)) # randomly draw validation sample for analysis
  efa1_A <- psych::fa(r=dat_dev_A[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_A <- rownames(efa1_A)[order(efa1_A,decreasing=T)][1:(10-ss$NMissing)]
  efa2_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10))],nFactors=1,rotate="oblimin")$loadings
  var_sf2_A <- rownames(efa2_A)[order(efa2_A,decreasing=T)][1:(10-ss$NMissing)]
  efa3_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10*2))],nFactors=1,rotate="oblimin")$loadings
  var_sf3_A <- rownames(efa3_A)[order(efa3_A,decreasing=T)][1:(10-ss$NMissing)]
  
  dat_val_sf_A <- dat_val # the entire val dataset
  dat_val_sf_A$c1 <- rowMeans(dat_val_sf_A[,var_sf1_A])
  dat_val_sf_A$c2 <- rowMeans(dat_val_sf_A[,var_sf2_A])
  dat_val_sf_A$c3 <- rowMeans(dat_val_sf_A[,var_sf3_A])
  cor_sf_A <- cor(dat_val_sf_A[,c("c1","c2","c3")])
  cor_sf_A <- c(cor_sf_A[1,2],
                cor_sf_A[1,3],
                cor_sf_A[2,3]) # SF A intercorrelations
  
  # Short Form B
  dat_dev_B <- dat_val[sample(c(1:nrow(dat_val)),size=ss$SampleSize/2,replace=F),] # half of the validation sample to develop SF
  dat_dev1_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),] # a proportion (determined by PMissing) of half the val data for SF
  dat_dev2_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  dat_dev3_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  efa1_B <- psych::fa(r=dat_dev1_B[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_B <- rownames(efa1_B)[order(efa1_B,decreasing=T)][1:(10-ss$NMissing)]
  efa2_B <- psych::fa(r=dat_dev2_B[,paste0("X",1:10+10)],nFactors=1,rotate="oblimin")$loadings
  var_sf2_B <- rownames(efa2_B)[order(efa2_B,decreasing=T)][1:(10-ss$NMissing)] 
  efa3_B <- psych::fa(r=dat_dev3_B[,paste0("X",1:10+10*2)],nFactors=1,rotate="oblimin")$loadings
  var_sf3_B <- rownames(efa3_B)[order(efa3_B,decreasing=T)][1:(10-ss$NMissing)] 
  
  dat_val_sf_B <- dat_val[!row.names(dat_val) %in% row.names(dat_dev_B),] # the half of val dataset that wasn't used for developing SF
  dat_val_sf_B$c1 <- rowMeans(dat_val_sf_B[,var_sf1_B])
  dat_val_sf_B$c2 <- rowMeans(dat_val_sf_B[,var_sf2_B])
  dat_val_sf_B$c3 <- rowMeans(dat_val_sf_B[,var_sf3_B])
  cor_sf_B <- cor(dat_val_sf_B[,c("c1","c2","c3")])
  cor_sf_B <- c(cor_sf_B[1,2],
                cor_sf_B[1,3],
                cor_sf_B[2,3])
  
  # Short Form C
  var_sf1_C <- paste("X",(1:10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)] ## randomly choose indicators for short form
  var_sf2_C <- paste("X",(1:10+10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  var_sf3_C <- paste("X",(1:10+10*2),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  
  dat_val_sf_C <- dat_val # the entire val dataset
  dat_val_sf_C$c1 <- rowMeans(dat_val_sf_C[,var_sf1_C])
  dat_val_sf_C$c2 <- rowMeans(dat_val_sf_C[,var_sf2_C])
  dat_val_sf_C$c3 <- rowMeans(dat_val_sf_C[,var_sf3_C])
  cor_sf_C <- cor(dat_val_sf_C[,c("c1","c2","c3")])
  cor_sf_C <- c(cor_sf_C[1,2],
                cor_sf_C[1,3],
                cor_sf_C[2,3])
  
  ## Full
  dat_val_full <- dat_val
  dat_val_full$c1 <- rowMeans(dat_val_full[,paste("X",1:10,sep="")])
  dat_val_full$c2 <- rowMeans(dat_val_full[,paste("X",(1:10)+10,sep="")])
  dat_val_full$c3 <- rowMeans(dat_val_full[,paste("X",(1:10)+10*2,sep="")])
  cor_full <- cor(dat_val_full[,c("c1","c2","c3")])
  cor_full <- c(cor_full[1,2],
                cor_full[1,3],
                cor_full[2,3])
  
  ## PM
  dat_val_pm <- dat_val
  dat_val_pm <- apply(dat_val_pm,1,function(x){
    index <- c(sample(c(1:10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10*2),size=ss$NMissing,replace=F)) ## Identify the missing variables for each row
    x[as.numeric(index)] <- NA
    data.frame(t(x))
  }
  )
  dat_val_pm <- do.call(rbind,dat_val_pm)
  imp <- mice::mice(dat_val_pm,method="pmm",m=40,maxit=5,seed=iteration)
  imp1 <- imp
  imp1_long <- mice::complete(imp1,action="long",include=TRUE) ## convert imp object into long format
  imp1_long$c1 <- rowMeans(imp1_long[,paste("X",1:10,sep="")]) # create composite scores
  imp1_long$c2 <- rowMeans(imp1_long[,paste("X",(1:10+10),sep="")])
  imp1_long$c3 <- rowMeans(imp1_long[,paste("X",(1:10+10*2),sep="")])
  imp1_pm <- mice::as.mids(imp1_long) ## convert back into mice object
  cor_mi <- miceadds::micombine.cor(mi.res=imp1_pm,variables=c("c1","c2","c3"))[1:3,"r"]
  
  unlist(c(ss,cor_pop,cor_full,cor_sf_A,cor_sf_B,cor_sf_C,cor_mi))
}

## Run SFvsPM_3 Function 
local_cluster <- makeCluster(120)
registerDoParallel(local_cluster)
for(iteration in 1:100){
  res_Pop_3 <- foreach(i=1:4480,.export=c("SFvsPM_3","spreadsheet_3"),.combine=rbind,.options.RNG=iteration) %dorng%
    SFvsPM_3(ss=spreadsheet_3[i,])
  colnames(res_Pop_3)[6:23] <- c("Pop12","Pop13","Pop23",
                                 "Full12","Full13","Full23",
                                 "SFA12","SFA13","SFA23",
                                 "SFB12","SFB13","SFB23",
                                 "SFC12","SFC13","SFC23",
                                 "MI12","MI13","MI23")
  write_csv(data.frame(res_Pop_3),paste0(here("data","Three-Factor Raw Results","i"),iteration,".csv"))
}
stopCluster(local_cluster)

#####################################################################################################################################
# four-factor
spreadsheet_4 <- data.frame("Pop_Condition"=character(),
                            "Sample_Condition"=character(),
                            "NFac"=integer(),
                            "SampleSize"=integer(),
                            "NMissing"=integer())
spreadsheet_4[1:4480,"Pop_Condition"] <- rep(paste0("P",65:120),each=80)
spreadsheet_4[1:4480,"Sample_Condition"] <- rep(paste0("S",1:80),56)
spreadsheet_4[1:4480,"NFac"] <- rep(4,640)
spreadsheet_4[1:4480,"SampleSize"] <- rep(c(100,200,300,400,500,600,700,800,900,1000),448)
spreadsheet_4[1:4480,"NMissing"] <- rep(rep(1:8,each=10),56)

SFvsPM_4 <- function(ss){
  pop_param <- P[ss$Pop_Condition][[1]] # population correlation matrix
  
  # Population truth
  cor_truth <- psych::cluster.cor(keys=matrix(c(rep(1,10),rep(0,30),
                                                rep(0,10),rep(1,10),rep(0,20),
                                                rep(0,20),rep(1,10),rep(0,10),
                                                rep(0,30),rep(1,10)),
                                              nrow=40,ncol=4),r.mat=P[ss$Pop_Condition][[1]])$cor 
  cor_pop <- c(cor_truth[1,2],
               cor_truth[1,3],
               cor_truth[1,4],
               cor_truth[2,3],
               cor_truth[2,4],
               cor_truth[3,4]) # true population intercorrelations
  
  # Short Form A
  dat_dev_A <- data.frame(MASS::mvrnorm(n=500,mu=rep(0,40),Sigma=pop_param,empirical=F)) # randomly draw a sufficiently large developmental sample to choose SF items based on EFA loadings
  dat_val <- data.frame(MASS::mvrnorm(n=ss$SampleSize,mu=rep(0,40),Sigma=pop_param,empirical=F)) # randomly draw validation sample for analysis
  efa1_A <- psych::fa(r=dat_dev_A[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_A <- rownames(efa1_A)[order(efa1_A,decreasing=T)][1:(10-ss$NMissing)]
  efa2_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10))],nFactors=1,rotate="oblimin")$loadings
  var_sf2_A <- rownames(efa2_A)[order(efa2_A,decreasing=T)][1:(10-ss$NMissing)]
  efa3_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10*2))],nFactors=1,rotate="oblimin")$loadings
  var_sf3_A <- rownames(efa3_A)[order(efa3_A,decreasing=T)][1:(10-ss$NMissing)]
  efa4_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10*3))],nFactors=1,rotate="oblimin")$loadings
  var_sf4_A <- rownames(efa4_A)[order(efa4_A,decreasing=T)][1:(10-ss$NMissing)]
  
  dat_val_sf_A <- dat_val # the entire val dataset
  dat_val_sf_A$c1 <- rowMeans(dat_val_sf_A[,var_sf1_A])
  dat_val_sf_A$c2 <- rowMeans(dat_val_sf_A[,var_sf2_A])
  dat_val_sf_A$c3 <- rowMeans(dat_val_sf_A[,var_sf3_A])
  dat_val_sf_A$c4 <- rowMeans(dat_val_sf_A[,var_sf4_A])
  cor_sf_A <- cor(dat_val_sf_A[,c("c1","c2","c3","c4")])
  cor_sf_A <- c(cor_sf_A[1,2],
                cor_sf_A[1,3],
                cor_sf_A[1,4],
                cor_sf_A[2,3],
                cor_sf_A[2,4],
                cor_sf_A[3,4]) # SF A intercorrelations
  
  # Short Form B
  dat_dev_B <- dat_val[sample(c(1:nrow(dat_val)),size=ss$SampleSize/2,replace=F),] # half of the validation sample to develop SF
  dat_dev1_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),] # a proportion (determined by PMissing) of half the val data for SF
  dat_dev2_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  dat_dev3_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  dat_dev4_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  efa1_B <- psych::fa(r=dat_dev1_B[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_B <- rownames(efa1_B)[order(efa1_B,decreasing=T)][1:(10-ss$NMissing)]
  efa2_B <- psych::fa(r=dat_dev2_B[,paste0("X",1:10+10)],nFactors=1,rotate="oblimin")$loadings
  var_sf2_B <- rownames(efa2_B)[order(efa2_B,decreasing=T)][1:(10-ss$NMissing)] 
  efa3_B <- psych::fa(r=dat_dev3_B[,paste0("X",1:10+10*2)],nFactors=1,rotate="oblimin")$loadings
  var_sf3_B <- rownames(efa3_B)[order(efa3_B,decreasing=T)][1:(10-ss$NMissing)] 
  efa4_B <- psych::fa(r=dat_dev4_B[,paste0("X",1:10+10*3)],nFactors=1,rotate="oblimin")$loadings
  var_sf4_B <- rownames(efa4_B)[order(efa4_B,decreasing=T)][1:(10-ss$NMissing)] 
  
  dat_val_sf_B <- dat_val[!row.names(dat_val) %in% row.names(dat_dev_B),] # the half of val dataset that wasn't used for developing SF
  dat_val_sf_B$c1 <- rowMeans(dat_val_sf_B[,var_sf1_B])
  dat_val_sf_B$c2 <- rowMeans(dat_val_sf_B[,var_sf2_B])
  dat_val_sf_B$c3 <- rowMeans(dat_val_sf_B[,var_sf3_B])
  dat_val_sf_B$c4 <- rowMeans(dat_val_sf_B[,var_sf4_B])
  cor_sf_B <- cor(dat_val_sf_B[,c("c1","c2","c3","c4")])
  cor_sf_B <- c(cor_sf_B[1,2],
                cor_sf_B[1,3],
                cor_sf_B[1,4],
                cor_sf_B[2,3],
                cor_sf_B[2,4],
                cor_sf_B[3,4])
  
  # Short Form C
  var_sf1_C <- paste("X",(1:10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)] ## randomly choose indicators for short form
  var_sf2_C <- paste("X",(1:10+10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  var_sf3_C <- paste("X",(1:10+10*2),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  var_sf4_C <- paste("X",(1:10+10*3),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  
  dat_val_sf_C <- dat_val # the entire val dataset
  dat_val_sf_C$c1 <- rowMeans(dat_val_sf_C[,var_sf1_C])
  dat_val_sf_C$c2 <- rowMeans(dat_val_sf_C[,var_sf2_C])
  dat_val_sf_C$c3 <- rowMeans(dat_val_sf_C[,var_sf3_C])
  dat_val_sf_C$c4 <- rowMeans(dat_val_sf_C[,var_sf4_C])
  cor_sf_C <- cor(dat_val_sf_C[,c("c1","c2","c3","c4")])
  cor_sf_C <- c(cor_sf_C[1,2],
                cor_sf_C[1,3],
                cor_sf_C[1,4],
                cor_sf_C[2,3],
                cor_sf_C[2,4],
                cor_sf_C[3,4])
  
  ## Full
  dat_val_full <- dat_val
  dat_val_full$c1 <- rowMeans(dat_val_full[,paste("X",1:10,sep="")])
  dat_val_full$c2 <- rowMeans(dat_val_full[,paste("X",(1:10)+10,sep="")])
  dat_val_full$c3 <- rowMeans(dat_val_full[,paste("X",(1:10)+10*2,sep="")])
  dat_val_full$c4 <- rowMeans(dat_val_full[,paste("X",(1:10)+10*3,sep="")])
  cor_full <- cor(dat_val_full[,c("c1","c2","c3","c4")])
  cor_full <- c(cor_full[1,2],
                cor_full[1,3],
                cor_full[1,4],
                cor_full[2,3],
                cor_full[2,4],
                cor_full[3,4])
  
  ## PM
  dat_val_pm <- dat_val
  dat_val_pm <- apply(dat_val_pm,1,function(x){
    index <- c(sample(c(1:10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10*2),size=ss$NMissing,replace=F),
               sample(c((1:10)+10*3),size=ss$NMissing,replace=F)) ## Identify the missing variables for each row
    x[as.numeric(index)] <- NA
    data.frame(t(x))
  }
  )
  dat_val_pm <- do.call(rbind,dat_val_pm)
  imp <- mice::mice(dat_val_pm,method="pmm",m=40,maxit=5,seed=iteration)
  imp1 <- imp
  imp1_long <- mice::complete(imp1,action="long",include=TRUE) ## convert imp object into long format
  imp1_long$c1 <- rowMeans(imp1_long[,paste("X",1:10,sep="")]) # create composite scores
  imp1_long$c2 <- rowMeans(imp1_long[,paste("X",(1:10+10),sep="")])
  imp1_long$c3 <- rowMeans(imp1_long[,paste("X",(1:10+10*2),sep="")])
  imp1_long$c4 <- rowMeans(imp1_long[,paste("X",(1:10+10*3),sep="")])
  imp1_pm <- mice::as.mids(imp1_long) ## convert back into mice object
  cor_mi <- miceadds::micombine.cor(mi.res=imp1_pm,variables=c("c1","c2","c3","c4"))[1:6,"r"]
  
  unlist(c(ss,cor_pop,cor_full,cor_sf_A,cor_sf_B,cor_sf_C,cor_mi))
}

## Run SFvsPM_4 Function 
local_cluster <- makeCluster(120)
registerDoParallel(local_cluster)
for(iteration in 1:100){
  res_Pop_4 <- foreach(i=1:4480,.export=c("SFvsPM_4","spreadsheet_4"),.combine=rbind,.options.RNG=iteration) %dorng%
    SFvsPM_4(ss=spreadsheet_4[i,])
  colnames(res_Pop_4)[6:41] <- c("Pop12","Pop13","Pop14","Pop23","Pop24","Pop34",
                                 "Full12","Full13","Full14","Full23","Full24","Full34",
                                 "SFA12","SFA13","SFA14","SFA23","SFA24","SFA34",
                                 "SFB12","SFB13","SFB14","SFB23","SFB24","SFB34",
                                 "SFC12","SFC13","SFC14","SFC23","SFC24","SFC34",
                                 "MI12","MI13","MI14","MI23","MI24","MI34")
  write_csv(data.frame(res_Pop_4),paste0(here("data","Four-Factor Raw Results","i"),iteration,".csv"))
}
stopCluster(local_cluster)

#####################################################################################################################################
# five-factor
spreadsheet_5 <- data.frame("Pop_Condition"=character(),
                            "Sample_Condition"=character(),
                            "NFac"=integer(),
                            "SampleSize"=integer(),
                            "NMissing"=integer())
spreadsheet_5[1:4480,"Pop_Condition"] <- rep(paste0("P",121:176),each=80)
spreadsheet_5[1:4480,"Sample_Condition"] <- rep(paste0("S",1:80),56)
spreadsheet_5[1:4480,"NFac"] <- rep(5,640)
spreadsheet_5[1:4480,"SampleSize"] <- rep(c(100,200,300,400,500,600,700,800,900,1000),448)
spreadsheet_5[1:4480,"NMissing"] <- rep(rep(1:8,each=10),56)

SFvsPM_5 <- function(ss){
  pop_param <- P[ss$Pop_Condition][[1]] # population correlation matrix
  
  # Population truth
  cor_truth <- psych::cluster.cor(keys=matrix(c(rep(1,10),rep(0,40),
                                                rep(0,10),rep(1,10),rep(0,30),
                                                rep(0,20),rep(1,10),rep(0,20),
                                                rep(0,30),rep(1,10),rep(0,10),
                                                rep(0,40),rep(1,10)),
                                              nrow=50,ncol=5),r.mat=P[ss$Pop_Condition][[1]])$cor 
  cor_pop <- c(cor_truth[1,2],
               cor_truth[1,3],
               cor_truth[1,4],
               cor_truth[1,5],
               cor_truth[2,3],
               cor_truth[2,4],
               cor_truth[2,5],
               cor_truth[3,4],
               cor_truth[3,5],
               cor_truth[4,5]) # true population intercorrelations
  
  # Short Form A
  dat_dev_A <- data.frame(MASS::mvrnorm(n=500,mu=rep(0,50),Sigma=pop_param,empirical=F)) # randomly draw a sufficiently large developmental sample to choose SF items based on EFA loadings
  dat_val <- data.frame(MASS::mvrnorm(n=ss$SampleSize,mu=rep(0,50),Sigma=pop_param,empirical=F)) # randomly draw validation sample for analysis
  efa1_A <- psych::fa(r=dat_dev_A[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_A <- rownames(efa1_A)[order(efa1_A,decreasing=T)][1:(10-ss$NMissing)]
  efa2_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10))],nFactors=1,rotate="oblimin")$loadings
  var_sf2_A <- rownames(efa2_A)[order(efa2_A,decreasing=T)][1:(10-ss$NMissing)]
  efa3_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10*2))],nFactors=1,rotate="oblimin")$loadings
  var_sf3_A <- rownames(efa3_A)[order(efa3_A,decreasing=T)][1:(10-ss$NMissing)]
  efa4_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10*3))],nFactors=1,rotate="oblimin")$loadings
  var_sf4_A <- rownames(efa4_A)[order(efa4_A,decreasing=T)][1:(10-ss$NMissing)]
  efa5_A <- psych::fa(r=dat_dev_A[,paste0("X",(1:10+10*4))],nFactors=1,rotate="oblimin")$loadings
  var_sf5_A <- rownames(efa5_A)[order(efa5_A,decreasing=T)][1:(10-ss$NMissing)]
  
  dat_val_sf_A <- dat_val # the entire val dataset
  dat_val_sf_A$c1 <- rowMeans(dat_val_sf_A[,var_sf1_A])
  dat_val_sf_A$c2 <- rowMeans(dat_val_sf_A[,var_sf2_A])
  dat_val_sf_A$c3 <- rowMeans(dat_val_sf_A[,var_sf3_A])
  dat_val_sf_A$c4 <- rowMeans(dat_val_sf_A[,var_sf4_A])
  dat_val_sf_A$c5 <- rowMeans(dat_val_sf_A[,var_sf5_A])
  cor_sf_A <- cor(dat_val_sf_A[,c("c1","c2","c3","c4","c5")])
  cor_sf_A <- c(cor_sf_A[1,2],
                cor_sf_A[1,3],
                cor_sf_A[1,4],
                cor_sf_A[1,5],
                cor_sf_A[2,3],
                cor_sf_A[2,4],
                cor_sf_A[2,5],
                cor_sf_A[3,4],
                cor_sf_A[3,5],
                cor_sf_A[4,5]) # SF A intercorrelations
  
  # Short Form B
  dat_dev_B <- dat_val[sample(c(1:nrow(dat_val)),size=ss$SampleSize/2,replace=F),] # half of the validation sample to develop SF
  dat_dev1_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),] # a proportion (determined by PMissing) of half the val data for SF
  dat_dev2_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  dat_dev3_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  dat_dev4_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  dat_dev5_B <- dat_dev_B[sample(c(1:nrow(dat_dev_B)),size=ss$SampleSize/2*(1-ss$NMissing/10),replace=F),]
  efa1_B <- psych::fa(r=dat_dev1_B[,paste0("X",1:10)],nFactors=1,rotate="oblimin")$loadings
  var_sf1_B <- rownames(efa1_B)[order(efa1_B,decreasing=T)][1:(10-ss$NMissing)]
  efa2_B <- psych::fa(r=dat_dev2_B[,paste0("X",1:10+10)],nFactors=1,rotate="oblimin")$loadings
  var_sf2_B <- rownames(efa2_B)[order(efa2_B,decreasing=T)][1:(10-ss$NMissing)] 
  efa3_B <- psych::fa(r=dat_dev3_B[,paste0("X",1:10+10*2)],nFactors=1,rotate="oblimin")$loadings
  var_sf3_B <- rownames(efa3_B)[order(efa3_B,decreasing=T)][1:(10-ss$NMissing)] 
  efa4_B <- psych::fa(r=dat_dev4_B[,paste0("X",1:10+10*3)],nFactors=1,rotate="oblimin")$loadings
  var_sf4_B <- rownames(efa4_B)[order(efa4_B,decreasing=T)][1:(10-ss$NMissing)] 
  efa5_B <- psych::fa(r=dat_dev5_B[,paste0("X",1:10+10*4)],nFactors=1,rotate="oblimin")$loadings
  var_sf5_B <- rownames(efa5_B)[order(efa5_B,decreasing=T)][1:(10-ss$NMissing)] 
  
  dat_val_sf_B <- dat_val[!row.names(dat_val) %in% row.names(dat_dev_B),] # the half of val dataset that wasn't used for developing SF
  dat_val_sf_B$c1 <- rowMeans(dat_val_sf_B[,var_sf1_B])
  dat_val_sf_B$c2 <- rowMeans(dat_val_sf_B[,var_sf2_B])
  dat_val_sf_B$c3 <- rowMeans(dat_val_sf_B[,var_sf3_B])
  dat_val_sf_B$c4 <- rowMeans(dat_val_sf_B[,var_sf4_B])
  dat_val_sf_B$c5 <- rowMeans(dat_val_sf_B[,var_sf5_B])
  cor_sf_B <- cor(dat_val_sf_B[,c("c1","c2","c3","c4","c5")])
  cor_sf_B <- c(cor_sf_B[1,2],
                cor_sf_B[1,3],
                cor_sf_B[1,4],
                cor_sf_B[1,5],
                cor_sf_B[2,3],
                cor_sf_B[2,4],
                cor_sf_B[2,5],
                cor_sf_B[3,4],
                cor_sf_B[3,5],
                cor_sf_B[4,5])
  
  # Short Form C
  var_sf1_C <- paste("X",(1:10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)] ## randomly choose indicators for short form
  var_sf2_C <- paste("X",(1:10+10),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  var_sf3_C <- paste("X",(1:10+10*2),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  var_sf4_C <- paste("X",(1:10+10*3),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  var_sf5_C <- paste("X",(1:10+10*4),sep="")[sample(1:10,size=(10-ss$NMissing),replace=F)]
  
  dat_val_sf_C <- dat_val # the entire val dataset
  dat_val_sf_C$c1 <- rowMeans(dat_val_sf_C[,var_sf1_C])
  dat_val_sf_C$c2 <- rowMeans(dat_val_sf_C[,var_sf2_C])
  dat_val_sf_C$c3 <- rowMeans(dat_val_sf_C[,var_sf3_C])
  dat_val_sf_C$c4 <- rowMeans(dat_val_sf_C[,var_sf4_C])
  dat_val_sf_C$c5 <- rowMeans(dat_val_sf_C[,var_sf5_C])
  cor_sf_C <- cor(dat_val_sf_C[,c("c1","c2","c3","c4","c5")])
  cor_sf_C <- c(cor_sf_C[1,2],
                cor_sf_C[1,3],
                cor_sf_C[1,4],
                cor_sf_C[1,5],
                cor_sf_C[2,3],
                cor_sf_C[2,4],
                cor_sf_C[2,5],
                cor_sf_C[3,4],
                cor_sf_C[3,5],
                cor_sf_C[4,5])
  
  ## Full
  dat_val_full <- dat_val
  dat_val_full$c1 <- rowMeans(dat_val_full[,paste("X",1:10,sep="")])
  dat_val_full$c2 <- rowMeans(dat_val_full[,paste("X",(1:10)+10,sep="")])
  dat_val_full$c3 <- rowMeans(dat_val_full[,paste("X",(1:10)+10*2,sep="")])
  dat_val_full$c4 <- rowMeans(dat_val_full[,paste("X",(1:10)+10*3,sep="")])
  dat_val_full$c5 <- rowMeans(dat_val_full[,paste("X",(1:10)+10*4,sep="")])
  cor_full <- cor(dat_val_full[,c("c1","c2","c3","c4","c5")])
  cor_full <- c(cor_full[1,2],
                cor_full[1,3],
                cor_full[1,4],
                cor_full[1,5],
                cor_full[2,3],
                cor_full[2,4],
                cor_full[2,5],
                cor_full[3,4],
                cor_full[3,5],
                cor_full[4,5])
  
  ## PM
  dat_val_pm <- dat_val
  dat_val_pm <- apply(dat_val_pm,1,function(x){
    index <- c(sample(c(1:10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10),size=ss$NMissing,replace=F),
               sample(c((1:10)+10*2),size=ss$NMissing,replace=F),
               sample(c((1:10)+10*3),size=ss$NMissing,replace=F),
               sample(c((1:10)+10*4),size=ss$NMissing,replace=F)) ## Identify the missing variables for each row
    x[as.numeric(index)] <- NA
    data.frame(t(x))
  }
  )
  dat_val_pm <- do.call(rbind,dat_val_pm)
  imp <- mice::mice(dat_val_pm,method="pmm",m=40,maxit=5,seed=iteration)
  imp1 <- imp
  imp1_long <- mice::complete(imp1,action="long",include=TRUE) ## convert imp object into long format
  imp1_long$c1 <- rowMeans(imp1_long[,paste("X",1:10,sep="")]) # create composite scores
  imp1_long$c2 <- rowMeans(imp1_long[,paste("X",(1:10+10),sep="")])
  imp1_long$c3 <- rowMeans(imp1_long[,paste("X",(1:10+10*2),sep="")])
  imp1_long$c4 <- rowMeans(imp1_long[,paste("X",(1:10+10*3),sep="")])
  imp1_long$c5 <- rowMeans(imp1_long[,paste("X",(1:10+10*4),sep="")])
  imp1_pm <- mice::as.mids(imp1_long) ## convert back into mice object
  cor_mi <- miceadds::micombine.cor(mi.res=imp1_pm,variables=c("c1","c2","c3","c4","c5"))[1:10,"r"]
  
  unlist(c(ss,cor_pop,cor_full,cor_sf_A,cor_sf_B,cor_sf_C,cor_mi))
}

## Run SFvsPM_5 Function 
local_cluster <- makeCluster(10)
registerDoParallel(local_cluster)
for(iteration in 1:100){
  res_Pop_5 <- foreach(i=1:4480,.export=c("SFvsPM_5","spreadsheet_5"),.combine=rbind,.options.RNG=iteration) %dorng%
    SFvsPM_5(ss=spreadsheet_5[i,])
  colnames(res_Pop_5)[6:65] <- c("Pop12","Pop13","Pop14","Pop15","Pop23","Pop24","Pop25","Pop34","Pop35","Pop45",
                                 "Full12","Full13","Full14","Full15","Full23","Full24","Full25","Full34","Full35","Full45",
                                 "SFA12","SFA13","SFA14","SFA15","SFA23","SFA24","SFA25","SFA34","SFA35","SFA45",
                                 "SFB12","SFB13","SFB14","SFB15","SFB23","SFB24","SFB25","SFB34","SFB35","SFB45",
                                 "SFC12","SFC13","SFC14","SFC15","SFC23","SFC24","SFC25","SFC34","SFC35","SFC45",
                                 "MI12","MI13","MI14","MI15","MI23","MI24","MI25","MI34","MI35","MI45")
  write_csv(data.frame(res_Pop_5),paste0(here("data","Five-Factor Raw Results","i"),iteration,".csv"))
}
stopCluster(local_cluster)