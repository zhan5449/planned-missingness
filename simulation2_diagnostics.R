# SF vs. PM Simulation
#library(tidyverse)#takes too long to download in MSI, calling separate libraries
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(MASS)
library(psych)
library(mice)
library(miceadds)
library(doParallel)
library(parallel)
library(doRNG)
library(rlist)
#SET ss=spreadhseet2[] to sample condition to test
#run inside of SFvsPM2 function to get dat_val_pm
P <- readRDS("Population Correlation Matrices.rds") # read in the population matrices object

# two-factor
# first create table into which results are going to be written
spreadsheet_2 <- data.frame("Pop_Condition"=character(), # identifier for population condition (i.e., P1 has a different true population correlation matrix than P2)
                            "Sample_Condition"=character(), # identifier for sample condition (i.e., S1 has a different NFac-SampleSize-NMissing combination than S2)
                            "NFac"=integer(), # number of factors/constructs
                            "SampleSize"=integer(), # sample size (below I use "validation sample" to differentiate from other draws from the population, but it's really just the sample of the condition)
                            "NMissing"=integer()) # number of items missing per factor (total number of items per factor is 10 so 1 missing = 9-item factor)
spreadsheet_2[1:640,"Pop_Condition"] <- rep(paste0("P",1:8),each=80)
spreadsheet_2[1:640,"Sample_Condition"] <- rep(paste0("S",1:80),8)
spreadsheet_2[1:640,"NFac"] <- rep(2,640) #always 2
spreadsheet_2[1:640,"SampleSize"] <- rep(c(100,200,300,400,500,600,700,800,900,1000),64)
spreadsheet_2[1:640,"NMissing"] <- rep(rep(1:8,each=10),8)
##reponse to reviewers point 2- increasing iterations not adding inferential validity
#showcase for pop condition of true cor .40 (P6), with highest 80% missing, at n=500
#how does increasing iterations affect 1) autocorr, 2) psr, 3) bias, 4) coverage
spreadsheet_2_cond <- rbind(spreadsheet_2[415,],
                            spreadsheet_2[435, ],
                            spreadsheet_2[455,],
                            spreadsheet_2[475, ])
#write_csv(spreadsheet_2_cond,"spreadsheet_2_cond.csv")
# ss = spreadsheet_2[475, ] #p6, n500, 80% missing
# ss= spreadsheet_2[455,] #p6, n500, 60%miss
# ss = spreadsheet_2[435, ] #p6, n500, 40% missing
# ss= spreadsheet_2[415,] #p6, n500, 20%miss

# ss=spreadsheet_2[475,]
#iteration=1

run_diagnostics <- function(ss) {

  #get the data
  pop_param <- P[ss$Pop_Condition][[1]] # get the population correlation matrix for this specific condition
  
  
  # Population truth--source of truth against which Short Form A, Short Form B, Short Form C, and PM are going to be compared against
  cor_truth <- psych::cluster.cor(keys=matrix(c(rep(1,10),rep(0,10),
                                                rep(0,10),rep(1,10)),
                                              nrow=20,ncol=2),r.mat=P[ss$Pop_Condition][[1]])$cor 
  cor_pop <- cor_truth[1,2]
  
  dat_val <- data.frame(MASS::mvrnorm(n=ss$SampleSize,mu=rep(0,20),Sigma=pop_param,empirical=F)) # randomly draw validation sample for analysis based on sample size in each condition
  
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
  
  
true.effect <- cor_pop
m <- sum(is.na(dat_val_pm)) / (length(dat_val_pm)*nrow(dat_val_pm)) 
m <- m*100
if(ss$NMissing >=6) iter <- 250 else iter <- 100

#empty tibble to save diagnostics
diagnostic_df <-
  data.frame(
    pop_condition = rep(ss$Pop_Condition,iter),
    sample_condition = rep(ss$Sample_Condition,iter),
    no_factors = rep(ss$NFac,iter),
    sample_size = rep(ss$SampleSize, iter),
    no_missing =rep(ss$NMissing, iter),
    cor_pop = rep(cor_pop,iter),
   iteration_no = integer(length = iter),
    cor_pooled = integer(length = iter),
     fmi = integer(length = iter),
     coverage = vector(length = iter), #coverage rate
     RB = integer(length = iter), #bias estimates, RB=raw bias, PB=percent bias
     PB = integer(length = iter),
     ac = integer(length = iter), #Autocorrelation
    psrf = integer(length = iter) #potential scale reduction factor
  )

 for (i in 1:iter) { 
  
   if (i == 1) {
    imp <- mice::mice(dat_val_pm,method = "pmm", m = m, maxit = 3, print = FALSE, seed = simulation)
  }   else {
    imp <- mice.mids(imp, maxit = 1, print = FALSE, seed = simulation)
    }
   #add one iteration for each loop **have to start at maxit=3 for convergence() to work

  diagnostic_df$iteration_no[i] <- imp$iteration
 
   #compute intercor
   imp_long <- mice::complete(imp,action="long",include=TRUE) ## convert imp object into long format, includes original data with NA values, .imp=0!
   imp_long$c1 <- rowMeans(imp_long[,paste("X",1:10,sep="")]) # create composite scores, first 100 are NA for og data
   imp_long$c2 <- rowMeans(imp_long[,paste("X",(1:10+10),sep="")])

   
   imp_pm <- mice::as.mids(imp_long) ## convert back into mice object, this loses the iteration info so i cannot plot convergence estimates for the composite scores
   #in imp1_pm, chain means and vars are lost. nmis for c1/c2= 100 which is ALL of the original data since these are means
 
   #get the pooled corr across m (this is what charlene saves)
   cor_mi_full <- miceadds::micombine.cor(mi.res=imp_pm,variables=c("c1","c2"))
   diagnostic_df$cor_pooled[i] <- cor_mi_full$r[1]
  
   diagnostic_df$fmi[i] <- cor_mi_full$fmi[1] #save the fmi for cor(c1,c2)

   

    #coverage rate from Oberman 2021 code
    diagnostic_df$coverage[i] <- cor_mi_full$lower95[1] < true.effect & true.effect < cor_mi_full$upper95[1] #coverage
   #will return logical result

    diagnostic_df$RB[i] <- diagnostic_df$cor_pooled[i] - true.effect
    diagnostic_df$PB[i] <- 100 * abs((diagnostic_df$cor_pooled[i]  - true.effect)/ true.effect) #inf when true is zero
    
    
    #convergence estimates
    conver <- mice::convergence(imp) #data frame
    conver <- dplyr::filter(conver, .it==diagnostic_df$iteration_no[i]) #take the ac and psrf for the current iteration
    
      
    diagnostic_df$ac[i] =  mean(conver$ac, na.rm = T)
 
    diagnostic_df$psrf[i] =  mean(conver$psrf, na.rm = T)


 }

return(diagnostic_df)

}


#run the four ss conditions in larger loop?
local_cluster <- makeCluster(15) #on local, set to 7 (nathan laptop, use 15). 
registerDoParallel(local_cluster)
for(simulation in 1:100){
  diag_imp_2 <- foreach(i=1:4,.export=c("run_diagnostics","spreadsheet_2_cond","mice.mids"),.combine=rbind,.options.RNG=simulation) %dorng%
  run_diagnostics(ss=spreadsheet_2_cond[i,]) #each row is for a condition
  write_csv(diag_imp_2,paste0("two-factor diagnostics/i",simulation,".csv"))
}
stopCluster(local_cluster)







