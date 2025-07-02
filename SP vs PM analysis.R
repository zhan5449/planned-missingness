# SF vs. PM Simulation Analysis
options(scipen=999)
library(openxlsx)
library(matrixStats)
library(sjPlot)
library(wesanderson)
library(ggplot2)
library(broom)
library(ggfortify)
library(plyr)
library(tidyverse)
library(DescTools)
library(here)

# Read and Clean Data
files2F <- ldply(list.files(path=here::here("data","Two-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T),read_csv) %>%
  mutate(iteration=rep(gsub("data/Two-Factor Raw Results/i|.csv","",list.files(path=here::here("data","Two-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T)),each=640))
# write_csv(files2F,here("data","Two-Factor i1-i100.csv")) #already in the data folder

files3F <- ldply(list.files(path=here::here("data","Three-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T),read_csv) %>%
  mutate(iteration=rep(gsub("data/Three-Factor Raw Results/i|.csv","",list.files(path=here::here("data","Three-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T)),each=4480))
# write_csv(files3F,here("data","Three-Factor i1-i100.csv"))

files4F <- ldply(list.files(path=here::here("data","Four-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T),read_csv) %>%
  mutate(iteration=rep(gsub("data/Four-Factor Raw Results/i|.csv","",list.files(path=here::here("data","Four-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T)),each=4480))
# write_csv(files4F,here("data","Four-Factor i1-i100.csv"))

files5F <- ldply(list.files(path=here::here("data","Five-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T),read_csv) %>%
  mutate(iteration=rep(gsub("data/Five-Factor Raw Results/i|.csv","",list.files(path=here::here("data","Five-Factor Raw Results"),pattern="i\\d+.*.csv",full.names=T)),each=4480))
# write_csv(files5F,here("data","Five-Factor i1-i100.csv"))

# files2F <- read_csv(here("data","Two-Factor i1-i100.csv"))
# files3F <- read_csv(here("data","Three-Factor i1-i100.csv"))
# files4F <- read_csv(here("data","Four-Factor i1-i100.csv"))
# files5F <- read_csv(here("data","Five-Factor i1-i100.csv"))

# Merge Data and Create Variables for Analysis
RowSD <- function(x) {
  sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
}
## 2-Factor
cor_res_2 <- files2F %>%
  mutate(PMissing=NMissing/10,
         PopR_mean=Pop12,
         PopR_sd=0) %>%
  # mutate(across(Pop12:MI12,FisherZ)) %>% # uncomment if trying fisher's z transformation on correlations
  mutate(FULL_TRUE=abs(Pop12-Full12),
         SFA_TRUE=abs(Pop12-SFA12),
         SFB_TRUE=abs(Pop12-SFB12),
         SFC_TRUE=abs(Pop12-SFC12),
         MI_TRUE=abs(Pop12-MI12),
         FULL_TRUE_r=Pop12-Full12,
         SFA_TRUE_r=Pop12-SFA12,
         SFB_TRUE_r=Pop12-SFB12,
         SFC_TRUE_r=Pop12-SFC12,
         MI_TRUE_r=Pop12-MI12) 
## 3-Factor
cor_res_3 <- files3F %>%
  mutate(MI12=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI23),NA,MI12),
         MI13=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI23),NA,MI13),
         MI23=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI23),NA,MI23)) %>%
  mutate(PMissing=NMissing/10,
         PopR_mean=rowMeans(cbind(abs(Pop12),abs(Pop13),abs(Pop23))),
         PopR_sd=RowSD(cbind(Pop12,Pop13,Pop23))) %>%
  # mutate(across(Pop12:MI23,FisherZ)) %>% # uncomment if trying fisher's z transformation on correlations
  mutate(FULL_TRUE_12=abs(Pop12-Full12),
         FULL_TRUE_13=abs(Pop13-Full13),
         FULL_TRUE_23=abs(Pop23-Full23),
         SFA_TRUE_12=abs(Pop12-SFA12),
         SFA_TRUE_13=abs(Pop13-SFA13),
         SFA_TRUE_23=abs(Pop23-SFA23),
         SFB_TRUE_12=abs(Pop12-SFB12),
         SFB_TRUE_13=abs(Pop13-SFB13),
         SFB_TRUE_23=abs(Pop23-SFB23),
         SFC_TRUE_12=abs(Pop12-SFC12),
         SFC_TRUE_13=abs(Pop13-SFC13),
         SFC_TRUE_23=abs(Pop23-SFC23),
         MI_TRUE_12=abs(Pop12-MI12),
         MI_TRUE_13=abs(Pop13-MI13),
         MI_TRUE_23=abs(Pop23-MI23),
         FULL_TRUE_r_12=Pop12-Full12,
         FULL_TRUE_r_13=Pop13-Full13,
         FULL_TRUE_r_23=Pop23-Full23,
         SFA_TRUE_r_12=Pop12-SFA12,
         SFA_TRUE_r_13=Pop13-SFA13,
         SFA_TRUE_r_23=Pop23-SFA23,
         SFB_TRUE_r_12=Pop12-SFB12,
         SFB_TRUE_r_13=Pop13-SFB13,
         SFB_TRUE_r_23=Pop23-SFB23,
         SFC_TRUE_r_12=Pop12-SFC12,
         SFC_TRUE_r_13=Pop13-SFC13,
         SFC_TRUE_r_23=Pop23-SFC23,
         MI_TRUE_r_12=Pop12-MI12,
         MI_TRUE_r_13=Pop13-MI13,
         MI_TRUE_r_23=Pop23-MI23) %>%
  mutate(FULL_TRUE=rowMeans(select(.,FULL_TRUE_12,FULL_TRUE_13,FULL_TRUE_23),na.rm=T),
         SFA_TRUE=rowMeans(select(.,SFA_TRUE_12,SFA_TRUE_13,SFA_TRUE_23),na.rm=T),
         SFB_TRUE=rowMeans(select(.,SFB_TRUE_12,SFB_TRUE_13,SFB_TRUE_23),na.rm=T),
         SFC_TRUE=rowMeans(select(.,SFC_TRUE_12,SFC_TRUE_13,SFC_TRUE_23),na.rm=T),
         MI_TRUE=rowMeans(select(.,MI_TRUE_12,MI_TRUE_13,MI_TRUE_23),na.rm=T), # no direction
         FULL_TRUE_r=rowMeans(select(.,FULL_TRUE_r_12,FULL_TRUE_r_13,FULL_TRUE_r_23),na.rm=T),
         SFA_TRUE_r=rowMeans(select(.,SFA_TRUE_r_12,SFA_TRUE_r_13,SFA_TRUE_r_23),na.rm=T),
         SFB_TRUE_r=rowMeans(select(.,SFB_TRUE_r_12,SFB_TRUE_r_13,SFB_TRUE_r_23),na.rm=T),
         SFC_TRUE_r=rowMeans(select(.,SFC_TRUE_r_12,SFC_TRUE_r_13,SFC_TRUE_r_23),na.rm=T),
         MI_TRUE_r=rowMeans(select(.,MI_TRUE_r_12,MI_TRUE_r_13,MI_TRUE_r_23),na.rm=T))

## 4-Factor
cor_res_4 <- files4F %>%
  mutate(MI12=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI23)|is.na(MI24)|is.na(MI34),NA,MI12),
         MI13=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI23)|is.na(MI24)|is.na(MI34),NA,MI13),
         MI14=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI23)|is.na(MI24)|is.na(MI34),NA,MI14),
         MI23=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI23)|is.na(MI24)|is.na(MI34),NA,MI23),
         MI24=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI23)|is.na(MI24)|is.na(MI34),NA,MI24),
         MI34=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI23)|is.na(MI24)|is.na(MI34),NA,MI34)) %>%
  mutate(PMissing=NMissing/10,
         PopR_mean=rowMeans(cbind(abs(Pop12),abs(Pop13),abs(Pop14),abs(Pop23),abs(Pop24),abs(Pop34))),
         PopR_sd=RowSD(cbind(Pop12,Pop13,Pop14,Pop23,Pop24,Pop34))) %>%
  # mutate(across(Pop12:MI34,FisherZ)) %>% # uncomment if trying fisher's z transformation on correlations
  mutate(FULL_TRUE_12=abs(Pop12-Full12),
         FULL_TRUE_13=abs(Pop13-Full13),
         FULL_TRUE_14=abs(Pop14-Full14),
         FULL_TRUE_23=abs(Pop23-Full23),
         FULL_TRUE_24=abs(Pop24-Full24),
         FULL_TRUE_34=abs(Pop34-Full34),
         SFA_TRUE_12=abs(Pop12-SFA12),
         SFA_TRUE_13=abs(Pop13-SFA13),
         SFA_TRUE_14=abs(Pop14-SFA14),
         SFA_TRUE_23=abs(Pop23-SFA23),
         SFA_TRUE_24=abs(Pop24-SFA24),
         SFA_TRUE_34=abs(Pop34-SFA34),
         SFB_TRUE_12=abs(Pop12-SFB12),
         SFB_TRUE_13=abs(Pop13-SFB13),
         SFB_TRUE_14=abs(Pop14-SFB14),
         SFB_TRUE_23=abs(Pop23-SFB23),
         SFB_TRUE_24=abs(Pop24-SFB24),
         SFB_TRUE_34=abs(Pop34-SFB34),
         SFC_TRUE_12=abs(Pop12-SFC12),
         SFC_TRUE_13=abs(Pop13-SFC13),
         SFC_TRUE_14=abs(Pop14-SFC14),
         SFC_TRUE_23=abs(Pop23-SFC23),
         SFC_TRUE_24=abs(Pop24-SFC24),
         SFC_TRUE_34=abs(Pop34-SFC34),
         MI_TRUE_12=abs(Pop12-MI12),
         MI_TRUE_13=abs(Pop13-MI13),
         MI_TRUE_14=abs(Pop14-MI14),
         MI_TRUE_23=abs(Pop23-MI23),
         MI_TRUE_24=abs(Pop24-MI24),
         MI_TRUE_34=abs(Pop34-MI34),
         FULL_TRUE_r_12=Pop12-Full12,
         FULL_TRUE_r_13=Pop13-Full13,
         FULL_TRUE_r_14=Pop14-Full14,
         FULL_TRUE_r_23=Pop23-Full23,
         FULL_TRUE_r_24=Pop24-Full24,
         FULL_TRUE_r_34=Pop34-Full34,
         SFA_TRUE_r_12=Pop12-SFA12,
         SFA_TRUE_r_13=Pop13-SFA13,
         SFA_TRUE_r_14=Pop14-SFA14,
         SFA_TRUE_r_23=Pop23-SFA23,
         SFA_TRUE_r_24=Pop24-SFA24,
         SFA_TRUE_r_34=Pop34-SFA34,
         SFB_TRUE_r_12=Pop12-SFB12,
         SFB_TRUE_r_13=Pop13-SFB13,
         SFB_TRUE_r_14=Pop14-SFB14,
         SFB_TRUE_r_23=Pop23-SFB23,
         SFB_TRUE_r_24=Pop24-SFB24,
         SFB_TRUE_r_34=Pop34-SFB34,
         SFC_TRUE_r_12=Pop12-SFC12,
         SFC_TRUE_r_13=Pop13-SFC13,
         SFC_TRUE_r_14=Pop14-SFC14,
         SFC_TRUE_r_23=Pop23-SFC23,
         SFC_TRUE_r_24=Pop24-SFC24,
         SFC_TRUE_r_34=Pop34-SFC34,
         MI_TRUE_r_12=Pop12-MI12,
         MI_TRUE_r_13=Pop13-MI13,
         MI_TRUE_r_14=Pop14-MI14,
         MI_TRUE_r_23=Pop23-MI23,
         MI_TRUE_r_24=Pop24-MI24,
         MI_TRUE_r_34=Pop34-MI34) %>%
  mutate(FULL_TRUE=rowMeans(select(.,FULL_TRUE_12:FULL_TRUE_34),na.rm=T),
         SFA_TRUE=rowMeans(select(.,SFA_TRUE_12:SFA_TRUE_34),na.rm=T),
         SFB_TRUE=rowMeans(select(.,SFB_TRUE_12:SFB_TRUE_34),na.rm=T),
         SFC_TRUE=rowMeans(select(.,SFC_TRUE_12:SFC_TRUE_34),na.rm=T),
         MI_TRUE=rowMeans(select(.,MI_TRUE_12:MI_TRUE_34),na.rm=T), # no direction
         FULL_TRUE_r=rowMeans(select(.,FULL_TRUE_r_12:FULL_TRUE_r_34),na.rm=T),
         SFA_TRUE_r=rowMeans(select(.,SFA_TRUE_r_12:SFA_TRUE_r_34),na.rm=T),
         SFB_TRUE_r=rowMeans(select(.,SFB_TRUE_r_12:SFB_TRUE_r_34),na.rm=T),
         SFC_TRUE_r=rowMeans(select(.,SFC_TRUE_r_12:SFC_TRUE_r_34),na.rm=T),
         MI_TRUE_r=rowMeans(select(.,MI_TRUE_r_12:MI_TRUE_r_34),na.rm=T))

## 5-Factor
cor_res_5 <- files5F %>%
  mutate(MI12=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI12),
         MI13=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI13),
         MI14=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI14),
         MI15=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI15),
         MI23=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI23),
         MI24=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI24),
         MI25=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI25),
         MI34=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI34),
         MI35=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI35),
         MI45=ifelse(is.na(MI12)|is.na(MI13)|is.na(MI14)|is.na(MI15)|is.na(MI23)|is.na(MI24)|is.na(MI25)|is.na(MI34)|is.na(MI35)|is.na(MI45),NA,MI45)) %>%
  mutate(PMissing=NMissing/10,
         PopR_mean=rowMeans(cbind(abs(Pop12),abs(Pop13),abs(Pop14),abs(Pop15),abs(Pop23),abs(Pop24),abs(Pop25),abs(Pop34),abs(Pop35),abs(Pop45))),
         PopR_sd=RowSD(cbind(Pop12,Pop13,Pop14,Pop15,Pop23,Pop24,Pop25,Pop34,Pop35,Pop45))) %>%
  # mutate(across(Pop12:MI45,FisherZ)) %>% # uncomment if trying fisher's z transformation on correlations
  mutate(FULL_TRUE_12=abs(Pop12-Full12),
         FULL_TRUE_13=abs(Pop13-Full13),
         FULL_TRUE_14=abs(Pop14-Full14),
         FULL_TRUE_15=abs(Pop15-Full15),
         FULL_TRUE_23=abs(Pop23-Full23),
         FULL_TRUE_24=abs(Pop24-Full24),
         FULL_TRUE_25=abs(Pop25-Full25),
         FULL_TRUE_34=abs(Pop34-Full34),
         FULL_TRUE_35=abs(Pop35-Full35),
         FULL_TRUE_45=abs(Pop45-Full45),
         SFA_TRUE_12=abs(Pop12-SFA12),
         SFA_TRUE_13=abs(Pop13-SFA13),
         SFA_TRUE_14=abs(Pop14-SFA14),
         SFA_TRUE_15=abs(Pop15-SFA15),
         SFA_TRUE_23=abs(Pop23-SFA23),
         SFA_TRUE_24=abs(Pop24-SFA24),
         SFA_TRUE_25=abs(Pop25-SFA25),
         SFA_TRUE_34=abs(Pop34-SFA34),
         SFA_TRUE_35=abs(Pop35-SFA35),
         SFA_TRUE_45=abs(Pop45-SFA45),
         SFB_TRUE_12=abs(Pop12-SFB12),
         SFB_TRUE_13=abs(Pop13-SFB13),
         SFB_TRUE_14=abs(Pop14-SFB14),
         SFB_TRUE_15=abs(Pop15-SFB15),
         SFB_TRUE_23=abs(Pop23-SFB23),
         SFB_TRUE_24=abs(Pop24-SFB24),
         SFB_TRUE_25=abs(Pop25-SFB25),
         SFB_TRUE_34=abs(Pop34-SFB34),
         SFB_TRUE_35=abs(Pop35-SFB35),
         SFB_TRUE_45=abs(Pop45-SFB45),
         SFC_TRUE_12=abs(Pop12-SFC12),
         SFC_TRUE_13=abs(Pop13-SFC13),
         SFC_TRUE_14=abs(Pop14-SFC14),
         SFC_TRUE_15=abs(Pop15-SFC15),
         SFC_TRUE_23=abs(Pop23-SFC23),
         SFC_TRUE_24=abs(Pop24-SFC24),
         SFC_TRUE_25=abs(Pop25-SFC25),
         SFC_TRUE_34=abs(Pop34-SFC34),
         SFC_TRUE_35=abs(Pop35-SFC35),
         SFC_TRUE_45=abs(Pop45-SFC45),
         MI_TRUE_12=abs(Pop12-MI12),
         MI_TRUE_13=abs(Pop13-MI13),
         MI_TRUE_14=abs(Pop14-MI14),
         MI_TRUE_15=abs(Pop15-MI15),
         MI_TRUE_23=abs(Pop23-MI23),
         MI_TRUE_24=abs(Pop24-MI24),
         MI_TRUE_25=abs(Pop25-MI25),
         MI_TRUE_34=abs(Pop34-MI34),
         MI_TRUE_35=abs(Pop35-MI35),
         MI_TRUE_45=abs(Pop45-MI45),
         FULL_TRUE_r_12=Pop12-Full12,
         FULL_TRUE_r_13=Pop13-Full13,
         FULL_TRUE_r_14=Pop14-Full14,
         FULL_TRUE_r_15=Pop15-Full15,
         FULL_TRUE_r_23=Pop23-Full23,
         FULL_TRUE_r_24=Pop24-Full24,
         FULL_TRUE_r_25=Pop25-Full25,
         FULL_TRUE_r_34=Pop34-Full34,
         FULL_TRUE_r_35=Pop35-Full35,
         FULL_TRUE_r_45=Pop45-Full45,
         SFA_TRUE_r_12=Pop12-SFA12,
         SFA_TRUE_r_13=Pop13-SFA13,
         SFA_TRUE_r_14=Pop14-SFA14,
         SFA_TRUE_r_15=Pop15-SFA15,
         SFA_TRUE_r_23=Pop23-SFA23,
         SFA_TRUE_r_24=Pop24-SFA24,
         SFA_TRUE_r_25=Pop25-SFA25,
         SFA_TRUE_r_34=Pop34-SFA34,
         SFA_TRUE_r_35=Pop35-SFA35,
         SFA_TRUE_r_45=Pop45-SFA45,
         SFB_TRUE_r_12=Pop12-SFB12,
         SFB_TRUE_r_13=Pop13-SFB13,
         SFB_TRUE_r_14=Pop14-SFB14,
         SFB_TRUE_r_15=Pop15-SFB15,
         SFB_TRUE_r_23=Pop23-SFB23,
         SFB_TRUE_r_24=Pop24-SFB24,
         SFB_TRUE_r_25=Pop25-SFB25,
         SFB_TRUE_r_34=Pop34-SFB34,
         SFB_TRUE_r_35=Pop35-SFB35,
         SFB_TRUE_r_45=Pop45-SFB45,
         SFC_TRUE_r_12=Pop12-SFC12,
         SFC_TRUE_r_13=Pop13-SFC13,
         SFC_TRUE_r_14=Pop14-SFC14,
         SFC_TRUE_r_15=Pop15-SFC15,
         SFC_TRUE_r_23=Pop23-SFC23,
         SFC_TRUE_r_24=Pop24-SFC24,
         SFC_TRUE_r_25=Pop25-SFC25,
         SFC_TRUE_r_34=Pop34-SFC34,
         SFC_TRUE_r_35=Pop35-SFC35,
         SFC_TRUE_r_45=Pop45-SFC45,
         MI_TRUE_r_12=Pop12-MI12,
         MI_TRUE_r_13=Pop13-MI13,
         MI_TRUE_r_14=Pop14-MI14,
         MI_TRUE_r_15=Pop15-MI15,
         MI_TRUE_r_23=Pop23-MI23,
         MI_TRUE_r_24=Pop24-MI24,
         MI_TRUE_r_25=Pop25-MI25,
         MI_TRUE_r_34=Pop34-MI34,
         MI_TRUE_r_35=Pop35-MI35,
         MI_TRUE_r_45=Pop45-MI45) %>%
  mutate(FULL_TRUE=rowMeans(select(.,FULL_TRUE_12:FULL_TRUE_45),na.rm=T),
         SFA_TRUE=rowMeans(select(.,SFA_TRUE_12:SFA_TRUE_45),na.rm=T),
         SFB_TRUE=rowMeans(select(.,SFB_TRUE_12:SFB_TRUE_45),na.rm=T),
         SFC_TRUE=rowMeans(select(.,SFC_TRUE_12:SFC_TRUE_45),na.rm=T),
         MI_TRUE=rowMeans(select(.,MI_TRUE_12:MI_TRUE_45),na.rm=T), # no direction
         FULL_TRUE_r=rowMeans(select(.,FULL_TRUE_r_12:FULL_TRUE_r_45),na.rm=T),
         SFA_TRUE_r=rowMeans(select(.,SFA_TRUE_r_12:SFA_TRUE_r_45),na.rm=T),
         SFB_TRUE_r=rowMeans(select(.,SFB_TRUE_r_12:SFB_TRUE_r_45),na.rm=T),
         SFC_TRUE_r=rowMeans(select(.,SFC_TRUE_r_12:SFC_TRUE_r_45),na.rm=T),
         MI_TRUE_r=rowMeans(select(.,MI_TRUE_r_12:MI_TRUE_r_45),na.rm=T))

overall_long <- data.frame(Pop=c(cor_res_2$Pop12,cor_res_3$Pop12,cor_res_3$Pop23,cor_res_3$Pop13,cor_res_4$Pop12,cor_res_4$Pop13,cor_res_4$Pop14,cor_res_4$Pop23,cor_res_4$Pop24,cor_res_4$Pop34,
                                     cor_res_5$Pop12,cor_res_5$Pop13,cor_res_5$Pop14,cor_res_5$Pop15,cor_res_5$Pop23,cor_res_5$Pop24,cor_res_5$Pop25,cor_res_5$Pop34,cor_res_5$Pop35,cor_res_5$Pop45),
                           Full=c(cor_res_2$Full12,cor_res_3$Full12,cor_res_3$Full23,cor_res_3$Full13,cor_res_4$Full12,cor_res_4$Full13,cor_res_4$Full14,cor_res_4$Full23,cor_res_4$Full24,cor_res_4$Full34,
                                 cor_res_5$Full12,cor_res_5$Full13,cor_res_5$Full14,cor_res_5$Full15,cor_res_5$Full23,cor_res_5$Full24,cor_res_5$Full25,cor_res_5$Full34,cor_res_5$Full35,cor_res_5$Full45),
                           SFA=c(cor_res_2$SFA12,cor_res_3$SFA12,cor_res_3$SFA23,cor_res_3$SFA13,cor_res_4$SFA12,cor_res_4$SFA13,cor_res_4$SFA14,cor_res_4$SFA23,cor_res_4$SFA24,cor_res_4$SFA34,
                                 cor_res_5$SFA12,cor_res_5$SFA13,cor_res_5$SFA14,cor_res_5$SFA15,cor_res_5$SFA23,cor_res_5$SFA24,cor_res_5$SFA25,cor_res_5$SFA34,cor_res_5$SFA35,cor_res_5$SFA45),
                           SFB=c(cor_res_2$SFB12,cor_res_3$SFB12,cor_res_3$SFB23,cor_res_3$SFB13,cor_res_4$SFB12,cor_res_4$SFB13,cor_res_4$SFB14,cor_res_4$SFB23,cor_res_4$SFB24,cor_res_4$SFB34,
                                 cor_res_5$SFB12,cor_res_5$SFB13,cor_res_5$SFB14,cor_res_5$SFB15,cor_res_5$SFB23,cor_res_5$SFB24,cor_res_5$SFB25,cor_res_5$SFB34,cor_res_5$SFB35,cor_res_5$SFB45),
                           SFC=c(cor_res_2$SFC12,cor_res_3$SFC12,cor_res_3$SFC23,cor_res_3$SFC13,cor_res_4$SFC12,cor_res_4$SFC13,cor_res_4$SFC14,cor_res_4$SFC23,cor_res_4$SFC24,cor_res_4$SFC34,
                                 cor_res_5$SFC12,cor_res_5$SFC13,cor_res_5$SFC14,cor_res_5$SFC15,cor_res_5$SFC23,cor_res_5$SFC24,cor_res_5$SFC25,cor_res_5$SFC34,cor_res_5$SFC35,cor_res_5$SFC45),
                           MI=c(cor_res_2$MI12,cor_res_3$MI12,cor_res_3$MI23,cor_res_3$MI13,cor_res_4$MI12,cor_res_4$MI13,cor_res_4$MI14,cor_res_4$MI23,cor_res_4$MI24,cor_res_4$MI34,
                                 cor_res_5$MI12,cor_res_5$MI13,cor_res_5$MI14,cor_res_5$MI15,cor_res_5$MI23,cor_res_5$MI24,cor_res_5$MI25,cor_res_5$MI34,cor_res_5$MI35,cor_res_5$MI45))
cor(overall_long,use="pairwise.complete.obs")

cor_res_2_5 <- select(cor_res_2,iteration,PMissing,SampleSize,NFac,PopR_mean,PopR_sd,FULL_TRUE,SFA_TRUE,SFB_TRUE,SFC_TRUE,MI_TRUE,FULL_TRUE_r,SFA_TRUE_r,SFB_TRUE_r,SFC_TRUE_r,MI_TRUE_r) %>%
  rbind(select(cor_res_3,iteration,PMissing,SampleSize,NFac,PopR_mean,PopR_sd,FULL_TRUE,SFA_TRUE,SFB_TRUE,SFC_TRUE,MI_TRUE,FULL_TRUE_r,SFA_TRUE_r,SFB_TRUE_r,SFC_TRUE_r,MI_TRUE_r)) %>%
  rbind(select(cor_res_4,iteration,PMissing,SampleSize,NFac,PopR_mean,PopR_sd,FULL_TRUE,SFA_TRUE,SFB_TRUE,SFC_TRUE,MI_TRUE,FULL_TRUE_r,SFA_TRUE_r,SFB_TRUE_r,SFC_TRUE_r,MI_TRUE_r)) %>%
  rbind(select(cor_res_5,iteration,PMissing,SampleSize,NFac,PopR_mean,PopR_sd,FULL_TRUE,SFA_TRUE,SFB_TRUE,SFC_TRUE,MI_TRUE,FULL_TRUE_r,SFA_TRUE_r,SFB_TRUE_r,SFC_TRUE_r,MI_TRUE_r))

overall_des <- cor_res_2_5 %>%
  summarize(FULL_TRUE=paste0(mean(.$FULL_TRUE,na.rm=T)," (",sd(.$FULL_TRUE,na.rm=T),")"),
            SFA_TRUE=paste0(mean(.$SFA_TRUE,na.rm=T)," (",sd(.$SFA_TRUE,na.rm=T),")"),
            SFB_TRUE=paste0(mean(.$SFB_TRUE,na.rm=T)," (",sd(.$SFB_TRUE,na.rm=T),")"),
            SFC_TRUE=paste0(mean(.$SFC_TRUE,na.rm=T)," (",sd(.$SFC_TRUE,na.rm=T),")"),
            MI_TRUE=paste0(mean(.$MI_TRUE,na.rm=T)," (",sd(.$MI_TRUE,na.rm=T),")"),
            FULL_TRUE_r=paste0(mean(.$FULL_TRUE_r,na.rm=T)," (",sd(.$FULL_TRUE_r,na.rm=T),")"),
            SFA_TRUE_r=paste0(mean(.$SFA_TRUE_r,na.rm=T)," (",sd(.$SFA_TRUE_r,na.rm=T),")"),
            SFB_TRUE_r=paste0(mean(.$SFB_TRUE_r,na.rm=T)," (",sd(.$SFB_TRUE_r,na.rm=T),")"),
            SFC_TRUE_r=paste0(mean(.$SFC_TRUE_r,na.rm=T)," (",sd(.$SFC_TRUE_r,na.rm=T),")"),
            MI_TRUE_r=paste0(mean(.$MI_TRUE_r,na.rm=T)," (",sd(.$MI_TRUE_r,na.rm=T),")"))

cor_res_all <- cbind(cor_res_2_5 %>%
                       select(iteration:MI_TRUE) %>%
                       pivot_longer(cols=FULL_TRUE:MI_TRUE,names_to="technique",values_to="abs_dev"),
                     cor_res_2_5 %>%
                       select(iteration:PopR_sd,FULL_TRUE_r:MI_TRUE_r) %>%
                       pivot_longer(cols=FULL_TRUE_r:MI_TRUE_r,names_to="technique",values_to="raw_dev") %>%
                       select(raw_dev))

# contrasts
cor_res_all <- within(cor_res_all,technique <- relevel(as.factor(technique),ref="FULL_TRUE"))
levels(cor_res_all$technique)
c1 <- c(-4,1,1,1,1)
c2 <- c(0,1,1,1,-3)
c3 <- c(0,1,-2,1,0)
c4 <- c(0,1,0,-1,0)
contrasts(cor_res_all$technique) <- cbind(c1,c2,c3,c4)
contrasts(cor_res_all$technique)

# regression-absolute deviations
reg_abs1 <- lm(abs_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique,data=cor_res_all)
reg_abs1_std <- lm(log(abs_dev)~scale(PMissing,scale=T)+scale(SampleSize/100,scale=T)+scale(NFac,scale=T)+scale(PopR_mean,scale=T)+scale(PopR_sd,scale=T)+technique,data=cor_res_all)

reg_abs2 <- lm(abs_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique+
                 scale(PMissing,scale=F)*technique+scale(SampleSize/100,scale=F)*technique+scale(NFac,scale=F)*technique+scale(PopR_mean,scale=F)*technique+scale(PopR_sd,scale=F)*technique,data=cor_res_all)
reg_abs2_std <- lm(abs_dev~scale(PMissing,scale=T)+scale(SampleSize/100,scale=T)+scale(NFac,scale=T)+scale(PopR_mean,scale=T)+scale(PopR_sd,scale=T)+technique+
                 scale(PMissing,scale=T)*technique+scale(SampleSize/100,scale=T)*technique+scale(NFac,scale=T)*technique+scale(PopR_mean,scale=T)*technique+scale(PopR_sd,scale=T)*technique,data=cor_res_all)

reg_abs3 <- lm(abs_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique+
                 scale(PMissing,scale=F)*technique+scale(SampleSize/100,scale=F)*technique+scale(NFac,scale=F)*technique+scale(PopR_mean,scale=F)*technique+scale(PopR_sd,scale=F)*technique+
                 scale(PMissing,scale=F)*scale(SampleSize/100,scale=F)+scale(PMissing,scale=F)*scale(NFac,scale=F)+scale(PMissing,scale=F)*scale(PopR_mean,scale=F)+scale(PMissing,scale=F)*scale(PopR_sd,scale=F)+
                 scale(SampleSize/100,scale=F)*scale(NFac,scale=F)+scale(SampleSize/100,scale=F)*scale(PopR_mean,scale=F)+scale(SampleSize/100,scale=F)*scale(PopR_sd,scale=F)+
                 scale(NFac,scale=F)*scale(PopR_mean,scale=F)+scale(NFac,scale=F)*scale(PopR_sd,scale=F)+scale(PopR_mean,scale=F)*scale(PopR_sd,scale=F),data=cor_res_all)
# par(mfrow = c(2, 2))
# plot(reg_abs3) # saved in distribution_diagnostics

reg_abs3_std <- lm(abs_dev~scale(PMissing,scale=T)+scale(SampleSize/100,scale=T)+scale(NFac,scale=T)+scale(PopR_mean,scale=T)+scale(PopR_sd,scale=T)+technique+
                 scale(PMissing,scale=T)*technique+scale(SampleSize/100,scale=T)*technique+scale(NFac,scale=T)*technique+scale(PopR_mean,scale=T)*technique+scale(PopR_sd,scale=T)*technique+
                 scale(PMissing,scale=T)*scale(SampleSize/100,scale=T)+scale(PMissing,scale=T)*scale(NFac,scale=T)+scale(PMissing,scale=T)*scale(PopR_mean,scale=T)+scale(PMissing,scale=T)*scale(PopR_sd,scale=T)+
                 scale(SampleSize/100,scale=T)*scale(NFac,scale=T)+scale(SampleSize/100,scale=T)*scale(PopR_mean,scale=T)+scale(SampleSize/100,scale=T)*scale(PopR_sd,scale=T)+
                 scale(NFac,scale=T)*scale(PopR_mean,scale=T)+scale(NFac,scale=T)*scale(PopR_sd,scale=T)+scale(PopR_mean,scale=T)*scale(PopR_sd,scale=T),data=cor_res_all)
# reg_abs3_g <- glm(abs_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique+
#                  scale(PMissing,scale=F)*technique+scale(SampleSize/100,scale=F)*technique+scale(NFac,scale=F)*technique+scale(PopR_mean,scale=F)*technique+scale(PopR_sd,scale=F)*technique+
#                  scale(PMissing,scale=F)*scale(SampleSize/100,scale=F)+scale(PMissing,scale=F)*scale(NFac,scale=F)+scale(PMissing,scale=F)*scale(PopR_mean,scale=F)+scale(PMissing,scale=F)*scale(PopR_sd,scale=F)+
#                  scale(SampleSize/100,scale=F)*scale(NFac,scale=F)+scale(SampleSize/100,scale=F)*scale(PopR_mean,scale=F)+scale(SampleSize/100,scale=F)*scale(PopR_sd,scale=F)+
#                  scale(NFac,scale=F)*scale(PopR_mean,scale=F)+scale(NFac,scale=F)*scale(PopR_sd,scale=F)+scale(PopR_mean,scale=F)*scale(PopR_sd,scale=F),data=cor_res_all,family=Gamma(link="log")) # tried gamme distribution with log link
# par(mfrow = c(2, 2))
# plot(reg_abs3_g) # saved in distribution_diagnostics

abs1 <- cbind(coef(reg_abs1),coef(reg_abs1_std)) 
abs1 <- rbind(abs1,c(summary(reg_abs1)$adj.r.squared,NA))

abs2 <- cbind(coef(reg_abs2),coef(reg_abs2_std)) 
abs2 <- rbind(abs2,c(summary(reg_abs2)$adj.r.squared,NA))

abs3 <- cbind(coef(reg_abs3),coef(reg_abs3_std)) 
abs3 <- rbind(abs3,c(summary(reg_abs3)$adj.r.squared,NA))

wb <- createWorkbook()
addWorksheet(wb,"absolute deviations")
writeData(wb,1,abs1,startCol=1,rowNames=T)
writeData(wb,1,abs2,startCol=5,rowNames=T)
writeData(wb,1,abs3,startCol=9,rowNames=T)
saveWorkbook(wb,here::here("analysis_results","table1.xlsx"),overwrite=T)


# predicted values
test_dat <- data.frame(PMissing=rep(c(.2,.4,.6,.8),20),
                       SampleSize=rep(c(200,400,600,800,1000),each=16),
                       NFac=5,
                       PopR_mean=rep(rep(c(.2,.4,.6,.8),each=4),5),
                       PopR_sd=.11) #mean
predicted <- test_dat %>%
  bind_cols(Full=predict(object=reg_abs3,
                        newdata=cbind(test_dat,data.frame(technique="FULL_TRUE")))) %>%
  bind_cols(SFA=predict(object=reg_abs3,
                        newdata=cbind(test_dat,data.frame(technique="SFA_TRUE")))) %>%
  bind_cols(SFB=predict(object=reg_abs3,
                        newdata=cbind(test_dat,data.frame(technique="SFB_TRUE")))) %>%
  bind_cols(SFC=predict(object=reg_abs3,
                        newdata=cbind(test_dat,data.frame(technique="SFC_TRUE")))) %>%
  bind_cols(PM=predict(object=reg_abs3,
                        newdata=cbind(test_dat,data.frame(technique="MI_TRUE")))) %>%
  select(-NFac,-PopR_sd)
predicted[,4:8] <- round(predicted[,4:8],3)
write_csv(predicted,here::here("analysis_results","table3.csv"))


# regression-raw deviations
reg_raw1 <- lm(raw_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique,data=cor_res_all)
reg_raw1_std <- lm(raw_dev~scale(PMissing,scale=T)+scale(SampleSize/100,scale=T)+scale(NFac,scale=T)+scale(PopR_mean,scale=T)+scale(PopR_sd,scale=T)+technique,data=cor_res_all)
reg_raw2 <- lm(raw_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique+
                 scale(PMissing,scale=F)*technique+scale(SampleSize/100,scale=F)*technique+scale(NFac,scale=F)*technique+scale(PopR_mean,scale=F)*technique+scale(PopR_sd,scale=F)*technique,data=cor_res_all)
reg_raw2_std <- lm(raw_dev~scale(PMissing,scale=T)+scale(SampleSize/100,scale=T)+scale(NFac,scale=T)+scale(PopR_mean,scale=T)+scale(PopR_sd,scale=T)+technique+
                 scale(PMissing,scale=T)*technique+scale(SampleSize/100,scale=T)*technique+scale(NFac,scale=T)*technique+scale(PopR_mean,scale=T)*technique+scale(PopR_sd,scale=T)*technique,data=cor_res_all)
reg_raw3 <- lm(raw_dev~scale(PMissing,scale=F)+scale(SampleSize/100,scale=F)+scale(NFac,scale=F)+scale(PopR_mean,scale=F)+scale(PopR_sd,scale=F)+technique+
                 scale(PMissing,scale=F)*technique+scale(SampleSize/100,scale=F)*technique+scale(NFac,scale=F)*technique+scale(PopR_mean,scale=F)*technique+scale(PopR_sd,scale=F)*technique+
                 scale(PMissing,scale=F)*scale(SampleSize/100,scale=F)+scale(PMissing,scale=F)*scale(NFac,scale=F)+scale(PMissing,scale=F)*scale(PopR_mean,scale=F)+scale(PMissing,scale=F)*scale(PopR_sd,scale=F)+
                 scale(SampleSize/100,scale=F)*scale(NFac,scale=F)+scale(SampleSize/100,scale=F)*scale(PopR_mean,scale=F)+scale(SampleSize/100,scale=F)*scale(PopR_sd,scale=F)+
                 scale(NFac,scale=F)*scale(PopR_mean,scale=F)+scale(NFac,scale=F)*scale(PopR_sd,scale=F)+scale(PopR_mean,scale=F)*scale(PopR_sd,scale=F),data=cor_res_all)
# par(mfrow = c(2, 2))
# plot(reg_raw3) # saved in distribution_diagnostics
reg_raw3_std <- lm(raw_dev~scale(PMissing,scale=T)+scale(SampleSize/100,scale=T)+scale(NFac,scale=T)+scale(PopR_mean,scale=T)+scale(PopR_sd,scale=T)+technique+
                 scale(PMissing,scale=T)*technique+scale(SampleSize/100,scale=T)*technique+scale(NFac,scale=T)*technique+scale(PopR_mean,scale=T)*technique+scale(PopR_sd,scale=T)*technique+
                 scale(PMissing,scale=T)*scale(SampleSize/100,scale=T)+scale(PMissing,scale=T)*scale(NFac,scale=T)+scale(PMissing,scale=T)*scale(PopR_mean,scale=T)+scale(PMissing,scale=T)*scale(PopR_sd,scale=T)+
                 scale(SampleSize/100,scale=T)*scale(NFac,scale=T)+scale(SampleSize/100,scale=T)*scale(PopR_mean,scale=T)+scale(SampleSize/100,scale=T)*scale(PopR_sd,scale=T)+
                 scale(NFac,scale=T)*scale(PopR_mean,scale=T)+scale(NFac,scale=T)*scale(PopR_sd,scale=T)+scale(PopR_mean,scale=T)*scale(PopR_sd,scale=T),data=cor_res_all)

raw1 <- cbind(coef(reg_raw1),coef(reg_raw1_std)) 
raw1 <- rbind(raw1,c(summary(reg_raw1)$adj.r.squared,NA))

raw2 <- cbind(coef(reg_raw2),coef(reg_raw2_std)) 
raw2 <- rbind(raw2,c(summary(reg_raw2)$adj.r.squared,NA))

raw3 <- cbind(coef(reg_raw3),coef(reg_raw3_std)) 
raw3 <- rbind(raw3,c(summary(reg_raw3)$adj.r.squared,NA))

wb <- createWorkbook()
addWorksheet(wb,"raw deviations")
writeData(wb,1,raw1,startCol=1,rowNames=T)
writeData(wb,1,raw2,startCol=5,rowNames=T)
writeData(wb,1,raw3,startCol=9,rowNames=T)
saveWorkbook(wb,here::here("analysis_results","table2.xlsx"),overwrite=T)


# figures
# figure 1 is generated manually
fig_2 <- ggplot(cor_res_all,aes(x=factor(plyr::round_any(SampleSize,100)),y=abs_dev,fill=ordered(technique,c("FULL_TRUE","SFA_TRUE","SFB_TRUE","SFC_TRUE","MI_TRUE")))) + 
  geom_boxplot(outlier.shape=NA) + 
  xlab("Sample Size") +
  ylab("Absolute deviation from true correlations") +
  ylim(0,.3) +
  theme_bw() +
  theme(text=element_text(size=24),
        legend.position="right") +
  scale_fill_manual(values=c("White",adjustcolor("#7a0019",alpha.f=.4),adjustcolor("#7a0019",alpha.f=.7),"#7a0019","Black"),
                    name="Method",labels=c("Full length","Short Form A","Short Form B","Short Form C","Planned Missingness")) 
ggsave(here::here("analysis_results","fig2.png"),fig_2,width=15,height=8,units="in",dpi=300)

fig_3 <- ggplot(cor_res_all,aes(x=factor(PMissing),y=abs_dev,fill=ordered(technique,c("FULL_TRUE","SFA_TRUE","SFB_TRUE","SFC_TRUE","MI_TRUE")))) + 
  geom_boxplot(outlier.shape=NA) + 
  xlab("Proportion of missingness") +
  ylab("Absolute deviation from true correlations") +
  ylim(0,.35) +
  theme_bw() +
  theme(text=element_text(size=24),
        legend.position="right") +
  scale_fill_manual(values=c("White",adjustcolor("#7a0019",alpha.f=.4),adjustcolor("#7a0019",alpha.f=.7),"#7a0019","Black"),
                    name="Method",labels=c("Full length","Short Form A","Short Form B","Short Form C","Planned Missingness")) 
ggsave(here::here("analysis_results","fig3.png"),fig_3,width=15,height=8,units="in",dpi=300)

fig_4 <- ggplot(cor_res_all,aes(x=factor(plyr::round_any(PopR_mean,.1)),y=abs_dev,fill=ordered(technique,c("FULL_TRUE","SFA_TRUE","SFB_TRUE","SFC_TRUE","MI_TRUE")))) + 
  geom_boxplot(outlier.shape=NA) + 
  xlab("Mean population intercorrelation") +
  ylab("Absolute deviation from true correlations") +
  ylim(0,.35) +
  theme_bw() +
  theme(text=element_text(size=24),
        legend.position="right") +
  scale_fill_manual(values=c("White",adjustcolor("#7a0019",alpha.f=.4),adjustcolor("#7a0019",alpha.f=.7),"#7a0019","Black"),
                    name="Method",labels=c("Full length","Short Form A","Short Form B","Short Form C","Planned Missingness")) 
ggsave(here::here("analysis_results","fig4.png"),fig_4,width=15,height=8,units="in",dpi=300)

fig_5 <- ggplot(cor_res_all,aes(x=factor(plyr::round_any(PopR_mean,.1)),y=abs_dev,fill=ordered(technique,c("FULL_TRUE","SFA_TRUE","SFB_TRUE","SFC_TRUE","MI_TRUE")))) + 
  geom_boxplot(outlier.shape=NA) + 
  xlab("Mean population intercorrelation") +
  ylab("Absolute deviation from true correlations") +
  ylim(0,.45) +
  theme_bw() +
  theme(text=element_text(size=24),
        legend.position="right") +
  scale_fill_manual(values=c("White",adjustcolor("#7a0019",alpha.f=.4),adjustcolor("#7a0019",alpha.f=.7),"#7a0019","Black"),
                    name="Method",labels=c("Full length","Short Form A","Short Form B","Short Form C","Planned Missingness")) +
  facet_wrap(~factor(round_any(PMissing,.1)),ncol=4,labeller=as_labeller(c(`0.1`="10% Missing",
                                                                           `0.2`="20% Missing",
                                                                           `0.3`="30% Missing",
                                                                           `0.4`="40% Missing",
                                                                           `0.5`="50% Missing",
                                                                           `0.6`="60% Missing",
                                                                           `0.7`="70% Missing",
                                                                           `0.8`="80% Missing"))) 
ggsave(here::here("analysis_results","fig5.png"),fig_5,width=18,height=8,units="in",dpi=300)

failure <- cor_res_2_5 %>%
  group_by(PMissing,SampleSize,NFac) %>%
  dplyr::summarize(total=n(),
                   failure=sum(is.na(MI_TRUE)),
                   convergence=(total-failure)/total*100)
failure %>%
  group_by(PMissing) %>%
  summarize(failure=sum(failure)) %>%
  ungroup()
failure %>%
  group_by(SampleSize) %>%
  summarize(failure=sum(failure)) %>%
  ungroup()

fig_6 <- ggplot(filter(failure,PMissing>.6),aes(x=SampleSize,y=convergence,group=factor(round_any(PMissing,.1)))) +
  geom_line(aes(color=factor(round_any(PMissing,.1)))) +
  geom_point(aes(color=factor(round_any(PMissing,.1)))) +
  scale_color_manual(values=c(adjustcolor("#7a0019",alpha.f=.3),adjustcolor("#7a0019",alpha.f=.8)),
                     name="Proportion Missing per Scale",labels=c(".7",".8")) +
  xlab("Sample Size") +
  ylab("Convergence Rate (%)") +
  theme_bw() +
  theme(text=element_text(size=12),
        legend.position="bottom") +
  facet_wrap(~NFac) +
  facet_grid(~NFac,labeller=as_labeller(c(`2`="2 Factors",
                                          `3`="3 Factors",
                                          `4`="4 Factors",
                                          `5`="5 Factors"))) 
ggsave(here::here("analysis_results","fig6.png"),fig_6,width=12,height=8,units="in",dpi=300)

