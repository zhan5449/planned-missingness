library(tidyverse)
##take output of simulation2_diagnostics.R to create plots
#import files i1-iX in Two-Factor Diagnostics
#each has iterations for 20,40,60,80% missing (all same pop cor-.4 and n-500)
#desired output:
#one plot per parameter- AC, psrf, Percent Bias, coverage rate, fmi
#x axis should be iterations 
#y is parameter
#4 lines for each level of missingness

##BEFORE RUNNING CODE BELOW, SET GLOBAL GGPLOT THEME PARAMETERS
#doing this globally so i don't need to repeat call to theme() and set same format for each plot
#open .Rprofile and type in code below from .First to end of this function
#then restart R session for theme to take effect
#see https://ggplot2tutor.com/tutorials/apa_theme_rprofile for steps
# .First <- function() {
#   # Load packages
#   library(tidyverse)
#   
#   # Set theme in .Rprofile- comment out when done
#   #usethis::edit_r_profile()
#   apa_theme <- theme(
#     plot.margin = unit(c(1, 1, 1, 1), "cm"),
#     plot.background = element_rect(fill = "white", color = NA),
#     plot.title = element_text(size = 22, face = "bold",
#                               hjust = 0.5,
#                               margin = margin(b = 15)),
#     axis.line = element_line(color = "black", size = .5),
#     axis.title = element_text(size = 18, color = "black",
#                               face = "bold"),
#     axis.text = element_text(size = 15, color = "black"),
#     axis.text.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10)),
#     axis.ticks = element_line(size = .5),
#     panel.grid = element_blank(),
#     legend.position = "bottom",
#     legend.background = element_rect(color = "black"),
#     legend.text = element_text(size = 15),
#     legend.margin = margin(t = 5, l = 5, r = 5, b = 5),
#     legend.key = element_rect(color = NA, fill = NA)
#   )
#   
#   theme_set(theme_minimal(base_size = 18) +
#               apa_theme)
# }

#the different i1-iX are just simulations of the same task with different seeds
#i.e., i want to average across iteration=3 from i1-iX and this will be one point on the plot
#should i compute SD across simulations? ignore for now for time

###CODE BELOW IS HASHED OUT AFTER FIRST RUN THROUGH TO IMPORT AND SAVE COMBINED FILE THAT AVERAGES SIMULATION RESULTS
##step 1- read in all files
#setwd("./Two-Factor Diagnostics/")
# temp <- list.files(pattern=".csv")
# sim_results <- lapply(temp, read_csv) #will this be a list object
# 
# ##step 2- collapse columns across each list object- take average for numeric columns
# #data.table::rbindlist(sim_results) #do i need this, i don't necessarily want one large large df
# 
# sim_averages<-bind_rows(sim_results) %>% 
#   group_by(no_missing,iteration_no) %>% 
#   summarise(across(c(cor_pooled:psrf), mean, .names="{.col}")) %>% 
#   mutate(cor_pop=0.398302,
#          sample_size=500) %>%  #P6
#   mutate(missing_prop=case_when(
#     no_missing == 2 ~ "20%",
#     no_missing == 4 ~ "40%",
#     no_missing == 6 ~ "60%",
#     no_missing == 8 ~ "80%"
#   ))
# write_rds(sim_averages,"./sim_averages.rds")

sim_averages <- read_rds("./Two-Factor Diagnostics/sim_averages.rds")

##step 3 - make plots *Note, my diagnostics code runs more iterations for 60% and 80% so x axis will be shorter for the 20% and 40%

###PERFORMANCE ESTIMATES
##Cor_pooled - reference the pop_cor
cor_pooled <- sim_averages %>% select(missing_prop, cor_pop, iteration_no, cor_pooled)
cor_pooled %>% ggplot(aes(x=iteration_no,y=cor_pooled, color=missing_prop)) +
  geom_line() +
  geom_line(aes(y=cor_pop),linetype="dashed",color="grey") +
  labs(x="Iteration",y="Imputed Factor Correlation",color="Missing Rate (%)") +
  ylim(c(0.3,0.5))
ggsave("cor_pooled.png")
##pb
PB <- sim_averages %>% select(missing_prop, cor_pop, iteration_no, PB)
PB %>% ggplot(aes(x=iteration_no,y=PB, color=missing_prop)) +
  geom_line() +
  geom_line(aes(y=5), linetype="dashed",color="grey") +
  geom_line(aes(y=0), linetype="dashed",color="grey") +
  scale_fill_grey() +
  # jtools::theme_apa() +
  labs(x="Iteration",y="Percent Bias(%)",color="Missing Rate(%)")
ylim(c(0.0,13))
ggsave("PB.png")

##fmi
fmi <- sim_averages %>% select(missing_prop, cor_pop, iteration_no, fmi)
mn_by_miss <- fmi %>% 
  group_by(missing_prop) %>% summarise(mn_fmi=mean(fmi))
check <- fmi %>% ungroup() %>% 
  filter(no_missing ==8 )
#THRESHOLD AT WHICH FMI/m = .01
#for 20%, that is .2, .4, .6, .8 etc. values here and lower good
fmi %>% ggplot(aes(x=iteration_no,y=fmi, group=missing_prop)) +
  geom_line(aes(color=missing_prop),linewidth=.6) + 
  scale_fill_grey() +
  labs(x="Iteration",y="Fraction of Missing Information",color="Missing Rate(%)")



### CONVERGENCE ESTIMATES
##AC
ac <- sim_averages %>% select(missing_prop, cor_pop, iteration_no, ac)
ac %>% ggplot(aes(x=iteration_no,y=ac, color=missing_prop)) +
  geom_line() +
  geom_line(aes(y=0),linetype="dashed",color="grey") +
  labs(x="Iteration",y="Autocorrelation",color="Missing Rate (%)")+
ylim(c(0,1))
ggsave("AC.png")

##psrf
psrf <- sim_averages %>% select(missing_prop, cor_pop, iteration_no, psrf)
psrf %>% ggplot(aes(x=iteration_no,y=psrf, color=missing_prop)) +
  geom_line() +
  geom_line(aes(y=1.01),linetype="dashed",color="grey") +
  geom_line(aes(y=1.10),linetype="dashed",color="grey") +
  labs(x="Iteration",y=expression(hat(R)),color="Missing Rate (%)") +
  ylim(c(1.0,1.8))
ggsave("psrf.png")


  


