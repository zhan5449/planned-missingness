FROM rocker/verse:4.4.3
RUN install2.r --error --skipinstalled\
  here tidyverse MASS psych rlist mice miceadds doParallel parallel doRNG openxlsx matrixStats sjPlot wesanderson ggplot2 broom ggfortify DescTools plyr
WORKDIR /home/rstudio
