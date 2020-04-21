# hw4_data_clean.R
# Noah Stafford

# libraries
library(tidyverse)
library(readxl)

# Setwd
setwd("C:/Users/woody/Google Drive/github_repos/826_Advanced_Data_Analysis/hw4")

# Clean data and save to .RData file
nh_diab <- read.csv('./rawData/nhanes_diabetes.csv', stringsAsFactors = FALSE)
