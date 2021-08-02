####################################################################
# Author: QY YU
# Date: 19 July 2021
# Purpose: Load required R packages or functions for data analyses
####################################################################


# 1. Load R packages ------------------------------------------------------
# install.packages('pacman')
pacman::p_load(data.table, tidyverse, rio, lubridate, reshape2, lmtest, stringi, sandwich, dbplyr, parallel, RSQLite, magrittr)


# 2. Function to combine raw data files -----------------------------------
mergefile <- function(path){
  filenames = list.files(path = path, pattern = '*.xlsx', full.names = T)
  datalist = lapply(filenames, readxl::read_excel)
  Reduce(function(x,y){merge(x,y, all = T)}, datalist)
}


# 3. Function to rename column names --------------------------------------
colrename<- function(x){
  tmp <- make.names(tolower(colnames(x)), unique = T)
  tmp1 <- gsub(pattern = "\\.+", replacement = "\\.", tmp)
  newcolnames <<- gsub(pattern = "\\.+$", replacement = "", tmp1)
  return(newcolnames)
}


# 4. Function to convert dob to age ---------------------------------------
age <- function(dob, age.day, units = "years", floor = TRUE) {
  calc.age = interval(dob, age.day) / duration(num = 1, units = units)
  if (floor) return(as.integer(floor(calc.age)))
  return(calc.age)
}


# 5. Load HK census population from 2001 to 2020 --------------------------
hk_pop_agegp <- readRDS('../8. HK_Population/hk_pop_agegp.rds')


# -------------------------------------------------------------------------



###############################################
# End of general R packages and functions.
# Final modified date on July 19, 2021.
###############################################




