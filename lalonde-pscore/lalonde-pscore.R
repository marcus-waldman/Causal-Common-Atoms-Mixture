rm(list = ls())

library(tidyverse)
library(MplusAutomation)

lalonde = read.table("data/Asst4Data_2014_Obs.txt", header=T, quote="\"") %>% 
  dplyr::relocate(RE78:HISPANIC,U74,U75)
names(lalonde) = tolower(names(lalonde))


X = lalonde %>% dplyr::select(educ:re75) %>% dplyr::mutate(across(everything(), scale))
names(X) = paste0("x",1:4)
#QR = qr(X)  
#Q = qr.Q(QR) %>% as.data.frame()
#names(Q) = paste0("q",1:ncol(Q))
#R = qr.R(QR)

dat = lalonde %>% dplyr::select(re78:hispanic) %>% dplyr::bind_cols(X)
MplusAutomation::prepareMplusData(df = dat, filename = "data/lalonde_obs.csv")
