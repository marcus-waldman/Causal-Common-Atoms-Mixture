rm(list = ls())

library(tidyverse)
library(MplusAutomation)

lalonde = read.table("data/Asst4Data_2014_Obs.txt", header=T, quote="\"") %>% 
  dplyr::relocate(RE78:HISPANIC,U74,U75)
names(lalonde) = tolower(names(lalonde))


X = lalonde %>% dplyr::select(educ:re75, -u74, -u75) %>% 
  dplyr::mutate(re74 = log(re74 + 1), re75 = log(re75+1)) %>% 
  dplyr::mutate(across(educ:age, scale)) %>% 
  dplyr::mutate(across(re74:re75, function(x){scale(x, center = F, scale = T)}))

#QR = qr(X)  
#Q = qr.Q(QR) %>% as.data.frame()
#names(Q) = paste0("q",1:ncol(Q))
#R = qr.R(QR)

# Let's remove observations without common support
dat = lalonde %>% dplyr::select(re78, treat) %>% dplyr::bind_cols(X)


fit = glm(treat~., data = dat %>% dplyr::select(-re78, -re74, -age), family = "binomial")
dat = dat %>% dplyr::mutate(e = qlogis(fit$fitted.values), e = scale(e))

ggplot(dat, aes(e)) + 
  geom_histogram() +
  facet_grid(treat~., scale="free_y")
  
# Indicator for outside of common support
e1_min = with(dat %>% dplyr::filter(treat==1), min(e))
e1_max = with(dat %>% dplyr::filter(treat==1), max(e))

dat = dat %>% 
  dplyr::mutate(lout = as.integer(e<e1_min), rout = as.integer(e>e1_max)) %>% 
  dplyr::mutate(u = ifelse(lout==1,0,1), 
                u = ifelse(rout==1,2,u))


mplus_dat = dat %>% dplyr::select(re78, treat, u, e, educ, re75)

# let's get the minimum
with(mplus_dat %>% dplyr::filter(e<e1_min), mean(e))

# # Quantiles
# mplus_dat = mplus_dat %>% 
#   dplyr::mutate(u75 = ifelse(re75==0,1,0), 
#                 re75 = ifelse(u75==1,NA,re75)
#   )

MplusAutomation::prepareMplusData(df = mplus_dat, filename = "data/lalonde_obs.dat")
