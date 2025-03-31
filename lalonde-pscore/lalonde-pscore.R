rm(list = ls())

library(tidyverse)
library(MplusAutomation)
library(Matching)
library(ggthemes)


data(lalonde)


glm1 <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
              hisp + married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
              u74 + u75, family=binomial, data=lalonde)

lalonde = lalonde %>% dplyr::mutate(e = glm1$fitted)


breaks = quantile(lalonde$e, probs = seq(0,1,len = 6))


lalonde = lalonde %>% dplyr::mutate(strata = cut(e, breaks = breaks, labels = LETTERS[1:5], include.lowest = T))


nas = lalonde %>% 
  dplyr::group_by(strata,treat) %>% 
  dplyr::summarise(n_as = n())

ns = lalonde  %>% dplyr::group_by(strata) %>% 
  dplyr::summarise(n_s = n())

N = nrow(lalonde)

na = lalonde %>% dplyr::group_by(treat) %>% 
  dplyr::summarize(n_a = n())

nas = nas %>% dplyr::left_join(ns) %>% dplyr::left_join(na) %>% 
  dplyr::mutate(mmws = (n_s/N) * (n_a/n_as) )


lalonde = lalonde %>% dplyr::left_join(nas %>%  dplyr::select(strata,treat, n_a, mmws))

lalonde = lalonde %>% dplyr::mutate(iptw = (n_a/N)/(treat+1-e), iptw = iptw/mean(iptw))


ggplot(lalonde %>% dplyr::mutate(Arm = ifelse(treat==1,"A=1","A=0")), ) +
  geom_density(aes(e, col = Arm, fill = Arm), alpha = .3, show.legend = F) +
  geom_point(aes(x=e, y = 0*e, col = Arm, fill = Arm, size = iptw), shape = 124, show.legend = F) +
  facet_grid(Arm~.) +
  theme_minimal() +
  labs(x = "Estimated Propensity Score", y = "") +
  coord_cartesian(xlim = c(0,1)) +
  ggthemes::scale_colour_colorblind() + ggthemes::scale_fill_colorblind()



ggplot(lalonde %>% dplyr::mutate(Arm = ifelse(treat==1,"A=1","A=0")), ) +
  geom_histogram(aes(e, y=..density.., col = Arm, fill = Arm), breaks = breaks, alpha = .3, show.legend = F) +
  geom_point(aes(x=e, y = 0*e, col = Arm, fill = Arm, size = mmws), shape = 124, show.legend = F) +
  facet_grid(Arm~.) +
  theme_minimal() +
  labs(x = "Estimated Propensity Score", y = "") +
  coord_cartesian(xlim = c(0,1)) +
  ggthemes::scale_colour_colorblind() + 
  ggthemes::scale_fill_colorblind() +
  geom_vline(xintercept = breaks[1], linetype = 2) +
  geom_vline(xintercept = breaks[2], linetype = 2) +
  geom_vline(xintercept = breaks[3], linetype = 2) +
  geom_vline(xintercept = breaks[4], linetype = 2) +
  geom_vline(xintercept = breaks[5], linetype = 2) +
  geom_vline(xintercept = breaks[6], linetype = 2)




ggplot(lalonde %>% tidyr::pivot_longer(mmws:iptw, names_to = "estimator", values_to = "w")) +
  geom_histogram(aes(w, y = ..density..)) +
  facet_grid(estimator~.) +
  theme_bw() + 
  labs(y = "", x = "Weight Estimate")


### Make mplus
X = model.matrix(~0 + poly(age,2) + poly(educ,2) + black +
                   hisp + married + nodegr + poly(re74,2) + poly(re75,2) +
                   u74 + u75, family=binomial, 
                 data=lalonde %>% dplyr::select(age,educ,re74,re75,u74,u75,black,hisp,married,nodegr) %>% 
                   dplyr::mutate(across(age:re75, scale)))

Xt = data.frame(X) %>% dplyr::mutate(treat = lalonde$treat)
glm2 = glm(treat~., data = Xt, family = "binomial")
logLik(glm2)-logLik(glm1) #identical

qr = qr(X)
Q = qr.Q(qr)

eigen(t(Q) %*% Q,symmetric = T, only.values = T) #Yes, returns unitary 

Q = data.frame(Q)
names(Q) = paste0("q",1:ncol(Q))


Qt = Q %>% dplyr::mutate(treat = lalonde$treat)
glm3 = glm(treat~.,data= Qt, family = "binomial")
logLik(glm3)-logLik(glm1) #identical


mplus_dat = lalonde %>% dplyr::select(age,educ,re74,re75,black,hisp,married,nodegr) %>%
  dplyr::mutate(re74 = log(re74+1), re75 = log(re75 + 1)) %>% 
  dplyr::rename(ln74 = re74, ln75 = re75) %>% 
  dplyr::bind_cols(Q) %>% 
  mutate(treat=lalonde$treat, lnyobs = log(lalonde$re78+1), lny0 = ifelse(treat==0,lnyobs,NA), lny1 = ifelse(treat==1,lnyobs,NA)) %>% 
  dplyr::mutate(age = scale(age), educ = scale(educ))


setwd("C:/Users/waldmanm/git-repositories/Causal-Common-Atoms-Mixture/lalonde-pscore/approach1/")

MplusAutomation::prepareMplusData(mplus_dat, filename = "lalondeQPO.dat")



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

median_re75 = with(mplus_dat %>% dplyr::filter(re75!=0 & treat==1), mean(re75))


# Quantiles
mplus_dat = mplus_dat %>%
  dplyr::mutate(re75 = ifelse(re75==0,0,re75), 
                re75 = ifelse(re75>0 & re75<=median_re75, .1, re75), 
                re75 = ifelse(re75>median_re75, .2, re75), 
                re75 = re75*10)

MplusAutomation::prepareMplusData(df = mplus_dat, filename = "data/lalonde_obs.dat")
