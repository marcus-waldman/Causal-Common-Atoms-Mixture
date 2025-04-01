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
Q = Q %>% dplyr::mutate(across(everything(),scale))


Qt = Q %>% dplyr::mutate(treat = lalonde$treat)
glm3 = glm(treat~.,data= Qt, family = "binomial")
logLik(glm3)-logLik(glm1) #identical

median_re75 = with(lalonde %>% dplyr::filter(re75!=0 & treat==1), median(re75))
median_re74 = with(lalonde %>% dplyr::filter(re74!=0 & treat==1), median(re74))
median_re78 = with(lalonde %>% dplyr::filter(re78!=0 & treat==1), median(re78))




mplus_dat = lalonde %>% 
  dplyr::mutate(ue74 = ifelse(re74==0,0, ifelse(re74>0 & re74<=median_re74,1,2)) ) %>% 
  dplyr::mutate(ue75 = ifelse(re75==0,0, ifelse(re75>0 & re75<=median_re75,1,2)) ) %>% 
  dplyr::mutate(
      ue780 = NA,
      ue780 = ifelse(treat==0 & re78==0, 0, ue780), 
      ue780 = ifelse(treat==0 & re78>0, 1, ue780), 
      ue780 = ifelse(treat==0 & re78>median_re78, 2, ue780)
  ) %>% 
  dplyr::mutate(
    ue781 = NA,
    ue781 = ifelse(treat==1 & re78==0, 0, ue781), 
    ue781 = ifelse(treat==1 & re78>0, 1, ue781), 
    ue781 = ifelse(treat==1 & re78>median_re78, 2, ue781)
  ) %>% 
  dplyr::select(treat, ue780, ue781,age,educ,ue74,ue75,black,hisp,married) %>%
  dplyr::bind_cols(Q) %>% 
  dplyr::mutate(age = scale(age), educ = scale(educ))


setwd("C:/Users/waldmanm/git-repositories/Causal-Common-Atoms-Mixture/lalonde-pscore/approach2/")

#MplusAutomation::prepareMplusData(mplus_dat, filename = "lalondeQPO.dat")


# Read in 5 class model

mod = MplusAutomation::readModels(target = getwd(), filefilter = "approach2-5")

tmp1 = mod$savedata %>% dplyr::mutate(rid = 1:n()) %>% 
  dplyr::filter(A==1) %>% 
  dplyr::select(rid,CPROB1,CPROB3,CPROB5,CPROB7,CPROB9) %>% 
  dplyr::rename(pi1 = CPROB1, pi2 = CPROB3, pi3 = CPROB5, pi4 = CPROB7, pi5 = CPROB9)

tmp2 = mod$savedata %>% dplyr::mutate(rid = 1:n()) %>% 
  dplyr::filter(A==2) %>% 
  dplyr::select(rid,CPROB2,CPROB4,CPROB6,CPROB8,CPROB10) %>% 
  dplyr::rename(pi1 = CPROB2, pi2 = CPROB4, pi3 = CPROB6, pi4 = CPROB8, pi5 = CPROB10)

pis = tmp1 %>% dplyr::bind_rows(tmp2) %>% dplyr::arrange(rid) %>% 
  dplyr::mutate(sum_pi = pi1+pi2+pi3+pi4+pi5) %>% 
  dplyr::mutate(pi1 = pi1/sum_pi, pi2 = pi2/sum_pi, pi3 = pi3/sum_pi, pi4 = pi4/sum_pi, pi5 = pi5/sum_pi) %>% 
  dplyr::mutate(sum_pi = pi1+pi2+pi3+pi4+pi5) 
  

lalonde = lalonde %>% dplyr::bind_cols(pis %>% dplyr::select(pi1:pi5))

gamma_hat = mod$parameters$unstandardized %>% 
  dplyr::filter(paramHeader == "A#1.ON") %>% 
  dplyr::bind_rows(
    mod$parameters$unstandardized %>% 
      dplyr::filter(param == "A#1")
  )


X1 = data.frame(c1=1,c2=0,c3 = 0, c4=0) %>% dplyr::bind_cols(Q) %>% dplyr::mutate(intercpt = 1) %>% as.matrix()
X2 = data.frame(c1=0,c2=1,c3 = 0, c4=0) %>% dplyr::bind_cols(Q)%>% dplyr::mutate(intercpt = 1) %>% as.matrix()
X3 = data.frame(c1=0,c2=0,c3 = 1, c4=0) %>% dplyr::bind_cols(Q) %>% dplyr::mutate(intercpt = 1) %>% as.matrix()
X4 = data.frame(c1=0,c2=0,c3 = 0, c4=1) %>% dplyr::bind_cols(Q) %>% dplyr::mutate(intercpt = 1) %>% as.matrix()
X5 = data.frame(c1=0,c2=0,c3 = 0, c4=0) %>% dplyr::bind_cols(Q) %>% dplyr::mutate(intercpt = 1) %>% as.matrix()

lalonde = lalonde %>% 
  dplyr::mutate(
    e1 = X1 %*% gamma_hat$est,
    e2 = X2 %*% gamma_hat$est,
    e3 = X3 %*% gamma_hat$est,
    e4 = X4 %*% gamma_hat$est,
    e5 = X5 %*% gamma_hat$est
  )

lalonde = lalonde %>% 
  dplyr::mutate(e_new =plogis( pi1*e1 + pi2*e2 + pi3*e3 + pi4*e4 + pi5*e5) )


ggplot(lalonde, aes(x= e, y = e_new)) + geom_point()

ggplot(lalonde %>% dplyr::mutate(Arm = ifelse(treat==1,"A=1","A=0")) %>% dplyr::select(Arm,e,e_new) %>% tidyr::pivot_longer(e:e_new, values_to = "e")) +
  geom_histogram(aes(e)) +
  facet_grid(Arm~name)
  
lalonde = lalonde %>% dplyr::mutate(
  iptw_new = (n_a/N)/(treat+1-e_new), iptw_new = iptw_new/mean(iptw_new)
)


ggplot(lalonde %>% dplyr::mutate(Arm = ifelse(treat==1,"A=1","A=0")) %>% dplyr::select(Arm,iptw,iptw_new) %>% tidyr::pivot_longer(iptw:iptw_new, values_to = "w")) +
  geom_histogram(aes(w)) +
  facet_grid(name~Arm, scales = "free_x") 


ggplot(lalonde %>% dplyr::select(e,e_new) %>% tidyr::pivot_longer(names_to = ""),)