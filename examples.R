#Examples

setwd("G:/Mi unidad/Investigación/R/Regresion CTP")

source("ctp.fit.R")
source("summary_ctp.R")
source("plots.R")
source("residuals.R")
source("expected.R")
source("utils.R")

library(GWRM)
library(COUNT)
library(DGLMExtPois)


# 1. Alemania -------------------------------------------------------------

Germany <- read.delim("Germany.txt")
formula<-Scored~Place+Type+Result+Rival
Germany$Place<-relevel(as.factor(Germany$Place),"home")

#CTP estimation
system.time(ctp.ale2 <- ctp.fit(formula, data=Germany))

summary.ctp(ctp.ale2)
#AIC: 3207

#Diagnosis plots
par(mfrow = c(1,2))
plot.ctp(ctp.ale2)

par(mfrow=c(1,1))
expec_ctp2<-expected.ctp(ctp.ale2)

#Goodness of fit measures
MPB(ctp.ale2)
MAD(ctp.ale2)
MSPE(ctp.ale2)
pearson(ctp.ale2)

#number of under-dispersed cases
sum(ctp.ale2$parameters[1]<(-ctp.ale2$fitted.values-1)/2)



#Poisson regression
pois.ale <- glm(formula, data = Germany, family = "poisson")
summary(pois.ale)
#AIC: 3286.4

#Goodness of fit measures
sum(pois.ale$weights*(pois.ale$y-pois.ale$fitted.values))/sum(pois.ale$weights)
# 0.002325247
sum(pois.ale$weights*abs(pois.ale$y-pois.ale$fitted.values))/sum(pois.ale$weights)
# 1.281843
sum(pois.ale$weights*abs(pois.ale$y-pois.ale$fitted.values)^2)/sum(pois.ale$weights)
#[1] 3.156906
sum(residuals(pois.ale)^2)
#[1] 980.8372

library(COUNT)
poisest<-poi.obs.pred(len=17, model=pois.ale)
#Dif
9.93*sum(abs(poisest$propObsv-poisest$propPred))
#[1] 177.1376
#Chi2
9.93*sum((poisest$propObsv-poisest$propPred)^2/poisest$propPred)
#[1] 4548.961



#Binomial Negative regression
library(MASS) #Si cargamos el paquete COUNT no hace falta
nb.ale <- glm.nb(formula, data=Germany)
summary(nb.ale)
#AIC: 3286.9

#Generalized Waring regression ¿Lo dejamos?
gw.ale <- gw(formula, data = Germany)
summary(gw.ale)
#AIC: 3289



#hP regression
library(DGLMExtPois)
hp.ale <- glm.hP(formula.mu = formula, formula.gamma = formula, data=Germany)
summary(hp.ale)
#AIC: 3265


#hP regression with constant dispersion
hp.ale1 <- glm.hP(formula.mu = formula, formula.gamma = ~1, data=Germany)
summary(hp.ale1)
#AIC: 3266

#hP model comparison
DGLMExtPois::lrtest(hp.ale,hp.ale1)
#p-value: 0.0376382


#hP regression with dispersion parameter depending on the significant covariates
system.time(hp.ale2 <- glm.hP(formula.mu = formula, formula.gamma = ~Type+Result, data=Germany))

summary(hp.ale2)
#AIC: 3260

#hP model comparison
DGLMExtPois::lrtest(hp.ale,hp.ale2)
#p-value: 0.6832639

DGLMExtPois::lrtest(hp.ale2,hp.ale1)
#p-value: 0.007856486

#Goodness of fit measures
hp.ale2$response<-hp.ale2$y
MAD(hp.ale2)#[1] 1.10363
MPB(hp.ale2)#[1] -0.005318106
MSPE(hp.ale2)#[1] 2.422297
sum(residuals(hp.ale2)^2)#[1] 1268.918

#number of under-dispersed cases
sum(hp.ale2$gammas < 1)#[1] 854

#Diagnosis plots
par(mfrow = c(1,2))
plot(hp.ale2)

par(mfrow=c(1,1))
expec_hp2<-hP_expected(hp.ale2)
barplot(expec_hp2$observed_freq,xlab = "Y", ylab = "Frequencies")
lines(expec_hp2$frequencies, x=c(1:18), col = "blue")
legend("topright", legend=c(paste("Dif = ", round(expec_hp2$dif, 3)), paste(expression(chi ^ 2), "=", round(37.75822, 3))))
#Dif
sum(abs(expec_hp2$frequencies-expec_hp2$observed_freq))
#[1] 126.124


#Figure 1 in paper
freq_agr<-matrix(c(expec_ctp2$observed_freq,expec_ctp2$frequencies,expec_hp2$frequencies),byrow=T,nrow=3)
colnames(freq_agr)<-0:17
rownames(freq_agr)<-c("Observed","CTP expected","hP expected")
color.names <- c("black","grey50","white")
barplot(freq_agr,beside = T, xlab = "Y", ylab = "Frequencies", col = color.names)
legend("topright",rownames(freq_agr), cex=0.9, fill = color.names, bty = "n")



#CMP regression
cmp.ale <- glm.CMP(formula.mu = formula, formula.nu = formula, data=Germany)
summary(cmp.ale)
#AIC: 12510


#CMP regression with constant dispersion
cmp.ale1 <- glm.CMP(formula.mu = formula, formula.nu = ~1, data=Germany)
summary(cmp.ale1)
#AIC: 3285


#CMP regression with dispersion parameter depending on the significant covariates
cmp.ale2 <- glm.CMP(formula.mu = formula, formula.nu = ~Type+Result, data=Germany)
summary(cmp.ale2)
#AIC: 3283

#cMP model comparison
DGLMExtPois::lrtest(cmp.ale2,cmp.ale1)
#p-value: 0.03799814 

sum(cmp.ale2$nu > 1)#[1] 717

#¿Lo dejamos?
cmp.ale3 <- glm.CMP(formula.mu = formula, formula.nu = ~Oficial, data=Germany)
summary(cmp.ale3)

DGLMExtPois::lrtest(cmp.ale3,cmp.ale2)
#p-value: 0.2788689

DGLMExtPois::lrtest(cmp.ale3,cmp.ale1)
#p-value: 0.01539249
sum(cmp.ale3$nu > 1)#[1] 578



# 2. Data mdvis (COUNT) -------------------------------------------------------------
data(mdvis)

#Poisson regression
md_p <- glm(numvisit ~ reform + factor(educ) + factor(agegrp), family=poisson, data=mdvis)
summary(md_p)
#AIC=13114



#Binomial Negative regression
md_nb <- glm.nb(numvisit ~ reform + factor(educ) + factor(agegrp), data=mdvis)
summary(md_nb)
#AIC=9393.1

#Goodness of fit measures
sum(md_nb$weights*(md_nb$y-md_nb$fitted.values))/sum(md_nb$weights)
#[1] 0.001848847
sum(md_nb$weights*abs(md_nb$y-md_nb$fitted.values))/sum(md_nb$weights)
#[1] 2.389566
sum(md_nb$weights*abs(md_nb$y-md_nb$fitted.values)^2)/sum(md_nb$weights)
#[1] 16.22662
sum(residuals(md_nb)^2)
#[1] 2389.509

library(COUNT)
nbest<-nb2.obs.pred(len=60, model=md_nb)
#Dif
22.27*sum(abs(nbest$propObsv-nbest$propPred))
#[1] 430.0494
#Chi2
22.27*sum((nbest$propObsv-nbest$propPred)^2/nbest$propPred)
#[1] 7861.859



#Generalized Waring regression
system.time(md_gw <- gw(numvisit ~ reform + factor(educ) + factor(agegrp), data=mdvis))

summary(md_gw)
#AIC=9338

#Goodness of fit measures
sum(md_gw$W*(md_gw$Y-md_gw$fitted.values))/sum(md_gw$W)
#[1] 0.02219918
sum(md_gw$W*abs(md_gw$Y-md_gw$fitted.values))/sum(md_gw$W)
#[1] 2.35340
sum(md_gw$W*abs(md_gw$Y-md_gw$fitted.values)^2)/sum(md_gw$W)
#[1] 15.91485
sum(residuals(md_gw)^2)
#[1] 2341.92

#Diagnosis plots
source("plot.glm_gw.R")
source("residuals_gw.R")
source("dgw.R")
par(mfrow = c(1,2))
plot.gw(md_gw)

source("expected_gw.R")
par(mfrow = c(1,1))
expec_gw<-expected.gw(md_gw)
sum(abs(expec_gw$frequencies-expec_gw$observed_freq))
#[1] 352.5629



#CTP regression
#Initial values
ctp<-fitctp(mdvis$numvisit, astart = 1, bstart = 1, gammastart = 3)
ctp$coefficients

system.time(md_ctp <- ctp.fit(numvisit ~ reform + factor(educ) + factor(agegrp), data=mdvis, astart = 2.5, bstart = 0.5))

summary.ctp(md_ctp)
#AIC: 9329


#Diagnosis plots
par(mfrow = c(1,2))
plot.ctp(md_ctp)

par(mfrow=c(1,1))
expec_ctp<-expected.ctp(md_ctp)

#Goodness of fit measures
MPB(md_ctp) #[1] 0.005063565
MAD(md_ctp) #[1] 2.357706
MSPE(md_ctp) #[1] 15.87876
pearson(md_ctp) #[1] 2052.401


#Figure 4 in paper
freq_agr<-matrix(c(expec_ctp$observed_freq[1:22],expec_ctp$frequencies[1:22],expec_gw$frequencies[1:22]),byrow=T,nrow=3)
colnames(freq_agr)<-0:21
rownames(freq_agr)<-c("Observed","CTP expected","GW expected")
color.names <- c("black","grey50","white")
barplot(freq_agr,beside = T, xlab = "Y", ylab = "Frequencies", col = color.names)
legend("topright",rownames(freq_agr), cex=0.9, fill = color.names, bty = "n")



#hP regression
md_hp <- glm.hP(numvisit ~ reform + factor(educ) + factor(agegrp), formula.gamma=~1, data=mdvis)
summary(md_hp)
#AIC: 14



#CMP regression
md_cmp <- glm.CMP(numvisit ~ reform + factor(educ) + factor(agegrp), formula.nu=~1, data=mdvis)
summary(md_cmp)
#AIC: 9415