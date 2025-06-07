# ------------------
# EXAMPLES
# ------------------

setwd(".../Scripts")

source("ctp.fit.R")
source("summary_ctp.R")
source("plots.R")
source("residuals.R")
source("expected.R")
source("utils.R")
source("predict_ctp.R")
source("residuals_gw.R")

library(GWRM)
library(COUNT)
library(DGLMExtPois)


# 1. Germany data
# ------------------------------------------------------------

Germany <- read.delim("Germany.txt")
formula <- Scored ~ Place + Type + Result + Rival
Germany$Place <- relevel(as.factor(Germany$Place),"home")

#CTP regression -----------------------------------------
system.time(ctp.ger2 <- ctp.fit(formula, data = Germany))

summary.ctp(ctp.ger2)
#AIC: 3207

#Diagnostic plots
par(mfrow = c(1,2))
plot.ctp(ctp.ger2)

par(mfrow=c(1,1))
residuals.ctp(ctp.ger2, type = "quantile", envelope = TRUE)

expec_ctp2<-expected.ctp(ctp.ger2)

#Goodness of fit measures
MPB(ctp.ger2) # -0.0214043
MAD(ctp.ger2) # 1.108523
MSPE(ctp.ger2) # 2.43015
pearson(ctp.ger2) # 862.7484

#number of under-dispersed cases
sum(ctp.ger2$parameters[1] < (- ctp.ger2$fitted.values - 1) / 2) # 421

#dispersion in terms of covariates
m <- aggregate(formula, Germany, mean)
v <- aggregate(formula, Germany, var)
l <- aggregate(formula, Germany, length)
estimated.means <- predict.CTP(ctp.ger2, m[,1:5], type ="response")
dispersion.table <- cbind(m[1:4], estimated.means, m[5], v[,5])
names(dispersion.table) <- c("Place", "Type", "Result", "Rival", "Expected mean", "Sample mean", "Sample variance")
dispersion.table <- dispersion.table[order(dispersion.table$`Expected mean`, decreasing = F), ]
dispersion.table

#Cross-validation CTP model
Germany2 <- read.csv("Germany2.txt", sep = "")
y.est <- predict.CTP(ctp.ger2, newdata = Germany2, type = "response")
sum((Germany2$Scored - y.est) ^ 2 ) / 25
#[1] 1.899919

#Poisson regression ----------------------------------------
pois.ger <- glm(formula, data = Germany, family = "poisson")
summary(pois.ger)
#AIC: 3286.4

#Goodness of fit measures
sum(pois.ger$weights * (pois.ger$y - pois.ger$fitted.values)) / sum(pois.ger$weights)
# 0.002325247
sum(pois.ger$weights * abs(pois.ger$y - pois.ger$fitted.values)) / sum(pois.ger$weights)
# 1.281843
sum(pois.ger$weights * abs(pois.ger$y - pois.ger$fitted.values)^2) / sum(pois.ger$weights)
# 3.156906
sum(residuals(pois.ger)^2)
# 980.8372

poisest <- poi.obs.pred(len = 17, model = pois.ger)
#Dif
9.93*sum(abs(poisest$propObsv - poisest$propPred))
# 177.1376
#Chi2
9.93*sum((poisest$propObsv - poisest$propPred)^2 / poisest$propPred)
# 4548.961


#Negative binomial regression -----------------------------------------
library(MASS) #If we load the COUNT package, we do not need to
nb.ger <- glm.nb(formula, data = Germany)
summary(nb.ger)
#AIC: 3286.9

#Generalized Waring regression ----------------------------------------
gw.ger <- gw(formula, data = Germany)
summary(gw.ger)
#AIC: 3289



#hP regression ----------------------------------------------------------------
hp.ger <- glm.hP(formula.mu = formula, formula.gamma = formula, data = Germany)
summary(hp.ger)
#AIC: 3265

#hP regression with constant dispersion ------------------------------------
hp.ger1 <- glm.hP(formula.mu = formula, formula.gamma = ~ 1, data = Germany)
summary(hp.ger1)
#AIC: 3266

#hP model comparison
DGLMExtPois::lrtest(hp.ger,hp.ger1)
#p-value: 0.0376382

#hP regression with dispersion parameter depending on the significant covariates
system.time(hp.ger2 <- glm.hP(formula.mu = formula, formula.gamma = ~ Type + Result, data = Germany))

summary(hp.ger2)
#AIC: 3260

#hP model comparison
DGLMExtPois::lrtest(hp.ger, hp.ger2)
#p-value: 0.6832639

DGLMExtPois::lrtest(hp.ger2, hp.ger1)
#p-value: 0.007856486

#Goodness of fit measures
hp.ger2$response <- hp.ger2$y
MAD(hp.ger2) # 1.10363
MPB(hp.ger2) # -0.005318106
MSPE(hp.ger2) # 2.422297
sum(residuals(hp.ger2)^2) # 1268.918

#number of under-dispersed cases
sum(hp.ger2$gammas < 1) # 854

#Diagnostic plots
par(mfrow = c(1,2))
plot(hp.ger2)

par(mfrow = c(1,1))
residuals(hp.ger2, type = "quantile", envelope = TRUE)

expec_hp2 <- hP_expected(hp.ger2)

barplot(expec_hp2$observed_freq, xlab = "Y", ylab = "Frequencies")
lines(expec_hp2$frequencies, x = c(1:18), col = "blue")
legend("topright", legend = c(paste("Dif = ", round(expec_hp2$dif, 3)), paste(expression(chi ^ 2), "=", round(37.75822, 3))))
#Dif
sum(abs(expec_hp2$frequencies - expec_hp2$observed_freq))
# 126.124


#Figure 1 in paper
freq_agr <- matrix(c(expec_ctp2$observed_freq, expec_ctp2$frequencies, expec_hp2$frequencies), byrow=T, nrow=3)
colnames(freq_agr) <- 0:17
rownames(freq_agr) <- c("Observed","CTP expected","hP expected")
color.names <- c("black","grey50","white")
barplot(freq_agr,beside = T, xlab = "Y", ylab = "Frequencies", col = color.names)
legend("topright",rownames(freq_agr), cex = 0.9, fill = color.names, bty = "n")



#CMP regression --------------------------------------------------------------
cmp.ger <- glm.CMP(formula.mu = formula, formula.nu = formula, data = Germany)
summary(cmp.ger)
#AIC: 12510


#CMP regression with constant dispersion ----------------------------------
cmp.ger1 <- glm.CMP(formula.mu = formula, formula.nu = ~ 1, data = Germany)
summary(cmp.ger1)
#AIC: 3285


#CMP regression with dispersion parameter depending on the significant covariates -----
cmp.ger2 <- glm.CMP(formula.mu = formula, formula.nu = ~ Type + Result, data = Germany)
summary(cmp.ger2)
#AIC: 3283

#cMP model comparison
DGLMExtPois::lrtest(cmp.ger2, cmp.ger1)
#p-value: 0.03799814 

sum(cmp.ger2$nu > 1) # 717



# 2. mdvis data (COUNT) 
# -------------------------------------------------------------

data(mdvis)

#Poisson regression --------------------------------------------------------------------------
md_p <- glm(numvisit ~ reform + factor(educ) + factor(agegrp), family = poisson, data = mdvis)
summary(md_p)
#AIC: 13114



#Negative binomial regression --------------------------------------------------
md_nb <- glm.nb(numvisit ~ reform + factor(educ) + factor(agegrp), data = mdvis)
summary(md_nb)
#AIC: 9393.1

#Goodness of fit measures
sum(md_nb$weights*(md_nb$y - md_nb$fitted.values)) / sum(md_nb$weights)
# 0.001848847
sum(md_nb$weights*abs(md_nb$y-md_nb$fitted.values))/sum(md_nb$weights)
# 2.389566
sum(md_nb$weights*abs(md_nb$y-md_nb$fitted.values)^2)/sum(md_nb$weights)
# 16.22662
sum(residuals(md_nb)^2)
# 2389.509

nbest<-nb2.obs.pred(len=60, model=md_nb)
#Dif
22.27 * sum(abs(nbest$propObsv-nbest$propPred))
# 430.0494
#Chi2
22.27 * sum((nbest$propObsv-nbest$propPred)^2/nbest$propPred)
# 7861.859



#Generalized Waring regression ----------------------------------------------------------
system.time(md_gw <- gw(numvisit ~ reform + factor(educ) + factor(agegrp), data = mdvis))

summary(md_gw)
#AIC: 9338

#Goodness of fit measures
sum(md_gw$W * (md_gw$Y - md_gw$fitted.values)) / sum(md_gw$W)
# 0.02219918
sum(md_gw$W*abs(md_gw$Y-md_gw$fitted.values))/sum(md_gw$W)
# 2.35340
sum(md_gw$W*abs(md_gw$Y-md_gw$fitted.values)^2)/sum(md_gw$W)
# 15.91485
sum(residuals(md_gw)^2)
# 2341.92

#Diagnostic plots
source("plot.glm_gw.R")
source("residuals_gw.R")
source("dgw.R")
source("pgw.R")
par(mfrow = c(1,2))
plot.gw(md_gw)


par(mfrow = c(1,1))
residuals.gw(md_gw, type = "quantile", envelope = TRUE)

source("expected_gw.R")
par(mfrow = c(1,1))
expec_gw <- expected.gw(md_gw)
sum(abs(expec_gw$frequencies - expec_gw$observed_freq))
# 352.5629



#CTP regression -----------------------------------------------------
#Initial values
ctp <- fitctp(mdvis$numvisit, astart = 1, bstart = 1, gammastart = 3)
ctp$coefficients

system.time(md_ctp <- ctp.fit(numvisit ~ reform + factor(educ) + factor(agegrp), data = mdvis, astart = 2.5, bstart = 0.5))

summary.ctp(md_ctp)
#AIC: 9329

#Diagnostic plots
par(mfrow = c(1,2))
plot.ctp(md_ctp)

par(mfrow = c(1,1))
residuals.ctp(md_ctp, type = "quantile", envelope = TRUE)

expec_ctp <- expected.ctp(md_ctp)

#Goodness of fit measures
MPB(md_ctp) # 0.005063565
MAD(md_ctp) # 2.357706
MSPE(md_ctp) # 15.87876
pearson(md_ctp) # 2052.401


#Figure 6 in paper
freq_agr <- matrix(c(expec_ctp$observed_freq[1:22], expec_ctp$frequencies[1:22], expec_gw$frequencies[1:22]), byrow = T, nrow = 3)
colnames(freq_agr) <- 0:21
rownames(freq_agr) <- c("Observed", "CTP expected", "GW expected")
color.names <- c("black","grey50","white")
barplot(freq_agr,beside = T, xlab = "Y", ylab = "Frequencies", col = color.names)
legend("topright",rownames(freq_agr), cex = 0.9, fill = color.names, bty = "n")

#Leave-one-out cross-validation for CTP model
library(doParallel)
registerDoParallel(cores=8)

s <- nrow(mdvis)
y_jack <- c()

LOOCV.mdvis <- function(k){
  datak <- mdvis[-k,]
  library(hypergeo)
  library(gsl)
  library(nloptr)
  library(pracma)
  library(cpd)
  md_ctp.fit <- ctp.fit(numvisit ~ reform + factor(educ) + factor(agegrp), data = datak, astart = 2.5, bstart = 0.5)
  y_jack <- predict.CTP(md_ctp.fit, newdata = mdvis, type = "response")[k, 1]
  list(ajuste = md_ctp.fit, y.jack = y_jack)
}

md_ctp.LOOCV <- foreach(i=1:s) %dopar% {
  LOOCV.mdvis(i) 
}

md_ctp.LOOCV

yjack <- c()
for (i in 1:2227){
  yjack[i] <- md_ctp.LOOCV[[i]]$y.jack
}

sum((mdvis$numvisit-yjack)^2)/s
#[1] 15.96577


#hP regression --------------------------------------------------------------------------------------
md_hp <- glm.hP(numvisit ~ reform + factor(educ) + factor(agegrp), formula.gamma = ~ 1, data = mdvis)
summary(md_hp)
#AIC: 14



#CMP regression ------------------------------------------------------------------------------------
md_cmp <- glm.CMP(numvisit ~ reform + factor(educ) + factor(agegrp), formula.nu = ~ 1, data = mdvis)
summary(md_cmp)
#AIC: 9415