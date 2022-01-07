# setwd("C:/Users/johan/Google Drive/UNI/UPC/Lifetime Data Analysis/Lifetime-Data-Analysis")

# Objective: To analyze the survival time in patients with HIV and
# low CD4 cell count. Study the association between
# the survival time and the five variables contained in the
# data set.

setwd("~/Downloads/Lab class 8 21.12.2021-20220104")
###### IMPORT ######
library(survival)
library(KMsurv)
library(FHtest)

## import data ##
names=list("id", "survival", "cens", "CD4_cells", "treatment","gender", "AIDS", "AZT")
aids=read.table("AIDS.txt", skip=28, col.names=names)
aids=subset(aids,select=-id)
aids$gender=as.factor(aids$gender)
aids$AIDS=as.factor(aids$AIDS)
aids$AZT=as.factor(aids$AZT)
aids$treatment=as.factor(aids$treatment)
attach(aids)


##### 2.DESCRIPTIVE ANALYSIS #####

### univariate analysis

summary(aids)

## survival
hist(survival, main="Distribution of Survival Times", xlab="Survival Times")
# distribution
qqplot(survival, rnorm(length(survival)))
shapiro.test(survival)
#clearly not normal
boxplot(survival, main="Boxplot of Survivaltimes")

#cens
plot(as.factor(cens), main="Censoring of the data")
#proportion of censoring
1-sum(aids$cens)/length(cens)

#CD4_cells
b=boxplot(CD4_cells, main="Boxplot of CD4 cells")
hist(CD4_cells)
lamb=mean(CD4_cells)
#distribution
qqplot(CD4_cells,rpois(length(CD4_cells),1/lamb))
qqplot(CD4_cells,rexp(length(CD4_cells), lamb))

#treatment
plot(treatment, main="Treatment Types")
#gender
plot(gender, main="Distribution of Gender")
#AIDS
plot(AIDS,main="Distribution of AIDS" )
#AZT
plot(AZT,main="Distribution of AZT" )

### bivariate analysis
# install.packages("corrplot")
# install.packages("PerformanceAnalytics")
library(corrplot)
library(PerformanceAnalytics)

aids.corr=aids
aids.corr$treatment=as.integer(aids.corr$treatment)-1
aids.corr$gender=as.integer(aids.corr$gender)-1
aids.corr$AZT=as.integer(aids.corr$AZT)-1
aids.corr$AIDS=as.integer(AIDS,levels=c("noAIDS","AIDS"))%%2
chart.Correlation(aids.corr, histogram=TRUE, pch=19)
corrplot(cor(aids.corr))

plot(survival,CD4_cells, main="Suvival times vs. CD4 cells")
boxplot(survival~gender, main="Boxplot of Survival time by Gender")
boxplot(survival~AIDS, main="Boxplot of Survival time by AIDS")
boxplot(survival~treatment, main="Boxplot of Survival time by Treatment")
boxplot(survival~AZT, main="Boxplot of Survival time by AZT")



##### 3.NON-PARAMETRIC ANALYIS #####

# The Surv object with right-censored data
with(aids, Surv(survival, cens))

### 3.1 overall survival function
svf <- survfit(Surv(survival, cens) ~ 1, aids) #KM estimator
svfNA= survfit(Surv(survival, cens) ~ 1, aids, type = "fleming") #NA estimator

# both the same

par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l")
plot(svf, col = 3, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25), conf.int=FALSE, main="Survival Functin AIDS patients \n (KM estimator)")
lines(svfNA, col = 2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25), conf.int=FALSE, main="Survival Functin AIDS patients \n (NA estimator)")
axis(1, seq(0, 25, 5))

#survival function wtr treatment groups, CD4_cells, gender, AIDS, AZT
svf_tr <- survfit(Surv(survival, cens) ~ treatment, aids)
svf_ge <- survfit(Surv(survival, cens) ~ gender, aids)
svf_AI <- survfit(Surv(survival, cens) ~ AIDS, aids)
svf_AZ <- survfit(Surv(survival, cens) ~ AZT, aids)

#plot
par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l", mfrow=c(2,2))
plot(svf_tr, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Functin AIDS patients")
legend("bottomleft", levels(aids$treatment), lwd = 3, col = 1:2, bty = "n")

plot(svf_ge, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Functin AIDS patients")
legend("bottomleft", levels(aids$gender), lwd = 3, col = 1:2, bty = "n")

plot(svf_AI, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Functin AIDS patients")
legend("bottomleft", levels(aids$AIDS), lwd = 3, col = 1:2, bty = "n")

plot(svf_AZ, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Functin AIDS patients")
legend("bottomleft", levels(aids$AZT), lwd = 3, col = 1:2, bty = "n")

### 3.2 estimatin of median survival time
svf
quantile(svf)
#median 19.1

#distribution function
plot(svf, col = 3, xlab = "Time to Death [Years]",
     ylab = expression(bolditalic(hat(F)(t))), xaxs = "i",xlim=c(0,25),ylim=c(0,0.6),
     lty = 1, conf.int=FALSE, fun = "event", lwd = 3, yaxs = "i", bty = "n", main="Cummulative Distribution of time to death")
abline(v=19.1, col=2, lwd=2, lty=1)
abline(h=0.5, col=2, lwd=1)

### 3.3 Comparison of survival functions by means of nonparametric tests such as the logrank test


#logrank test
svflr <- with(aids, Surv(survival, cens) ~ treatment)
lrt=survdiff(svflr)
lrt

#other non-parametric tests
survdiff(svflr, rho = 1)
survdiff(svflr, rho = -1)

#Fleming-Harrington class
FHtestrcc(svflr)
FHtestrcc(svflr, rho = 1)
FHtestrcc(svflr, lambda = 1)
FHtestrcc(svflr, rho = 1, lambda = 1)
FHtestrcc(svflr, rho = 0.5, lambda = 2)

#stratified test
#TODO: is balanced study?
staids <- with(aids, Surv(survival, cens) ~ treatment + AZT)
survfit(staids)
plot(staids)

##### 4.PARAMETRIC SURVIVAL MODEL #####

### 4.1: Fit of a Weibull, log-logistic, or lognormal model.

#weibull models
weimod1 <- survreg(Surv(survival, cens) ~ 1, dist = "weibull")
weimod2 <- survreg(Surv(survival, cens) ~ treatment, dist = "weibull")
weimod3 <- survreg(Surv(survival, cens) ~ AZT, dist = "weibull")
weimod4 <- survreg(Surv(survival, cens) ~ gender, dist = "weibull")
weimod5 <- survreg(Surv(survival, cens) ~ AIDS, dist = "weibull")
weimod6 <- survreg(Surv(survival, cens) ~ CD4_cells, dist = "weibull")

summary(weimod1)
summary(weimod2)
summary(weimod3)
summary(weimod4)
summary(weimod5)
summary(weimod6)

anova(weimod3)

#weibull: check fit
wm3pred <- predict(weimod3, type = "linear")
wm3pred
resids <- (log(aids$survival) - wm3pred) / weimod3$scale
resids

par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(resids, aids$cens) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the Weibull regression model")

# Graphical comparison with the theoretical survival function
survgumb <- function(x) {
  return(exp(-exp(x)))
}

curve(survgumb(x), from = min(resids), to = max(resids), col = 2, lwd = 3,
      add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Gumbel Distribution"),
       col = c(1, 1, 2), lty = c(1, 2, 1), lwd = 3, bty = "n")


#loglogistic models
loglomod1 <- survreg(Surv(survival, cens) ~ 1, dist = "loglogistic")
loglomod2 <- survreg(Surv(survival, cens) ~ treatment, dist = "loglogistic")
loglomod3 <- survreg(Surv(survival, cens) ~ AZT, dist = "loglogistic")
loglomod4 <- survreg(Surv(survival, cens) ~ gender, dist = "loglogistic")
loglomod5 <- survreg(Surv(survival, cens) ~ AIDS, dist = "loglogistic")
loglomod6 <- survreg(Surv(survival, cens) ~ CD4_cells, dist = "loglogistic")

summary(loglomod1)
summary(loglomod2)
summary(loglomod3)
summary(loglomod4)
summary(loglomod5)
summary(loglomod6)

#loglogistic: model fit

lnopred <- predict(loglomod1, type = "linear")
residsLN <- (log(aids$survival) - lnopred) / loglomod1$scale
residsLN

#plot residual
par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(residsLN, aids$cens) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the lognormal regression model")

# Adding the theoretical survival function

curve(pnorm(x, lower.tail = FALSE), from = min(residsLN), to = max(residsLN),
      col = 2, lwd = 3, add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Normal Distribution"),
       col = c(1, 1, 2), lty = c(1, 2, 1), lwd = 3, bty = "n")


#lognormal models
lognomod1 <- survreg(Surv(survival, cens) ~ 1, dist = "lognormal")
lognomod2 <- survreg(Surv(survival, cens) ~ treatment, dist = "lognormal")
lognomod3 <- survreg(Surv(survival, cens) ~ AZT, dist = "lognormal")
lognomod4 <- survreg(Surv(survival, cens) ~ gender, dist = "lognormal")
lognomod5 <- survreg(Surv(survival, cens) ~ AIDS, dist = "lognormal")
lognomod6 <- survreg(Surv(survival, cens) ~ CD4_cells, dist = "lognormal")

summary(lognomod1)
summary(lognomod2)
summary(lognomod3)
summary(lognomod4)
summary(lognomod5)
summary(lognomod6)


#lognormal: check fit

lnopred <- predict(lognomod3, type = "linear")
residsLN <- (log(aids$survival) - lnopred) / lognomod3$scale
residsLN

#plot residual
par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(residsLN, aids$cens) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the lognormal regression model")

# Adding the theoretical survival function
curve(pnorm(x, lower.tail = FALSE), from = min(residsLN), to = max(residsLN),
      col = 2, lwd = 3, add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Normal Distribution"),
       col = c(1, 1, 2), lty = c(1, 2, 1), lwd = 3, bty = "n")

##### 5. SEMI-PARAMETRIC MODEL #####
##### 5.1 Fit of the proportional hazards models.
## The fit of a Cox model
svfc <- with(myeloma, Surv(survival, cens))
coxm <- coxph(svfc ~ CD4_cells + treatment + gender + AIDS + AZT, aids)
summary(coxm)

coxm$var

#####  5.2 Interpretation of the model fit.

# In here we try to interpret the model.

#####  5.3 Interpretation of the model parameters in terms of relative hazards.
# which is the ratio of the hazard at time t to the baseline hazard, to the risk factors
library(Epi)
ci.lin(coxm)
round(ci.lin(coxm, Exp = TRUE), 3)
(ctmat <- matrix(c(1, 0, 0, 0, 0,
                   0, 1, 0, 0, 0,
                   0, 0, 1, 0, 0,
                   0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), byrow = TRUE, nr = 5))
round(ci.lin(coxm, ctr.mat = ctmat, Exp = TRUE), 3)

# A somewhat nicer presentation
HRmat1 <- round(ci.lin(coxm, ctr.mat = ctmat, Exp = TRUE), 3)[, c(1, 5:7)]
rownames(HRmat1) <- c("HR| CD4_cells", "HR| treatment ddI", "HR| gender male", "HR| noAIDS", "HR| AZT intolerance")
colnames(HRmat1) <- c("logHR", "HR", "Lower 95%", "Upper 95%")
HRmat1

par(font = 2, font.lab = 4, font.axis = 2, las = 1)
plot(survfit(coxm), lwd = 2, col = 3)

#####  5.4 Analysis of the residuals to check the model’s goodness-of-fit.
# What are the residuals for the Cox model?
# (a) generalized (Cox-Snell)
# (b) Schoenfeld
# (c) martingale

residuals(coxm)

resids1 <- residuals(update(coxm,  ~ . - CD4_cells))


par(mfrow = c(2, 3), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
plot(resids1 ~ aids$CD4_cells, xlab = "CD4_cells", ylab = "Residuals", pch = 19)
abline(h = 0, lwd = 2, lty = 2)
lines(lowess(aids$CD4_cells, resids1), lwd = 3)

## The dfbeta residuals: looking for influential observations.

dfbet <- residuals(coxm, type = "dfbeta")
dim(dfbet)
par(mfrow = c(2,3), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
for (i in 1:5) {
  plot(dfbet[, i], pch = 16, ylab = "")
  title(names(coef(coxm))[i])
  axis(1, at = seq(5, 45, 5))
  abline(h=0)
}

## Closer look at the influential observations on the estimation of a
## couple of parameters.
## Variable CD4_cells:
par(font = 2, font.lab = 4, font.axis = 2, las = 1, cex.lab = 1.3,
    cex.axis = 1.2)
plot(dfbet[, 3], pch = 16, ylab = "")
title(names(coef(coxm))[3])
axis(1, at = seq(5, 45, 5))
identify(dfbet[, 3],
         labels = paste0("Row: ", rownames(aids), "; Time: ", aids$survival,
                         "; Cens: ", aids$cens, "; CD4_cells: ", aids$CD4_cells))

# ID 38 is an influential observation:
aids[263, ]     
# Estimation of beta_3 excluding ID263
coxm1 <- coxph(Surv(survival, cens) ~ CD4_cells + treatment + gender + AIDS + AZT, data = aids[-263, ])
summary(coxm1)
## Variable gender:
par(font = 2, font.lab = 4, font.axis = 2, las = 1, cex.lab = 1.3,
    cex.axis = 1.2, mfrow = c(1,1))
plot(dfbet[, 4], pch = 16, ylab = "")
title(names(coef(coxm))[4])
axis(1, at = seq(5, 45, 5))
identify(dfbet[, 4],
         labels = paste0("Row: ", rownames(aids), "; Time: ", aids$survival,
                         "; Cens: ", aids$cens, "; gender: ", aids$gender))

# ID 13 is an influential observation
aids[c(14, 49, 261), ]

# Estimation of beta_4 excluding ID13
coxm2 <- coxph(Surv(survival, cens) ~ CD4_cells + treatment + gender + AIDS + AZT, data = aids[-c(14, 49, 261), ])
summary(coxm2)

## The Schoenfeld residuals:
## Can we assume the assumption of proportional hazards hold?

residuals(coxm, "schoenfeld")
## (i) Use of function cox.zph
cox.zph(coxm)

## (ii) Use of function plot.cox.zph
par(mfrow = c(2, 3), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
plot(cox.zph(coxm), lwd = 2)

library(tidyverse)

aids$resid_mart <- residuals(coxm, type = "martingale")

ggplot(data = aids, mapping = aes(x = CD4_cells, y = resid_mart)) +
  geom_point() +
  geom_smooth() +
  labs(title = "CD4_cells") +
  theme_bw() + theme(legend.key = element_blank())
