setwd("C:/Users/johan/Google Drive/UNI/UPC/Lifetime Data Analysis/Lifetime-Data-Analysis")

# Objective: To analyze the survival time in patients with HIV and
# low CD4 cell count. Study the association between
# the survival time and the five variables contained in the
# data set.

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
#install.packages("corrplot")
#install.packages("PerformanceAnalytics")
library(corrplot)
library(PerformanceAnalytics)

aids.corr=aids
aids.corr$treatment=as.integer(aids.corr$treatment)-1
aids.corr$gender=as.integer(aids.corr$gender)-1
aids.corr$AZT=as.integer(aids.corr$AZT)-1
aids.corr$AIDS=(as.integer(aids.corr$AIDS))%%2
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
quantile(svf)

# both the same

par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l")
plot(svf, col = 3, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25), conf.int=TRUE, main="Survival Function AIDS patients \n (KM estimator)")
lines(svfNA, col = 2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25), conf.int=FALSE, main="Survival Function AIDS patients \n (NA estimator)")
axis(1, seq(0, 25, 5))

#survival function wtr treatment groups, CD4_cells, gender, AIDS, AZT

#bin the variable CD4 count according to the quantiles
CD4_cat <- cut(CD4_cells, breaks=c(-Inf, 11, 37,109,Inf), labels=c("0-11", "12-37","38-109","110-307"))

svf_tr <- survfit(Surv(survival, cens) ~ treatment, aids)
svf_ge <- survfit(Surv(survival, cens) ~ gender, aids)
svf_AI <- survfit(Surv(survival, cens) ~ AIDS, aids)
svf_AZ <- survfit(Surv(survival, cens) ~ AZT, aids)
svf_CD <- survfit(Surv(survival, cens) ~ CD4_cat, aids)

svf_tr
svf_ge
svf_AI
svf_AZ
svf_CD

#plot
par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l",mfrow=c(1,2))
plot(svf_tr, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Function AIDS patients")
legend("bottomleft", levels(aids$treatment), lwd = 3, col = 1:2, bty = "n")

plot(svf_ge, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Function AIDS patients")
legend("bottomleft", levels(aids$gender), lwd = 3, col = 1:2, bty = "n")

plot(svf_AI, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Function AIDS patients")
legend("bottomleft", levels(aids$AIDS), lwd = 3, col = 1:2, bty = "n")

plot(svf_AZ, col = 1:2, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Function AIDS patients")
legend("bottomleft", levels(aids$AZT), lwd = 3, col = 1:2, bty = "n")

plot(svf_CD, col = 1:4, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25),main="Survival Function CD4 cell counts")
legend("bottomleft", levels(CD4_cat), lwd = 3, col = 1:4, bty = "n")


### 3.2 estimatin of median survival time
svf
quantile(svf)
#median 19.1
quantile(svf_tr)
quantile(svf_ge)
summary(svf_AZ[2])$surv
quantile(svf_AI)
quantile(svf_CD)

#distribution function
plot(svf, col = 3, xlab = "Time to Death [Years]",
     ylab = expression(bolditalic(hat(F)(t))), xaxs = "i",xlim=c(0,25),ylim=c(0,0.6),
     lty = 1, conf.int=FALSE, fun = "event", lwd = 3, yaxs = "i", bty = "n", main="Cummulative Distribution of time to death")
abline(v=19.1, col=2, lwd=2, lty=1)
abline(h=0.5, col=2, lwd=1)

### 3.3 Comparison of survival functions by means of nonparametric tests such as the logrank test


#logrank test
svflr_tr <- with(aids, Surv(survival, cens) ~ treatment)
svflr_ge <- with(aids, Surv(survival, cens) ~ gender)
svflr_AI <- with(aids, Surv(survival, cens) ~ AIDS)
svflr_CD <- with(aids, Surv(survival, cens) ~ CD4_cat)
svflr_AZ <- with(aids, Surv(survival, cens) ~ AZT)
lrt_tr=survdiff(svflr_tr)
lrt_ge=survdiff(svflr_ge)
lrt_AI=survdiff(svflr_AI)
lrt_CD=survdiff(svflr_CD)
lrt_AZ=survdiff(svflr_AZ)
lrt_tr #not reject
lrt_ge #not reject
lrt_AI #reject
lrt_CD #reject
lrt_AZ #reject
#not reject

#other non-parametric tests
survdiff(svflr, rho = 1)
survdiff(svflr, rho = -1)

#Fleming-Harrington class
FHtestrcc(svflr)
FHtestrcc(svflr, rho = 1)
FHtestrcc(svflr, lambda = 1)
FHtestrcc(svflr, rho = 1, lambda = 1)
FHtestrcc(svflr, rho = 0.5, lambda = 2)

##### 4.PARAMETRIC SURVIVAL MODEL #####

### 4.1: Fit of a Weibull, log-logistic, or lognormal model.

#weibull models
weimod1 <- survreg(Surv(survival, cens) ~ treatment+AIDS+CD4_cells+AZT+gender, dist = "weibull")
weimod2 <- survreg(Surv(survival, cens) ~ AIDS, dist = "weibull")
weimod3 <- survreg(Surv(survival, cens) ~ CD4_cat, dist = "weibull")
weimod4 <- survreg(Surv(survival, cens) ~ AZT, dist = "weibull")

summary(weimod1)
summary(weimod2)
summary(weimod3)
summary(weimod1)

anova(weimod1)

#weibull: check fit
wm1pred <- predict(weimod1, type = "linear")
wm1pred
resids <- (log(aids$survival) - wm1pred) / weimod1$scale
resids

par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),mfrow=c(1,3),
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
loglomod1 <- survreg(Surv(survival, cens) ~ AZT+AIDS+gender+treatment+CD4_cells, dist = "loglogistic")
loglomod2 <- survreg(Surv(survival, cens) ~ AIDS, dist = "loglogistic")
loglomod3 <- survreg(Surv(survival, cens) ~ CD4_cells, dist = "loglogistic")
loglomod4 <- survreg(Surv(survival, cens) ~ AZT, dist = "loglogistic")


summary(loglomod1)
summary(loglomod2)
summary(loglomod3)
summary(loglomod4)

anova(loglomod1)
#loglogistic: model fit

lnopred <- predict(loglomod1, type = "linear")
residsLN <- (log(aids$survival) - lnopred) / loglomod1$scale
residsLN

#plot residual
par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(residsLN, aids$cens) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the loglogistic regression model")

# Adding the theoretical survival function

curve(pnorm(x, lower.tail = FALSE), from = min(residsLN), to = max(residsLN),
      col = 2, lwd = 3, add = TRUE)
legend("bottomleft", c("KM estimate", "95% - CI", "Stand. Normal Distribution"),
       col = c(1, 1, 2), lty = c(1, 2, 1), lwd = 3, bty = "n")


#lognormal models
lognomod1 <- survreg(Surv(survival, cens) ~ AZT+AIDS+treatment+gender+CD4_cells, dist = "lognormal")
lognomod2 <- survreg(Surv(survival, cens) ~ AIDS, dist = "lognormal")
lognomod3 <- survreg(Surv(survival, cens) ~ CD4_cells, dist = "lognormal")
lognomod4 <- survreg(Surv(survival, cens) ~ AZT, dist = "lognormal")

summary(lognomod1)
summary(lognomod2)
summary(lognomod3)
summary(lognomod4)
anova(lognomod1)


#lognormal: check fit
lnopred <- predict(lognomod1, type = "linear")
residsLN <- (log(aids$survival) - lnopred) / lognomod1$scale
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

### 4.3

#acceleration factor: AIDS
exp(-weimod2$coefficient[2])
#odds ratio: AIDS
exp(-weimod2$coefficient[2] / weimod2$scale)

#acceleration factor: CD4 cells
exp(-weimod3$coefficient[2:4])
#odds ratio: CD4 cells
exp(-weimod3$coefficient[2:4] / weimod3$scale)

##### 5. SEMI-PARAMETRIC MODEL
