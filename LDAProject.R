setwd("C:/Users/johan/Google Drive/UNI/UPC/Lifetime Data Analysis/Lifetime-Data-Analysis")

# Objective: To analyze the survival time in patients with HIV and
# low CD4 cell count. Study the association between
# the survival time and the five variables contained in the
# data set.


###### IMPORT ######
library(survival)

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
summary(aids)
#num_aids=subset(aids, select =c(survival_times, CD4_cells))
hist(survival, main="Distribution of Survival Times", xlab="Survival Times")


##### 3.NON-PARAMETRIC ANALYIS ####

# The Surv object with right-censored data
with(aids, Surv(survival_times, cens))

### 3.1 overall survival function
svf <- survfit(Surv(survival, cens) ~ 1, aids)

par(las = 1, font = 2, font.axis = 2, font.lab = 4, bty = "l")
plot(svf, col = 1:3, lwd = 3, xlab = "Survival time [Years]",
     ylab = expression(bolditalic(hat(S)(t))), xaxs = "i", yaxs = "i",
     xlim = c(0, 25), conf.int=FALSE, main="Survival Functin AIDS patients")
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

axis(1, seq(0, 25, 5))
legend("bottomleft", levels(aids$treatment), lwd = 3, col = 1:2, bty = "n")
title("Survival functions of time to death due to vaginal cancer")

### 3.2 estimatin of median survival time
svf
#median 19.1

#distribution function
plot(svf, col = 3, xlab = "Time to Death [Years]",
     ylab = expression(bolditalic(hat(F)(t))), xaxs = "i",xlim=c(0,25),ylim=c(0,0.6),
     lty = 1, conf.int=FALSE, fun = "event", lwd = 3, yaxs = "i", bty = "n", main="Cummulative Distribution of time to death")
abline(v=19.1, col=2, lwd=2, lty=1)
abline(h=0.5, col=2, lwd=1)
