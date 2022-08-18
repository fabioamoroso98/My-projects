#Data Analysis Report

library(datasetsICR)
data(seeds)
head(seeds)

str(seeds)
summary(seeds)

library(gamlss)
library(gamlss.mx)

#Unibariate Analysis: Variety

summary(seeds$variety)
str(seeds$variety)

#Univariate Analysis: Area

summary(seeds$area)
sd(seeds$area)
var(seeds$area)
boxplot(seeds$area, main = "Boxplot of Area")
boxplot(seeds$area)$out
hist(seeds$area)

#Exponential Distribution
Area.EXP <- histDist(seeds$area, family = EXP, nbins = 30, main = "Area Exponential Distribution")
Area.EXP$df.fit
fitted(Area.EXP, "mu")[1]
logLik(Area.EXP)
AIC(Area.EXP) 
Area.EXP$sbc 

#Gamma Distribution
Area.GA <- histDist(seeds$area, family=GA, nbins = 30, main="Area Gamma distribution")
Area.GA$df.fit 
fitted(Area.GA, "mu")[1]
fitted(Area.GA, "sigma")[1] 
logLik(Area.GA)
AIC(Area.GA) 
Area.GA$sbc 

#Inverse Gaussian distribution
Area.IG <- histDist(seeds$area, family=IG, nbins = 30, main="Area Inverse Gaussian distribution")
Area.IG$df.fit
fitted(Area.IG, "mu")[1] 
fitted(Area.IG, "sigma")[1] 
logLik(Area.IG)
AIC(Area.IG) 
Area.IG$sbc  

#Log-Normal distribution
Area.LOGNO <- histDist(seeds$area, family=LOGNO, nbins = 30, main="Area Log-Normal distribution")
Area.LOGNO$df.fit
fitted(Area.LOGNO, "mu")[1] 
fitted(Area.LOGNO, "sigma")[1] 
logLik(Area.LOGNO)
AIC(Area.LOGNO) 
Area.LOGNO$sbc  

#Weibull distribution
Area.WEI <- histDist(seeds$area, family=WEI, nbins = 30, main="Area Weibull distribution")
Area.WEI$df.fit
fitted(Area.WEI, "mu")[1] 
fitted(Area.WEI, "sigma")[1] 
logLik(Area.WEI)
AIC(Area.WEI) 
Area.WEI$sbc 

#Comparing Distribution
Comparing.Area <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(Area.EXP), AIC(Area.GA), AIC(Area.IG), AIC(Area.LOGNO), AIC(Area.WEI)), SBC=c(Area.EXP$sbc, Area.GA$sbc, Area.IG$sbc, Area.LOGNO$sbc, Area.WEI$sbc), LL=c(logLik(Area.EXP), logLik(Area.GA), logLik(Area.IG), logLik(Area.LOGNO), logLik(Area.WEI)))
Comparing.Area

#Likelihood Ratio Test
LR.test(Area.EXP, Area.IG)

#Mixture Distribution with k = 2
Area.IG.2 <- gamlssMXfits(n = 5, seeds$area~1, family = IG, K = 2, data = NULL)
Area.IG.2
Area.IG.2$aic
Area.IG.2$sbc

#Estimates of "mu" and "sigma"
A.mu2.1 <- exp(Area.IG.2[["models"]][[1]][["mu.coefficients"]])
A.sigma2.1 <- exp(Area.IG.2[["models"]][[1]][["sigma.coefficients"]])

A.mu2.2 <- exp(Area.IG.2[["models"]][[2]][["mu.coefficients"]]) 
A.sigma2.2 <- exp(Area.IG.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 2
hist(seeds$area, breaks = 50,freq = FALSE, xlab = "Area", main = "Area: mixture of two Inverse Gaussian Distribution")
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),Area.IG.2[["prob"]][1]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu2.1, sigma = A.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),Area.IG.2[["prob"]][2]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu2.2, sigma = A.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),
      Area.IG.2[["prob"]][1]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu2.1, sigma = A.sigma2.1) +
        Area.IG.2[["prob"]][2]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu2.2, sigma = A.sigma2.2),
          lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
Area.IG.3 <- gamlssMXfits(n = 5, seeds$area~1, family = IG, K = 3, data = NULL)
Area.IG.3
Area.IG.3$aic
Area.IG.3$sbc

#Estimates of "mu" and "sigma"
A.mu3.1 <- exp(Area.IG.3[["models"]][[1]][["mu.coefficients"]])
A.sigma3.1 <- exp(Area.IG.3[["models"]][[1]][["sigma.coefficients"]])

A.mu3.2 <- exp(Area.IG.3[["models"]][[2]][["mu.coefficients"]]) 
A.sigma3.2 <- exp(Area.IG.3[["models"]][[2]][["sigma.coefficients"]])

A.mu3.3 <- exp(Area.IG.3[["models"]][[3]][["mu.coefficients"]]) 
A.sigma3.3 <- exp(Area.IG.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 3
hist(seeds$area, breaks = 50,freq = FALSE, xlab = "Area", main = "Area: mixture of three Inverse Gaussian Distribution")
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),Area.IG.3[["prob"]][1]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu3.1, sigma = A.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),Area.IG.3[["prob"]][2]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu3.2, sigma = A.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),Area.IG.3[["prob"]][3]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu3.3, sigma = A.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)),
      Area.IG.3[["prob"]][1]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu3.1, sigma = A.sigma3.1) +
        Area.IG.3[["prob"]][2]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu3.2, sigma = A.sigma3.2) +
          Area.IG.3[["prob"]][3]*dIG(seq(min(seeds$area),max(seeds$area),length=length(seeds$area)), mu = A.mu3.3, sigma = A.sigma3.3) ,
            lty = 1, lwd = 3, col = 1)

#Univariate Analysis: Perimeter

summary(seeds$perimeter)
sd(seeds$perimeter)
var(seeds$perimeter)
boxplot(seeds$perimeter, main = "Boxplot of Perimeter")
boxplot(seeds$perimeter)$out
hist(seeds$perimeter)

#Exponential Distribution
Perimeter.EXP <- histDist(seeds$perimeter, family = EXP, nbins = 30, main = "Perimeter Exponential Distribution")
Perimeter.EXP$df.fit
fitted(Perimeter.EXP, "mu")[1]
logLik(Perimeter.EXP)
AIC(Perimeter.EXP) 
Perimeter.EXP$sbc 

#Gamma Distribution
Perimeter.GA <- histDist(seeds$perimeter, family=GA, nbins = 30, main="Perimeter Gamma distribution")
Perimeter.GA$df.fit 
fitted(Perimeter.GA, "mu")[1]
fitted(Perimeter.GA, "sigma")[1] 
logLik(Perimeter.GA)
AIC(Perimeter.GA) 
Perimeter.GA$sbc 

#Inverse Gaussian distribution
Perimeter.IG <- histDist(seeds$perimeter, family=IG, nbins = 30, main="Perimeter Inverse Gaussian distribution")
Perimeter.IG$df.fit
fitted(Perimeter.IG, "mu")[1] 
fitted(Perimeter.IG, "sigma")[1] 
logLik(Perimeter.IG)
AIC(Perimeter.IG) 
Perimeter.IG$sbc  

#Log-Normal distribution
Perimeter.LOGNO <- histDist(seeds$perimeter, family=LOGNO, nbins = 30, main="Perimeter Log-Normal distribution")
Perimeter.LOGNO$df.fit
fitted(Perimeter.LOGNO, "mu")[1] 
fitted(Perimeter.LOGNO, "sigma")[1] 
logLik(Perimeter.LOGNO)
AIC(Perimeter.LOGNO) 
Perimeter.LOGNO$sbc  

#Weibull distribution
Perimeter.WEI <- histDist(seeds$perimeter, family=WEI, nbins = 30, main="Perimeter Weibull distribution")
Perimeter.WEI$df.fit
fitted(Perimeter.WEI, "mu")[1] 
fitted(Perimeter.WEI, "sigma")[1] 
logLik(Perimeter.WEI)
AIC(Perimeter.WEI) 
Perimeter.WEI$sbc 

#Comparing Distribution
Comparing.Perimeter <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(Perimeter.EXP), AIC(Perimeter.GA), AIC(Perimeter.IG), AIC(Perimeter.LOGNO), AIC(Perimeter.WEI)), SBC=c(Perimeter.EXP$sbc, Perimeter.GA$sbc, Perimeter.IG$sbc, Perimeter.LOGNO$sbc, Perimeter.WEI$sbc), LL=c(logLik(Perimeter.EXP), logLik(Perimeter.GA), logLik(Perimeter.IG), logLik(Perimeter.LOGNO), logLik(Perimeter.WEI)))
Comparing.Perimeter

#Likelihood Ratio Test
LR.test(Perimeter.EXP, Perimeter.IG)

#Mixture Distribution with k = 2
Perimeter.IG.2 <- gamlssMXfits(n = 5, seeds$perimeter~1, family = IG, K = 2, data = NULL)
Perimeter.IG.2
Perimeter.IG.2$aic
Perimeter.IG.2$sbc

#Estimates of "mu" and "sigma"
P.mu2.1 <- exp(Perimeter.IG.2[["models"]][[1]][["mu.coefficients"]])
P.sigma2.1 <- exp(Perimeter.IG.2[["models"]][[1]][["sigma.coefficients"]])

P.mu2.2 <- exp(Perimeter.IG.2[["models"]][[2]][["mu.coefficients"]]) 
P.sigma2.2 <- exp(Perimeter.IG.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 2
hist(seeds$perimeter, breaks = 50,freq = FALSE, xlab = "Perimeter", main = "Perimeter: mixture of two Inverse Gaussian Distribution")
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),Perimeter.IG.2[["prob"]][1]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu2.1, sigma = P.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),Perimeter.IG.2[["prob"]][2]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu2.2, sigma = P.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),
      Perimeter.IG.2[["prob"]][1]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu2.1, sigma = P.sigma2.1) +
        Perimeter.IG.2[["prob"]][2]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu2.2, sigma = P.sigma2.2),
      lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
Perimeter.IG.3 <- gamlssMXfits(n = 5, seeds$perimeter~1, family = IG, K = 3, data = NULL)
Perimeter.IG.3
Perimeter.IG.3$aic
Perimeter.IG.3$sbc

#Estimates of "mu" and "sigma"
P.mu3.1 <- exp(Perimeter.IG.3[["models"]][[1]][["mu.coefficients"]])
P.sigma3.1 <- exp(Perimeter.IG.3[["models"]][[1]][["sigma.coefficients"]])

P.mu3.2 <- exp(Perimeter.IG.3[["models"]][[2]][["mu.coefficients"]]) 
P.sigma3.2 <- exp(Perimeter.IG.3[["models"]][[2]][["sigma.coefficients"]])

P.mu3.3 <- exp(Perimeter.IG.3[["models"]][[3]][["mu.coefficients"]]) 
P.sigma3.3 <- exp(Perimeter.IG.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 3
hist(seeds$perimeter, breaks = 50,freq = FALSE, xlab = "Perimeter", main = "Perimeter: mixture of three Inverse Gaussian Distribution")
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),Perimeter.IG.3[["prob"]][1]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu3.1, sigma = P.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),Perimeter.IG.3[["prob"]][2]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu3.2, sigma = P.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),Perimeter.IG.3[["prob"]][3]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu3.3, sigma = P.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)),
      Perimeter.IG.3[["prob"]][1]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu3.1, sigma = P.sigma3.1) +
        Perimeter.IG.3[["prob"]][2]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu3.2, sigma = P.sigma3.2) +
          Perimeter.IG.3[["prob"]][3]*dIG(seq(min(seeds$perimeter),max(seeds$perimeter),length=length(seeds$perimeter)), mu = P.mu3.3, sigma = P.sigma3.3) ,
            lty = 1, lwd = 3, col = 1)

#Univariate Analysis: Compactness

summary(seeds$compactness)
sd(seeds$compactness)
var(seeds$compactness)
boxplot(seeds$compactness, main = "Boxplot of Compactness")
boxplot(seeds$compactness)$out
hist(seeds$compactness)

#Exponential Distribution
Compactness.EXP <- histDist(seeds$compactness, family = EXP, nbins = 30, main = "Compactness Exponential Distribution")
Compactness.EXP$df.fit
fitted(Compactness.EXP, "mu")[1]
logLik(Compactness.EXP)
AIC(Compactness.EXP) 
Compactness.EXP$sbc 

#Gamma Distribution
Compactness.GA <- histDist(seeds$compactness, family=GA, nbins = 30, main="Compactness Gamma distribution")
Compactness.GA$df.fit 
fitted(Compactness.GA, "mu")[1]
fitted(Compactness.GA, "sigma")[1] 
logLik(Compactness.GA)
AIC(Compactness.GA) 
Compactness.GA$sbc 

#Inverse Gaussian distribution
Compactness.IG <- histDist(seeds$compactness, family=IG, nbins = 30, main="Compactness Inverse Gaussian distribution")
Compactness.IG$df.fit
fitted(Compactness.IG, "mu")[1] 
fitted(Compactness.IG, "sigma")[1] 
logLik(Compactness.IG)
AIC(Compactness.IG) 
Compactness.IG$sbc  

#Log-Normal distribution
Compactness.LOGNO <- histDist(seeds$compactness, family=LOGNO, nbins = 30, main="Compactness Log-Normal distribution")
Compactness.LOGNO$df.fit
fitted(Compactness.LOGNO, "mu")[1] 
fitted(Compactness.LOGNO, "sigma")[1] 
logLik(Compactness.LOGNO)
AIC(Compactness.LOGNO) 
Compactness.LOGNO$sbc  

#Weibull distribution
Compactness.WEI <- histDist(seeds$compactness, family=WEI, nbins = 30, main="Compactness Weibull distribution")
Compactness.WEI$df.fit
fitted(Compactness.WEI, "mu")[1] 
fitted(Compactness.WEI, "sigma")[1] 
logLik(Compactness.WEI)
AIC(Compactness.WEI) 
Compactness.WEI$sbc 

#Comparing Distribution
Comparing.Compactness <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(Compactness.EXP), AIC(Compactness.GA), AIC(Compactness.IG), AIC(Compactness.LOGNO), AIC(Compactness.WEI)), SBC=c(Compactness.EXP$sbc, Compactness.GA$sbc, Compactness.IG$sbc, Compactness.LOGNO$sbc, Compactness.WEI$sbc), LL=c(logLik(Compactness.EXP), logLik(Compactness.GA), logLik(Compactness.IG), logLik(Compactness.LOGNO), logLik(Compactness.WEI)))
Comparing.Compactness

#Likelihood Ratio Test
LR.test(Compactness.EXP, Compactness.WEI)

#Mixture Distribution with k = 2
Compactness.WEI.2 <- gamlssMXfits(n = 5, seeds$compactness~1, family = WEI, K = 2, data = NULL)
Compactness.WEI.2
Compactness.WEI.2$aic
Compactness.WEI.2$sbc

#Estimates of "mu" and "sigma"
C.mu2.1 <- exp(Compactness.WEI.2[["models"]][[1]][["mu.coefficients"]])
C.sigma2.1 <- exp(Compactness.WEI.2[["models"]][[1]][["sigma.coefficients"]])

C.mu2.2 <- exp(Compactness.WEI.2[["models"]][[2]][["mu.coefficients"]]) 
C.sigma2.2 <- exp(Compactness.WEI.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Weibull Mixture with k = 2
hist(seeds$compactness, breaks = 50,freq = FALSE, xlab = "Compactness", main = "Compactness: mixture of two Weibull Distribution")
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),Compactness.WEI.2[["prob"]][1]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu2.1, sigma = C.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),Compactness.WEI.2[["prob"]][2]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu2.2, sigma = C.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),
      Compactness.WEI.2[["prob"]][1]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu2.1, sigma = C.sigma2.1) +
        Compactness.WEI.2[["prob"]][2]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu2.2, sigma = C.sigma2.2),
      lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
Compactness.WEI.3 <- gamlssMXfits(n = 5, seeds$compactness~1, family = WEI, K = 3, data = NULL)
Compactness.WEI.3
Compactness.WEI.3$aic
Compactness.WEI.3$sbc

#Estimates of "mu" and "sigma"
C.mu3.1 <- exp(Compactness.WEI.3[["models"]][[1]][["mu.coefficients"]])
C.sigma3.1 <- exp(Compactness.WEI.3[["models"]][[1]][["sigma.coefficients"]])

C.mu3.2 <- exp(Compactness.WEI.3[["models"]][[2]][["mu.coefficients"]]) 
C.sigma3.2 <- exp(Compactness.WEI.3[["models"]][[2]][["sigma.coefficients"]])

C.mu3.3 <- exp(Compactness.WEI.3[["models"]][[3]][["mu.coefficients"]]) 
C.sigma3.3 <- exp(Compactness.WEI.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Weibull Mixture with k = 3
hist(seeds$compactness, breaks = 50,freq = FALSE, xlab = "Compactness", main = "Compactness: mixture of three Weibull Distribution")
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),Compactness.WEI.3[["prob"]][1]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu3.1, sigma = C.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),Compactness.WEI.3[["prob"]][2]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu3.2, sigma = C.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),Compactness.WEI.3[["prob"]][3]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu3.3, sigma = C.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)),
      Compactness.WEI.3[["prob"]][1]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu3.1, sigma = C.sigma3.1) +
        Compactness.WEI.3[["prob"]][2]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu3.2, sigma = C.sigma3.2) +
        Compactness.WEI.3[["prob"]][3]*dWEI(seq(min(seeds$compactness),max(seeds$compactness),length=length(seeds$compactness)), mu = C.mu3.3, sigma = C.sigma3.3) ,
            lty = 1, lwd = 3, col = 1)

#Univariate Analysis: Length of kernel

summary(seeds$`length of kernel`)
sd(seeds$`length of kernel`)
var(seeds$`length of kernel`)
boxplot(seeds$`length of kernel`, main = "Boxplot of Length of kernel")
boxplot(seeds$`length of kernel`)$out
hist(seeds$`length of kernel`)

#Exponential Distribution
LengthKernel.EXP <- histDist(seeds$`length of kernel`, family = EXP, nbins = 30, main = "Length of kernel Exponential Distribution")
LengthKernel.EXP$df.fit
fitted(LengthKernel.EXP, "mu")[1]
logLik(LengthKernel.EXP)
AIC(LengthKernel.EXP) 
LengthKernel.EXP$sbc 

#Gamma Distribution
LengthKernel.GA <- histDist(seeds$`length of kernel`, family=GA, nbins = 30, main="Length of kernel Gamma distribution")
LengthKernel.GA$df.fit 
fitted(LengthKernel.GA, "mu")[1]
fitted(LengthKernel.GA, "sigma")[1] 
logLik(LengthKernel.GA)
AIC(LengthKernel.GA) 
LengthKernel.GA$sbc 

#Inverse Gaussian distribution
LengthKernel.IG <- histDist(seeds$`length of kernel`, family=IG, nbins = 30, main="`Length of kernel Inverse Gaussian distribution")
LengthKernel.IG$df.fit
fitted(LengthKernel.IG, "mu")[1] 
fitted(LengthKernel.IG, "sigma")[1] 
logLik(LengthKernel.IG)
AIC(LengthKernel.IG) 
LengthKernel.IG$sbc  

#Log-Normal distribution
LengthKernel.LOGNO <- histDist(seeds$`length of kernel`, family=LOGNO, nbins = 30, main="Length of kernel Log-Normal distribution")
LengthKernel.LOGNO$df.fit
fitted(LengthKernel.LOGNO, "mu")[1] 
fitted(LengthKernel.LOGNO, "sigma")[1] 
logLik(LengthKernel.LOGNO)
AIC(LengthKernel.LOGNO) 
LengthKernel.LOGNO$sbc  

#Weibull distribution
LengthKernel.WEI <- histDist(seeds$`length of kernel`, family=WEI, nbins = 30, main="Length of kernel Weibull distribution")
LengthKernel.WEI$df.fit
fitted(LengthKernel.WEI, "mu")[1] 
fitted(LengthKernel.WEI, "sigma")[1] 
logLik(LengthKernel.WEI)
AIC(LengthKernel.WEI) 
LengthKernel.WEI$sbc 

#Comparing Distribution
Comparing.LengthKernel <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(LengthKernel.EXP), AIC(LengthKernel.GA), AIC(LengthKernel.IG), AIC(LengthKernel.LOGNO), AIC(LengthKernel.WEI)), SBC=c(LengthKernel.EXP$sbc, LengthKernel.GA$sbc, LengthKernel.IG$sbc, LengthKernel.LOGNO$sbc, LengthKernel.WEI$sbc), LL=c(logLik(LengthKernel.EXP), logLik(LengthKernel.GA), logLik(LengthKernel.IG), logLik(LengthKernel.LOGNO), logLik(LengthKernel.WEI)))
Comparing.LengthKernel

#Likelihood Ratio Test
LR.test(LengthKernel.EXP, LengthKernel.IG)

#Mixture Distribution with k = 2
LengthKernel.IG.2 <- gamlssMXfits(n = 5, seeds$`length of kernel`~1, family = IG, K = 2, data = NULL)
LengthKernel.IG.2
LengthKernel.IG.2$aic
LengthKernel.IG.2$sbc

#Estimates of "mu" and "sigma"
LK.mu2.1 <- exp(LengthKernel.IG.2[["models"]][[1]][["mu.coefficients"]])
LK.sigma2.1 <- exp(LengthKernel.IG.2[["models"]][[1]][["sigma.coefficients"]])

LK.mu2.2 <- exp(LengthKernel.IG.2[["models"]][[2]][["mu.coefficients"]]) 
LK.sigma2.2 <- exp(LengthKernel.IG.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 2
hist(seeds$`length of kernel`, breaks = 50,freq = FALSE, xlab = "Length of kernel", main = "Length of kernel: mixture of two Inverse Gaussian Distribution")
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),LengthKernel.IG.2[["prob"]][1]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu2.1, sigma = LK.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),LengthKernel.IG.2[["prob"]][2]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu2.2, sigma = LK.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),
      LengthKernel.IG.2[["prob"]][1]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu2.1, sigma = LK.sigma2.1) +
        LengthKernel.IG.2[["prob"]][2]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu2.2, sigma = LK.sigma2.2),
      lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
LengthKernel.IG.3 <- gamlssMXfits(n = 5, seeds$`length of kernel`~1, family = IG, K = 3, data = NULL)
LengthKernel.IG.3
LengthKernel.IG.3$aic
LengthKernel.IG.3$sbc

#Estimates of "mu" and "sigma"
LK.mu3.1 <- exp(LengthKernel.IG.3[["models"]][[1]][["mu.coefficients"]])
LK.sigma3.1 <- exp(LengthKernel.IG.3[["models"]][[1]][["sigma.coefficients"]])

LK.mu3.2 <- exp(LengthKernel.IG.3[["models"]][[2]][["mu.coefficients"]]) 
LK.sigma3.2 <- exp(LengthKernel.IG.3[["models"]][[2]][["sigma.coefficients"]])

LK.mu3.3 <- exp(LengthKernel.IG.3[["models"]][[3]][["mu.coefficients"]]) 
LK.sigma3.3 <- exp(LengthKernel.IG.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 3
hist(seeds$`length of kernel`, breaks = 50,freq = FALSE, xlab = "Length of kernel", main = "Length of kernel: mixture of three Inverse Gaussian Distribution")
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),LengthKernel.IG.3[["prob"]][1]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu3.1, sigma = LK.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),LengthKernel.IG.3[["prob"]][2]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu3.2, sigma = LK.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),LengthKernel.IG.3[["prob"]][3]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu3.3, sigma = LK.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)),
      LengthKernel.IG.3[["prob"]][1]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu3.1, sigma = LK.sigma3.1) +
        LengthKernel.IG.3[["prob"]][2]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu3.2, sigma = LK.sigma3.2) +
          LengthKernel.IG.3[["prob"]][3]*dIG(seq(min(seeds$`length of kernel`),max(seeds$`length of kernel`),length=length(seeds$`length of kernel`)), mu = LK.mu3.3, sigma = LK.sigma3.3) ,
            lty = 1, lwd = 3, col = 1)

#Univariate Analysis: Width of kernel

summary(seeds$`width of kernel`)
sd(seeds$`width of kernel`)
var(seeds$`width of kernel`)
boxplot(seeds$`width of kernel`, main = "Boxplot of Width of kernel")
boxplot(seeds$`width of kernel`)$out
hist(seeds$`width of kernel`)

#Exponential Distribution
WidthKernel.EXP <- histDist(seeds$`width of kernel`, family = EXP, nbins = 30, main = "Width of kernel Exponential Distribution")
WidthKernel.EXP$df.fit
fitted(WidthKernel.EXP, "mu")[1]
logLik(WidthKernel.EXP)
AIC(WidthKernel.EXP) 
WidthKernel.EXP$sbc 

#Gamma Distribution
WidthKernel.GA <- histDist(seeds$`width of kernel`, family=GA, nbins = 30, main="`Width of kernel Gamma distribution")
WidthKernel.GA$df.fit 
fitted(WidthKernel.GA, "mu")[1]
fitted(WidthKernel.GA, "sigma")[1] 
logLik(WidthKernel.GA)
AIC(WidthKernel.GA) 
WidthKernel.GA$sbc 

#Inverse Gaussian distribution
WidthKernel.IG <- histDist(seeds$`width of kernel`, family=IG, nbins = 30, main="`Width of kernel Inverse Gaussian distribution")
WidthKernel.IG$df.fit
fitted(WidthKernel.IG, "mu")[1] 
fitted(WidthKernel.IG, "sigma")[1] 
logLik(WidthKernel.IG)
AIC(WidthKernel.IG) 
WidthKernel.IG$sbc  

#Log-Normal distribution
WidthKernel.LOGNO <- histDist(seeds$`width of kernel`, family=LOGNO, nbins = 30, main="Width of kernel Log-Normal distribution")
WidthKernel.LOGNO$df.fit
fitted(WidthKernel.LOGNO, "mu")[1] 
fitted(WidthKernel.LOGNO, "sigma")[1] 
logLik(WidthKernel.LOGNO)
AIC(WidthKernel.LOGNO) 
WidthKernel.LOGNO$sbc  

#Weibull distribution
WidthKernel.WEI <- histDist(seeds$`width of kernel`, family=WEI, nbins = 30, main="Width of kernel Weibull distribution")
WidthKernel.WEI$df.fit
fitted(WidthKernel.WEI, "mu")[1] 
fitted(WidthKernel.WEI, "sigma")[1] 
logLik(WidthKernel.WEI)
AIC(WidthKernel.WEI) 
WidthKernel.WEI$sbc 

#Comparing Distribution
Comparing.WidthKernel <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(WidthKernel.EXP), AIC(WidthKernel.GA), AIC(WidthKernel.IG), AIC(WidthKernel.LOGNO), AIC(WidthKernel.WEI)), SBC=c(WidthKernel.EXP$sbc, WidthKernel.GA$sbc, WidthKernel.IG$sbc, WidthKernel.LOGNO$sbc, WidthKernel.WEI$sbc), LL=c(logLik(WidthKernel.EXP), logLik(WidthKernel.GA), logLik(WidthKernel.IG), logLik(WidthKernel.LOGNO), logLik(WidthKernel.WEI)))
Comparing.WidthKernel

#Likelihood Ratio Test
LR.test(WidthKernel.EXP, WidthKernel.IG)

#Mixture Distribution with k = 2
WidthKernel.IG.2 <- gamlssMXfits(n = 5, seeds$`width of kernel`~1, family = IG, K = 2, data = NULL)
WidthKernel.IG.2
WidthKernel.IG.2$aic
WidthKernel.IG.2$sbc

#Estimates of "mu" and "sigma"
WK.mu2.1 <- exp(WidthKernel.IG.2[["models"]][[1]][["mu.coefficients"]])
WK.sigma2.1 <- exp(WidthKernel.IG.2[["models"]][[1]][["sigma.coefficients"]])

WK.mu2.2 <- exp(WidthKernel.IG.2[["models"]][[2]][["mu.coefficients"]]) 
WK.sigma2.2 <- exp(WidthKernel.IG.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 2
hist(seeds$`width of kernel`, breaks = 50,freq = FALSE, xlab = "Width of kernel", main = "Width of kernel: mixture of two Inverse Gaussian Distribution")
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),WidthKernel.IG.2[["prob"]][1]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu2.1, sigma = WK.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),WidthKernel.IG.2[["prob"]][2]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu2.2, sigma = WK.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),
      WidthKernel.IG.2[["prob"]][1]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu2.1, sigma = WK.sigma2.1) +
        WidthKernel.IG.2[["prob"]][2]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu2.2, sigma = WK.sigma2.2),
      lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
WidthKernel.IG.3 <- gamlssMXfits(n = 5, seeds$`width of kernel`~1, family = IG, K = 3, data = NULL)
WidthKernel.IG.3
WidthKernel.IG.3$aic
WidthKernel.IG.3$sbc

#Estimates of "mu" and "sigma"
WK.mu3.1 <- exp(WidthKernel.IG.3[["models"]][[1]][["mu.coefficients"]])
WK.sigma3.1 <- exp(WidthKernel.IG.3[["models"]][[1]][["sigma.coefficients"]])

WK.mu3.2 <- exp(WidthKernel.IG.3[["models"]][[2]][["mu.coefficients"]]) 
WK.sigma3.2 <- exp(WidthKernel.IG.3[["models"]][[2]][["sigma.coefficients"]])

WK.mu3.3 <- exp(WidthKernel.IG.3[["models"]][[3]][["mu.coefficients"]]) 
WK.sigma3.3 <- exp(WidthKernel.IG.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 3
hist(seeds$`width of kernel`, breaks = 50,freq = FALSE, xlab = "Width of kernel", main = "Width of kernel: mixture of three Inverse Gaussian Distribution")
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),WidthKernel.IG.3[["prob"]][1]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu3.1, sigma = WK.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),WidthKernel.IG.3[["prob"]][2]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu3.2, sigma = WK.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),WidthKernel.IG.3[["prob"]][3]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu3.3, sigma = WK.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)),
      WidthKernel.IG.3[["prob"]][1]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu3.1, sigma = WK.sigma3.1) +
        WidthKernel.IG.3[["prob"]][2]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu3.2, sigma = WK.sigma3.2) +
        WidthKernel.IG.3[["prob"]][3]*dIG(seq(min(seeds$`width of kernel`),max(seeds$`width of kernel`),length=length(seeds$`width of kernel`)), mu = WK.mu3.3, sigma = WK.sigma3.3) ,
      lty = 1, lwd = 3, col = 1)

#Univariate Analysis: Asymmetry coefficient

summary(seeds$`asymmetry coefficient`)
sd(seeds$`asymmetry coefficient`)
var(seeds$`asymmetry coefficient`)
boxplot(seeds$`asymmetry coefficient`, main = "Boxplot of Asymmetry coefficient")
boxplot(seeds$`asymmetry coefficient`)$out
hist(seeds$`asymmetry coefficient`)

#Exponential Distribution
AsymmetryCoefficient.EXP <- histDist(seeds$`asymmetry coefficient`, family = EXP, nbins = 30, main = "Asymmetry coefficient Exponential Distribution")
AsymmetryCoefficient.EXP$df.fit
fitted(AsymmetryCoefficient.EXP, "mu")[1]
logLik(AsymmetryCoefficient.EXP)
AIC(AsymmetryCoefficient.EXP) 
AsymmetryCoefficient.EXP$sbc 

#Gamma Distribution
AsymmetryCoefficient.GA <- histDist(seeds$`asymmetry coefficient`, family=GA, nbins = 30, main="Asymmetry coefficient Gamma distribution")
AsymmetryCoefficient.GA$df.fit 
fitted(AsymmetryCoefficient.GA, "mu")[1]
fitted(AsymmetryCoefficient.GA, "sigma")[1] 
logLik(AsymmetryCoefficient.GA)
AIC(AsymmetryCoefficient.GA) 
AsymmetryCoefficient.GA$sbc 

#Inverse Gaussian distribution
AsymmetryCoefficient.IG <- histDist(seeds$`asymmetry coefficient`, family=IG, nbins = 30, main="Asymmetry coefficient Inverse Gaussian distribution")
AsymmetryCoefficient.IG$df.fit
fitted(AsymmetryCoefficient.IG, "mu")[1] 
fitted(AsymmetryCoefficient.IG, "sigma")[1] 
logLik(AsymmetryCoefficient.IG)
AIC(AsymmetryCoefficient.IG) 
AsymmetryCoefficient.IG$sbc  

#Log-Normal distribution
AsymmetryCoefficient.LOGNO <- histDist(seeds$`asymmetry coefficient`, family=LOGNO, nbins = 30, main="Asymmetry coefficient Log-Normal distribution")
AsymmetryCoefficient.LOGNO$df.fit
fitted(AsymmetryCoefficient.LOGNO, "mu")[1] 
fitted(AsymmetryCoefficient.LOGNO, "sigma")[1] 
logLik(AsymmetryCoefficient.LOGNO)
AIC(AsymmetryCoefficient.LOGNO) 
AsymmetryCoefficient.LOGNO$sbc  

#Weibull distribution
AsymmetryCoefficient.WEI <- histDist(seeds$`asymmetry coefficient`, family=WEI, nbins = 30, main="Asymmetry coefficient Weibull distribution")
AsymmetryCoefficient.WEI$df.fit
fitted(AsymmetryCoefficient.WEI, "mu")[1] 
fitted(AsymmetryCoefficient.WEI, "sigma")[1] 
logLik(AsymmetryCoefficient.WEI)
AIC(AsymmetryCoefficient.WEI) 
AsymmetryCoefficient.WEI$sbc 

#Comparing Distribution
Comparing.AsymmetryCoefficient <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(AsymmetryCoefficient.EXP), AIC(AsymmetryCoefficient.GA), AIC(AsymmetryCoefficient.IG), AIC(AsymmetryCoefficient.LOGNO), AIC(AsymmetryCoefficient.WEI)), SBC=c(AsymmetryCoefficient.EXP$sbc, AsymmetryCoefficient.GA$sbc, AsymmetryCoefficient.IG$sbc, AsymmetryCoefficient.LOGNO$sbc, AsymmetryCoefficient.WEI$sbc), LL=c(logLik(AsymmetryCoefficient.EXP), logLik(AsymmetryCoefficient.GA), logLik(AsymmetryCoefficient.IG), logLik(AsymmetryCoefficient.LOGNO), logLik(AsymmetryCoefficient.WEI)))
Comparing.AsymmetryCoefficient

#Likelihood Ratio Test
LR.test(AsymmetryCoefficient.EXP, AsymmetryCoefficient.WEI)

#Mixture Distribution with k = 2
AsymmetryCoefficient.WEI.2 <- gamlssMXfits(n = 5, seeds$`asymmetry coefficient`~1, family = WEI, K = 2, data = NULL)
AsymmetryCoefficient.WEI.2
AsymmetryCoefficient.WEI.2$aic
AsymmetryCoefficient.WEI.2$sbc

#Estimates of "mu" and "sigma"
AC.mu2.1 <- exp(AsymmetryCoefficient.WEI.2[["models"]][[1]][["mu.coefficients"]])
AC.sigma2.1 <- exp(AsymmetryCoefficient.WEI.2[["models"]][[1]][["sigma.coefficients"]])

AC.mu2.2 <- exp(AsymmetryCoefficient.WEI.2[["models"]][[2]][["mu.coefficients"]]) 
AC.sigma2.2 <- exp(AsymmetryCoefficient.WEI.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Weibull Mixture with k = 2
hist(seeds$`asymmetry coefficient`, breaks = 50,freq = FALSE, xlab = "Asymmetry Coefficient", main = "Asymmetry Coefficient: mixture of two Weibull Distribution")
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),AsymmetryCoefficient.WEI.2[["prob"]][1]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu2.1, sigma = AC.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),AsymmetryCoefficient.WEI.2[["prob"]][2]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu2.2, sigma = AC.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),
      AsymmetryCoefficient.WEI.2[["prob"]][1]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu2.1, sigma = AC.sigma2.1) +
        AsymmetryCoefficient.WEI.2[["prob"]][2]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu2.2, sigma = AC.sigma2.2),
      lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
AsymmetryCoefficient.WEI.3 <- gamlssMXfits(n = 5, seeds$`asymmetry coefficient`~1, family = WEI, K = 3, data = NULL)
AsymmetryCoefficient.WEI.3
AsymmetryCoefficient.WEI.3$aic
AsymmetryCoefficient.WEI.3$sbc

#Estimates of "mu" and "sigma"
AC.mu3.1 <- exp(AsymmetryCoefficient.WEI.3[["models"]][[1]][["mu.coefficients"]])
AC.sigma3.1 <- exp(AsymmetryCoefficient.WEI.3[["models"]][[1]][["sigma.coefficients"]])

AC.mu3.2 <- exp(AsymmetryCoefficient.WEI.3[["models"]][[2]][["mu.coefficients"]]) 
AC.sigma3.2 <- exp(AsymmetryCoefficient.WEI.3[["models"]][[2]][["sigma.coefficients"]])

AC.mu3.3 <- exp(AsymmetryCoefficient.WEI.3[["models"]][[3]][["mu.coefficients"]]) 
AC.sigma3.3 <- exp(AsymmetryCoefficient.WEI.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Weibull Mixture with k = 3
hist(seeds$`asymmetry coefficient`, breaks = 50,freq = FALSE, xlab = "Asymmetry coefficient", main = "Asymmetry coefficient: mixture of three Weibull Distribution")
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),AsymmetryCoefficient.WEI.3[["prob"]][1]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu3.1, sigma = AC.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),AsymmetryCoefficient.WEI.3[["prob"]][2]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu3.2, sigma = AC.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),AsymmetryCoefficient.WEI.3[["prob"]][3]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu3.3, sigma = AC.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)),
      AsymmetryCoefficient.WEI.3[["prob"]][1]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu3.1, sigma = AC.sigma3.1) +
        AsymmetryCoefficient.WEI.3[["prob"]][2]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu3.2, sigma = AC.sigma3.2) +
        AsymmetryCoefficient.WEI.3[["prob"]][3]*dWEI(seq(min(seeds$`asymmetry coefficient`),max(seeds$`asymmetry coefficient`),length=length(seeds$`asymmetry coefficient`)), mu = AC.mu3.3, sigma = AC.sigma3.3) ,
      lty = 1, lwd = 3, col = 1)

#Univariate Analysis: Length of kernel groove

summary(seeds$`length of kernel groove`)
sd(seeds$`length of kernel groove`)
var(seeds$`length of kernel groove`)
boxplot(seeds$`length of kernel groove`, main = "Boxplot of Length of kernel groove")
boxplot(seeds$`length of kernel groove`)$out
hist(seeds$`length of kernel groove`)

#Exponential Distribution
KernelGroove.EXP <- histDist(seeds$`length of kernel groove`, family = EXP, nbins = 30, main = "Length of kernel groove Exponential Distribution")
KernelGroove.EXP$df.fit
fitted(KernelGroove.EXP, "mu")[1]
logLik(KernelGroove.EXP)
AIC(KernelGroove.EXP) 
KernelGroove.EXP$sbc 

#Gamma Distribution
KernelGroove.GA <- histDist(seeds$`length of kernel groove`, family=GA, nbins = 30, main="Length of kernel groove Gamma distribution")
KernelGroove.GA$df.fit 
fitted(KernelGroove.GA, "mu")[1]
fitted(KernelGroove.GA, "sigma")[1] 
logLik(KernelGroove.GA)
AIC(KernelGroove.GA) 
KernelGroove.GA$sbc 

#Inverse Gaussian distribution
KernelGroove.IG <- histDist(seeds$`length of kernel groove`, family=IG, nbins = 30, main="Length of kernel groove Inverse Gaussian distribution")
KernelGroove.IG$df.fit
fitted(KernelGroove.IG, "mu")[1] 
fitted(KernelGroove.IG, "sigma")[1] 
logLik(KernelGroove.IG)
AIC(KernelGroove.IG) 
KernelGroove.IG$sbc  

#Log-Normal distribution
KernelGroove.LOGNO <- histDist(seeds$`length of kernel groove`, family=LOGNO, nbins = 30, main="Length of kernel groove Log-Normal distribution")
KernelGroove.LOGNO$df.fit
fitted(KernelGroove.LOGNO, "mu")[1] 
fitted(KernelGroove.LOGNO, "sigma")[1] 
logLik(KernelGroove.LOGNO)
AIC(KernelGroove.LOGNO) 
KernelGroove.LOGNO$sbc  

#Weibull distribution
KernelGroove.WEI <- histDist(seeds$`length of kernel groove`, family=WEI, nbins = 30, main="Length of kernel groovE Weibull distribution")
KernelGroove.WEI$df.fit
fitted(KernelGroove.WEI, "mu")[1] 
fitted(KernelGroove.WEI, "sigma")[1] 
logLik(KernelGroove.WEI)
AIC(KernelGroove.WEI) 
KernelGroove.WEI$sbc 

#Comparing Distribution
Comparing.KernelGroove <- data.frame(row.names = c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", "Weibull"), AIC=c(AIC(KernelGroove.EXP), AIC(KernelGroove.GA), AIC(KernelGroove.IG), AIC(KernelGroove.LOGNO), AIC(KernelGroove.WEI)), SBC=c(KernelGroove.EXP$sbc, KernelGroove.GA$sbc, KernelGroove.IG$sbc, KernelGroove.LOGNO$sbc, KernelGroove.WEI$sbc), LL=c(logLik(KernelGroove.EXP), logLik(KernelGroove.GA), logLik(KernelGroove.IG), logLik(KernelGroove.LOGNO), logLik(KernelGroove.WEI)))
Comparing.KernelGroove

#Likelihood Ratio Test
LR.test(KernelGroove.EXP, KernelGroove.IG)

#Mixture Distribution with k = 2
KernelGroove.IG.2 <- gamlssMXfits(n = 5, seeds$`length of kernel groove`~1, family = IG, K = 2, data = NULL)
KernelGroove.IG.2
KernelGroove.IG.2$aic
KernelGroove.IG.2$sbc

#Estimates of "mu" and "sigma"
KG.mu2.1 <- exp(KernelGroove.IG.2[["models"]][[1]][["mu.coefficients"]])
KG.sigma2.1 <- exp(KernelGroove.IG.2[["models"]][[1]][["sigma.coefficients"]])

KG.mu2.2 <- exp(KernelGroove.IG.2[["models"]][[2]][["mu.coefficients"]]) 
KG.sigma2.2 <- exp(KernelGroove.IG.2[["models"]][[2]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 2
hist(seeds$`length of kernel groove`, breaks = 50,freq = FALSE, xlab = "Length of kernel groove", main = "Length of kernel groove: mixture of two Inverse Gaussian Distribution")
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),KernelGroove.IG.2[["prob"]][1]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu2.1, sigma = KG.sigma2.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),KernelGroove.IG.2[["prob"]][2]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu2.2, sigma = KG.sigma2.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),
      KernelGroove.IG.2[["prob"]][1]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu2.1, sigma = KG.sigma2.1) +
        KernelGroove.IG.2[["prob"]][2]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu2.2, sigma = KG.sigma2.2),
          lty = 1, lwd = 3, col = 1)

#Mixture Distribution with k = 3
KernelGroove.IG.3 <- gamlssMXfits(n = 5, seeds$`length of kernel groove`~1, family = IG, K = 3, data = NULL)
KernelGroove.IG.3
KernelGroove.IG.3$aic
KernelGroove.IG.3$sbc

#Estimates of "mu" and "sigma"
KG.mu3.1 <- exp(KernelGroove.IG.3[["models"]][[1]][["mu.coefficients"]])
KG.sigma3.1 <- exp(KernelGroove.IG.3[["models"]][[1]][["sigma.coefficients"]])

KG.mu3.2 <- exp(KernelGroove.IG.3[["models"]][[2]][["mu.coefficients"]]) 
KG.sigma3.2 <- exp(KernelGroove.IG.3[["models"]][[2]][["sigma.coefficients"]])

KG.mu3.3 <- exp(KernelGroove.IG.3[["models"]][[3]][["mu.coefficients"]]) 
KG.sigma3.3 <- exp(KernelGroove.IG.3[["models"]][[3]][["sigma.coefficients"]])

#Histogram of the Inverse Gaussian Mixture with k = 3
hist(seeds$`length of kernel groove`, breaks = 50,freq = FALSE, xlab = "Length of kernel groove", main = "Length of kernel groove: mixture of three Inverse Gaussian Distribution")
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),KernelGroove.IG.3[["prob"]][1]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu3.1, sigma = KG.sigma3.1),lty=2,lwd=3,col=2)
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),KernelGroove.IG.3[["prob"]][2]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu3.2, sigma = KG.sigma3.2),lty=2,lwd=3,col=3)
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),KernelGroove.IG.3[["prob"]][3]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu3.3, sigma = KG.sigma3.3),lty=2,lwd=3,col=4)
lines(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)),
      KernelGroove.IG.3[["prob"]][1]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu3.1, sigma = KG.sigma3.1) +
        KernelGroove.IG.3[["prob"]][2]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu3.2, sigma = KG.sigma3.2) +
        KernelGroove.IG.3[["prob"]][3]*dIG(seq(min(seeds$`length of kernel groove`),max(seeds$`length of kernel groove`),length=length(seeds$`length of kernel groove`)), mu = KG.mu3.3, sigma = KG.sigma3.3) ,
      lty = 1, lwd = 3, col = 1)

#Principal Component Analysis
numerical.seeds <- seeds[,-8]
head(numerical.seeds)
library(tidyverse)
library(gridExtra)
library(ggplot2)

pairs(numerical.seeds, gap=0, pch=16, asp=1)
cor(numerical.seeds)

scaled_seeds <- scale(numerical.seeds)
head(scaled_seeds)

cov_seeds <- cov(scaled_seeds)
head(cov_seeds)

eigen_seeds <- eigen(cov_seeds)
head(eigen_seeds)
str(eigen_seeds)

loadings_seeds <- eigen_seeds$vectors
length(loadings_seeds[loadings_seeds>0])
length(loadings_seeds[loadings_seeds<0])

rownames(loadings_seeds) <- colnames(scaled_seeds)
colnames(loadings_seeds) <- c("PC1", "PC2", "PC3", "PC4","PC5","PC6","PC7")

head(loadings_seeds)

PVE <- eigen_seeds$values/sum(eigen_seeds$values)
PVE
cumsum(PVE)
plot.PVE <- plot(PVE, main = "Scree Plot", xlab = "Principal Component", ylab = "Proportion of Variance Explained", ylim = c(0,1), type = "b")
plot.cumsumPVE <- plot(cumsum(PVE),main = "Cumulative Scree Plot",  xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", ylim = c(0,1), type = "b")

eigen_seeds$values[eigen_seeds$values>1]

PC1 <- scaled_seeds %*% loadings_seeds[,1]
PC2 <- scaled_seeds %*% loadings_seeds[,2]

PC <- data.frame(PC1, PC2)

ggplot(PC, aes(PC1, PC2)) +
 modelr::geom_ref_line(h = 0) +
 modelr::geom_ref_line(v = 0) +
 geom_point(pch=16) +
 xlab("PC 1") +
 ylab("PC 2")

PCA.seeds <- ?prcomp(numerical.seeds, scale = T)
head(PCA.seeds)

biplot(PCA.seeds, cex = 0.6)
abline(h=0)
abline(v=0)

#Cluster Analysis

head(scaled_seeds)
library(factoextra)
library(cluster)
library(NbClust)

#Cluster Tendency
library(clustertend)
hopkins(scaled_seeds, n = nrow(scaled_seeds) -1)

fviz_dist(dist(scaled_seeds), show_labels = FALSE) +
  labs(title = "VAT Seeds Dataset")

#Choose the number of clusters: silhouette, wss, gap stat

fviz_nbclust(scaled_seeds, kmeans, method = "silhouette") + 
  labs(subtitle = "K-means: Silhouette method")

fviz_nbclust(scaled_seeds, kmeans, method = "wss") +
  labs(subtitle = "K-means: Elbow method")

fviz_nbclust(scaled_seeds, kmeans, method = "gap_stat") +
  labs(subtitle = "K-means: Gap Stat method")

nk.km <- NbClust(scaled_seeds, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 5, method = "kmeans")

fviz_nbclust(nk.km) +
  labs(subtitle = "K-means")


#k=3 according to Elbow and Gap State methods for K-means

fviz_nbclust(scaled_seeds, pam, method = "silhouette") +
  labs(subtitle = "K-medoids: Silhouette method")

fviz_nbclust(scaled_seeds, pam, method = "wss") +
  labs(subtitle = "K-medoids: Elbow method")

fviz_nbclust(scaled_seeds, pam, method = "gap_stat") +
  labs(subtitle = "K-medoids: Gap Stat method")

#k=3 according to Elbow and Gap State methods for K-medoids

fviz_nbclust(scaled_seeds, hcut, method = "silhouette") +
  labs(subtitle = "Hierarchical Clustering: Silhouette method")

fviz_nbclust(scaled_seeds, hcut, method = "wss") +
  labs(subtitle = "Hierarchical Clustering: Elbow method")

fviz_nbclust(scaled_seeds, hcut, method = "gap_stat") +
  labs(subtitle = "Hierarchical Clustering: Gap Stat method")

nk.hc <- NbClust(scaled_seeds, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 5, method = "ward.D2")

fviz_nbclust(nk.hc) +
  labs(subtitle = "Hierarchical Clustering")

#k=3 according to Elbow and Gap State methods for Hierarchical Clustering

#Generating K-means Clusters

set.seed(123)
km.res2 <- kmeans(scaled_seeds, 2, nstart = 25)
km.res3 <- kmeans(scaled_seeds, 3, nstart = 25)

print(km.res2)
print(km.res3)

aggregate(numerical.seeds, by=list(cluster=km.res2$cluster), mean)
aggregate(numerical.seeds, by=list(cluster=km.res3$cluster), mean)

seeds.km.2 <- cbind(seeds, cluster = km.res2$cluster)
head(seeds.km.2)

seeds.km.3 <- cbind(seeds, cluster = km.res3$cluster)
head(seeds.km.3)

km.res2$size
km.res3$size

km.res2$centers
km.res3$centers

table(seeds$variety, km.res2$cluster)
table(seeds$variety, km.res3$cluster)

#Plotting clusters in the original space

km.cl.2 <- km.res2$cluster
pairs(scaled_seeds, gap=0, pch=km.cl.2, col=c("red", "navy")[km.cl.2])

km.cl.3 <- km.res3$cluster
pairs(scaled_seeds, gap=0, pch=km.cl.3, col=c("red", "navy", "orange")[km.cl.3])

#Plotting cluster on the first 2PCs

fviz_cluster(km.res2, 
             data = scaled_seeds,
             palette = c("red", "navy"),
             ellipse.type = "euclid",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal()
)

fviz_cluster(km.res3, 
             data = scaled_seeds,
             palette = c("red", "navy", "orange"),
             ellipse.type = "euclid",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal()
)

#Generating K-medoids clusters

pam.res2 <- pam(scaled_seeds, 2)
print(pam.res2)

pam.res3 <- pam(scaled_seeds, 3)
print(pam.res3)

seeds.pam.2 <- cbind(seeds, cluster = pam.res2$cluster)
head(seeds.pam.2, n = 8)

seeds.pam.3 <- cbind(seeds, cluster = pam.res3$cluster)
head(seeds.pam.3, n = 8)

pam.res2$medoids
pam.res3$medoids

table(seeds$variety, pam.res2$cluster)
table(seeds$variety, pam.res3$cluster)

#Plotting cluster in the original space

pam.cl.2 <- pam.res2$clustering
pairs(scaled_seeds, gap=0, pch=pam.cl.2, col=c("red", "navy")[pam.cl.2])

pam.cl.3 <- pam.res3$clustering
pairs(scaled_seeds, gap=0, pch=pam.cl.3, col=c("red", "navy", "orange")[pam.cl.3])

#Plotting clusters on the first 2PCs

fviz_cluster(pam.res2,
             palette = c("red", "navy"), 
             ellipse.type = "t", 
             repel = TRUE, 
             ggtheme = theme_minimal()
)

fviz_cluster(pam.res3,
             palette = c("red", "navy", "orange"), 
             ellipse.type = "t", 
             repel = TRUE, 
             ggtheme = theme_minimal()
)

#Generate Hierarchical Clusters

#Dissimilarity Matrix

diss.seeds <- dist(scaled_seeds, method = "euclidean")

#Find the best Cophenetic Distance

ward.d2.seeds <- hclust(d = diss.seeds, method = "ward.D2")

fviz_dend(ward.d2.seeds, show_labels = F)

coph.ward.seeds <- cophenetic(ward.d2.seeds)

cor(diss.seeds, coph.ward.seeds)

comp.seeds <- hclust(d = diss.seeds, method = "complete")

fviz_dend(comp.seeds, show_labels = F)

coph.comp.seeds <- cophenetic(comp.seeds)

cor(diss.seeds, coph.comp.seeds)

single.seeds <- hclust(d = diss.seeds, method = "single")

fviz_dend(single.seeds, show_labels = F)

coph.single.seeds <- cophenetic(single.seeds)

cor(diss.seeds, coph.single.seeds)

avg.seeds <- hclust(d = diss.seeds, method = "average")

fviz_dend(avg.seeds, show_labels = F)

coph.avg.seeds <- cophenetic(avg.seeds)

cor(diss.seeds, coph.avg.seeds)

#The best method is the Ward method

#Hierarchical Clustering with k = 2

hc.seeds.2 <- cutree(ward.d2.seeds, k = 2)

table(hc.seeds.2)

fviz_dend(ward.d2.seeds, 
          k = 2, 
          show_labels = F,
          k_colors = c("red", "navy"),
          rect = T)

table(seeds$variety, hc.seeds.2)

#Plotting Clusters in the original space

pairs(scaled_seeds, gap=0, pch=hc.seeds.2, col=c("red", "navy")[hc.seeds.2])

#Plotting Clusters in the first two PCs

fviz_cluster(list(data = scaled_seeds, cluster = hc.seeds.2),
             palette = c("red", "navy"), 
             ellipse.type = "convex", 
             repel = TRUE, 
             ggtheme = theme_minimal()
)

#Hierarchical Clustering with k = 3

hc.seeds.3 <- cutree(ward.d2.seeds, k = 3)

table(hc.seeds.3)

fviz_dend(ward.d2.seeds, 
          k = 3, 
          show_labels = F,
          k_colors = c("red", "navy", "orange"),
          rect = T)

table(seeds$variety, hc.seeds.3)

#Plotting Clusters in the original space

pairs(scaled_seeds, gap=0, pch=hc.seeds.3, col=c("red", "navy", "orange")[hc.seeds.3])

#Plotting Clusters in the first two PCs

fviz_cluster(list(data = scaled_seeds, cluster = hc.seeds.3),
             palette = c("red", "navy", "orange"), 
             ellipse.type = "convex", 
             repel = TRUE, 
             ggtheme = theme_minimal()
)
#Model Based Clustering

library(mclust)

#Fit the Parsimonious Gaussian Mixture to our data

model.seeds <- Mclust(scaled_seeds, G = 1:9, modelNames = NULL)

summary(model.seeds$BIC)

summary(model.seeds)

table(seeds$variety, model.seeds$classification)

adjustedRandIndex(seeds$variety, model.seeds$classification)

# BIC values used for choosing the number of clusters

fviz_mclust(model.seeds, "BIC", palette = "jco")

#Plotting the clusters in the original space

pairs(scaled_seeds, gap=0, pch = 16, col = model.seeds$classification)

#Plotting the clusters in the first two PCs

fviz_mclust(model.seeds, "classification", geom = "point", pointsize = 1.5, palette = "jco")

#Classification uncertainty

fviz_mclust(model.seeds, "uncertainty", palette = "jco")

#Internal Cluster Validation

library(clValid)

clmethods <- c("kmeans", "pam", "hierarchical", "model")

ClusterValidation <- clValid(scaled_seeds, nClust = 2:5, clMethods = clmethods, validation = "internal")

summary(ClusterValidation)

#External Cluster Validation

library(fpc)

variety <- as.numeric(seeds$variety)

km.stats2 <- cluster.stats(d = dist(numerical.seeds), variety, km.res2$cluster)
km.stats3 <- cluster.stats(d = dist(numerical.seeds), variety, km.res3$cluster)

km.stats2$corrected.rand
km.stats2$vi

km.stats3$corrected.rand
km.stats3$vi

pam.stats2 <- cluster.stats(d = dist(numerical.seeds), variety, pam.res2$cluster)
pam.stats3 <- cluster.stats(d = dist(numerical.seeds), variety, pam.res3$cluster)

pam.stats2$corrected.rand
pam.stats2$vi

pam.stats3$corrected.rand
pam.stats3$vi

hc.stats2 <- cluster.stats(d = dist(numerical.seeds), variety, hc.seeds.2)
hc.stats3 <- cluster.stats(d = dist(numerical.seeds), variety, hc.seeds.3)

hc.stats2$corrected.rand
hc.stats2$vi

hc.stats3$corrected.rand
hc.stats3$vi

external.validation <- data.frame(c(km.stats2$corrected.rand, km.stats2$vi), c(km.stats3$corrected.rand, km.stats3$vi), c(pam.stats2$corrected.rand, pam.stats2$vi), c(pam.stats3$corrected.rand, pam.stats3$vi), c(hc.stats2$corrected.rand, hc.stats2$vi), c(hc.stats3$corrected.rand, hc.stats3$vi))
colnames(external.validation) <- c("K-means 2", "K-means 3", "K-medoids 2", "K-medoids 3", "HC 2", "HCl 3")
rownames(external.validation) <- c("Rand", "VI")
