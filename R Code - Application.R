### ======================================================================= ###
###                             APPLICATION                                 ###  
### ----------------------------------------------------------------------- ###
### Estimation and diagnostic Semiparametric additive beta regression model ###
### with varying precision parameter                                        ###
###                                                                         ###
###  Authors: G. Ibacache; Caamaño-Carrillo, C and Conde, F.                ###
### ======================================================================= ###

rm(list = ls())
###############################################################################
###                Caschool (Test score data in California)                 ###
###############################################################################
library(pracma)
library(MASS)
library("betareg")
library(gamlss)
library(caret)

setwd("C:/Users/chcaa/OneDrive/Escritorio/Papers/Beta Semiparametric/Aplicación")
source("Model_semip_beta_fixed.R")
source("Model_semip_beta_varying.R")
source("GCV.R")

Caschool <- read.csv("dataset.csv",header=T,sep = ",",dec = "." )
head(Caschool)
summary(Caschool)
attach(Caschool)

Y <- (mathscr-min(mathscr))/(max(mathscr)-min(mathscr))
summary(Y)
data <- cbind(Caschool,Y)
data1<-subset(data,Y>0 & Y<1)
Y<-data1$Y

###  Figure 1 of paper

par(mfrow=c(1,2))
hist(Y,main = "Histogram",xlab="y",col=4,freq = FALSE)
boxplot(Y,main = "Boxplot ",boxfill = 4)

readscr<-data1$readscr
avginc<-data1$avginc
mealpct<-data1$mealpct
calwpct<-data1$calwpct

###  Figure 2 of paper
par(mfrow=c(2,2))
plot(readscr,Y,ylab = "y",main = "a)")  
plot(avginc,Y,ylab = "y",main = "b)")   
plot(mealpct,Y,ylab = "y",main = "c)")  
plot(calwpct,Y,ylab = "y",main = "d)")

###############################################################################
###             Fitting statistical models to the dataset                   ###
###############################################################################

### MODEL 1: application paper 

model_1 <- betareg(Y ~ readscr + avginc + mealpct + calwpct)
summary(model_1)
AIC(model_1,k=log(418))
### cross-validation 
data1 <- data.frame(Y, readscr, avginc, mealpct, calwpct)
GCV_BR(Y ~ readscr + avginc + mealpct + calwpct,p=6, data1,k=20)


### MODEL 2: application paper 

model_2 <- betareg(Y ~ readscr + avginc | mealpct + calwpct)
summary(model_2)
AIC(model_2,k=log(418))
GCV_BR(Y ~ readscr + avginc | mealpct + calwpct,p=6, data1,k=20)

### MODEL 3: application paper 

model_3<-betareg(Y ~ readscr + I(avginc^2) | mealpct + I(calwpct^2))
summary(model_3)
AIC(model_3,k=log(418))
GCV_BR(Y ~ readscr + I(avginc^2) | mealpct + I(calwpct^2),p=6,data1,k=20)


### MODEL 4: application paper

model_4<- betareg(Y ~ readscr + avginc + I(avginc^2) | mealpct + calwpct + I(calwpct^2))
summary(model_4)
AIC(model_4,k=log(418))
GCV_BR(Y ~ readscr + avginc + I(avginc^2) | mealpct + calwpct + I(calwpct^2),p=8,data1,k=20)


### MODEL 5: application paper

Y <- Y
X <- readscr
T1 <- avginc

fit.fixed <- model.fixed(X=readscr,T1=avginc,Y=Y,lambda1=0.002)


fit.fixed$beta
fit.fixed$SD_beta
fit.fixed$phi
fit.fixed$SD_phi
fit.fixed$GCV
fit.fixed$Like
fit.fixed$AIC
fit.fixed$R2_adj



### MODEL 6: application paper 

model_6 <- gamlss(formula = Y ~ readscr+cs(avginc),
                  sigma.formula = ~mealpct+cs(calwpct),
                  family=BEo(mu.link = "log", sigma.link = "log"))
summary(model_6)
summary(model_6)$terms
logLik(model_6)
AIC(model_6,k=log(418))
R2gamlss(model_6)
gcvgamlss<-gcv_gamlss(data1,k = 20, p = 4, seed = 123)


###############################################################################
### PROPOSED MODEL: Semiparametric additive beta regression model           ###
###                 with varying precision parameter                        ###  
###############################################################################

Y <- Y
x1 <- readscr
z1 <- mealpct
T1 <- avginc
t1 <- calwpct

fit.varying <- model.varying(x1=readscr,z1=mealpct,T1=avginc,t1=calwpct,
                             delta1=2100,lambda1=500,tol=10^-5)

fit.varying$beta
fit.varying$alpha
fit.varying$SD_BETA_n
fit.varying$SD_Alpha_n
fit.varying$GCV
fit.varying$Like
fit.varying$AIC
fit.varying$R2_adj


###  Figure 3 of paper  

par(mfrow=c(1,2))

plot(fit.varying$T1_0,fit.varying$F1_en,type = "l",ylim = c(0.3,1.1),xlab = "avginc",ylab = "f1(avginc)",
     main="a)")
lines(fit.varying$T1_0,fit.varying$F1_en+2*fit.varying$SD_F1_n,type = "l", lty=2 ,col="red")
lines(fit.varying$T1_0,fit.varying$F1_en-2*fit.varying$SD_F1_n,type = "l", lty=2 ,col="red")

plot(fit.varying$t1_0,fit.varying$F1star_en ,type = "l",ylim = c(-3.5,1.5),xlab ="calwpct" ,ylab ="f1*(calwpct)",
     main="b)")
lines(fit.varying$t1_0,fit.varying$F1star_en+2*fit.varying$SD_F1star_n,type = "l",lty=2,col="red")
lines(fit.varying$t1_0,fit.varying$F1star_en-2*fit.varying$SD_F1star_n,type = "l",lty=2,col="red")

###  Figure 4 of paper

plot(fit.varying$residuals,ylab = "residuals",xlim = c(0,425),ylim = c(-4,5.5),cex = 0.7,bg="red",pch=21)
abline(h=2.7*sd(fit.varying$residuals),lwd=1,col="black",lty=2)
abline(h=-2.7*sd(fit.varying$residuals),lwd=1,col="black",lty=2)
identify(fit.varying$residuals, cex=1,lwd=2, n = 4)

###  Figure 5 of paper: Partial Residuals

par(mfrow=c(1,2))
plot(avginc,fit.varying$par_res1,ylim = c(-4,5),xlab="avginc",
     ylab = expression(r[p1]),main="a)")
lines(fit.varying$T1_0,fit.varying$F1_en,type = "l", lty=2,lwd=2 ,col="red")

plot(calwpct,fit.varying$par_res2,ylim = c(-4,4),xlab="calwpct",
     ylab =expression(r[p2]),main ="b)")
lines(fit.varying$t1_0,fit.varying$F1star_en,type = "l",lty=2,lwd=2 ,col="red")

###  Figure 6 of paper: Leverage

court <- 3*mean(fit.varying$leverage)
aa=1:length(Y)

plot(aa,fit.varying$leverage,type='b',xlim = c(0,425),ylim=c(0,max(1.1*court,max(fit.varying$leverage)))
     ,xlab='Index',ylab='GL',font.lab=2,cex = 0.7,bg="magenta",pch=21)
lines(aa,court*rep(1,length(Y)),lty=2,lwd=2,col='black')
identify(aa,fit.varying$leverage,n=1)

###############################################################################
### PROPOSED MODEL: Local influence analysis                                ###
###############################################################################


### Figure 7 of paper: Case 1: Case weighting

### Parametric component (mean)
par(mfrow=c(2,2))
plot(1:length(abs(fit.varying$wp1)), abs(fit.varying$wp1), ylim = c(0,1),cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="a)")
abline(h=1.2*fit.varying$lwp1,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$wp1)), abs(fit.varying$wp1), n = 1)

### Parametric component (precision)
plot(1:length(abs(fit.varying$wp2)), abs(fit.varying$wp2), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="b)")
abline(h=0.8*fit.varying$lwp2,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$wp2)), abs(fit.varying$wp2), n = 1 )

### Nonparametric component (mean)
plot(1:length(abs(fit.varying$wp3)), abs(fit.varying$wp3), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="c)")
abline(h=1.2*fit.varying$lwp3,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$wp3)), abs(fit.varying$wp3), n = 1)

### Nonparametric component (precision)
plot(1:length(abs(fit.varying$wp4)), abs(fit.varying$wp4), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="d)")
abline(h=1.2*fit.varying$lwp4,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$wp4)), abs(fit.varying$wp4), n = 1)

### Figure 8 of paper: Case 2: Response perturbation

### Parametric component (mean)
par(mfrow=c(2,2))
plot(1:length(abs(fit.varying$rp1)),abs(fit.varying$rp1),ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="a)")
abline(h=14*fit.varying$lrp1,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$rp1)), abs(fit.varying$rp1), n = 2)

### Parametric component (precision)
plot(1:length(abs(fit.varying$rp2)), abs(fit.varying$rp2), ylim = c(0,1),cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="b)")
abline(h=16*fit.varying$lrp2,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$rp2)), abs(fit.varying$rp2), n = 2)

### Nonparametric component (mean)
plot(1:length(abs(fit.varying$rp3)), abs(fit.varying$rp3), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="c)")
abline(h=15*fit.varying$lrp3,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$rp3)), abs(fit.varying$rp3), n = 1)

### Nonparametric component (precision)
plot(1:length(abs(fit.varying$rp4)), abs(fit.varying$rp4), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="d)")
abline(h=24*fit.varying$lrp4,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$rp4)), abs(fit.varying$rp4), n = 1)


### Figure 9 of paper: Case 3: Perturbation of predictor X

### Parametric component (mean)
par(mfrow=c(2,2))
plot(1:length(abs(fit.varying$pp1)), abs(fit.varying$pp1), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="a)")
abline(h=0.4*fit.varying$lpp1,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$pp1)), abs(fit.varying$pp1), n = 0)

### Parametric component (precision)

plot(1:length(abs(fit.varying$pp2)), abs(fit.varying$pp2), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="b)")
abline(h=0.24*fit.varying$lpp2,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$pp2)), abs(fit.varying$pp2), n = 0)

### Nonparametric component (mean)
plot(1:length(abs(fit.varying$pp3)), abs(fit.varying$pp3), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="c)")
abline(h=0.55*fit.varying$lpp3,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$pp3)), abs(fit.varying$pp3), n = 0)

### Nonparametric component (precision)
plot(1:length(abs(fit.varying$pp4)), abs(fit.varying$pp4), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="d)")
abline(h=0.3*fit.varying$lpp4,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$pp4)), abs(fit.varying$pp4), n = 1)


### Figure 10 of paper: Case 4:  Perturbation of predictor Z

### Parametric component (mean) 
par(mfrow=c(2,2))
plot(1:length(fit.varying$zp1), abs(fit.varying$zp1), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="a)")
abline(h=fit.varying$lzp1/150,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$zp1)), abs(fit.varying$zp1), n = 1)

### Parametric component  (precision)
plot(1:length(abs(fit.varying$zp2)), abs(fit.varying$zp2), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="b)")
abline(h=fit.varying$lzp2/22,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$zp2)), abs(fit.varying$zp2), n = 1)

### Nonparametric component (mean) 
plot(1:length(abs(fit.varying$zp3)), abs(fit.varying$zp3), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="c)")
abline(h=fit.varying$lzp3/150,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$zp3)), abs(fit.varying$zp3), n = 1)

### Nonparametric component (precision)
plot(1:length(abs(fit.varying$zp4)), abs(fit.varying$zp4), ylim = c(0,1), cex = 0.7,
     bg="red",pch=21,xlab = "Index",ylab =expression(B[i]),main ="d)")
abline(h=fit.varying$lzp4/130,lwd=1,col="black",lty=2)
identify(1:length(abs(fit.varying$zp4)), abs(fit.varying$zp4), n = 1)













