setwd("C:/Users/pprad/Desktop/Regression Analysis/Final Project")
library(car)

#------------------dataset-------------------------------------------------------------------------------
df <- read.csv('streamflow.csv')
df1 <- df[-1]
colnames(df1)<- c('Y','X1','X2', 'X3', 'X4', 'X5', 'X6', 'X7')


y <- df[2]
y

x <-df[-2]
X <- as.matrix(x)
inv.XX <- solve(t(X)%*%X)
H <- X%*%inv.XX%*%t(X)
lev <- diag(H)
#-------------------------------------------------------------------------------------------------------

#------------------Initial plots for summary-------------------------------------------------------------

boxplot(y, xlab='90th percentile of the annual maximum dtreamflow', main='Boxplot of Response variable Y')

boxplot(df$max90, xlab='90th percentile of the annual maximum dtreamflow', main='Boxplot of Response variable Y')

hist(df$max90, main ="Histogram for Y", xlab="90th percentile of the annual maximum dtreamflow")
hist(df$DRAIN_SQKM, main ="Histogram for X1", xlab="The drainage area")
hist(df$PPTAVG_BASIN, main ="Histogram for X2", xlab="The average basin precipitation")
hist(df$T_AVG_BASIN, main ="Histogram for X3", xlab="The average basin temperature")
hist(df$T_AVG_SITE, main ="Histogram for X4", xlab="The average temperature at the stream location")
hist(df$RH_BASIN, main ="Histogram for X5", xlab="The average relative humidity across the basin")
hist(df$MAR_PPT7100_CM, main ="Histogram for X6", xlab="The average March precipitation")
hist(df$RRMEDIAN, main ="Histogram for X7", xlab="The median relief ratio")

boxplot(df$max90, xlab='90th percentile of the annual maximum dtreamflow', main='Boxplot of Response variable Y')
boxplot(df$DRAIN_SQKM, main ="Boxplot for X1", xlab="The drainage area")
boxplot(df$PPTAVG_BASIN, main ="Boxplot for X2", xlab="The average basin precipitation")
boxplot(df$T_AVG_BASIN, main ="Boxplot for X3", xlab="The average basin temperature")
boxplot(df$T_AVG_SITE, main ="Boxplot for X4", xlab="The average temperature at the stream location")
boxplot(df$RH_BASIN, main ="Boxplot for X5", xlab="The average relative humidity across the basin")
boxplot(df$MAR_PPT7100_CM, main ="Boxplot for X6", xlab="The average March precipitation")
boxplot(df$RRMEDIAN, main ="Boxplot for X7", xlab="The median relief ratio")

pairs(df1)

library(car)
library(tidyselect)
linear_model <- lm(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, data =df1)
avPlots(linear_model)

cor(df1)

mean(df1$X1)
mean(df1$X2)
mean(df1$X3)
mean(df1$X4)
mean(df1$X5)
mean(df1$X6)
mean(df1$X7)

sd(df1$X1)
sd(df1$X2)
sd(df1$X3)
sd(df1$X4)
sd(df1$X5)
sd(df1$X6)
sd(df1$X7)

summary(linear_model)

###--------------------packages and libraries required-----------------------------------------
install.packages("leaps")
install.packages("HH")
install.packages("StepReg")
install.packages('/Library/gurobi901/mac64/R/gurobi_9.0-1_R_3.6.1.tgz', repos=NULL)
install.packages("slam", repos = "https://cloud.r-project.org")

#### Load HH, leaps, and StepReg packages ####

library(leaps)
library(HH)
library(StepReg)
library(olsrr)



#--------------Standardization and centralization------------------------------------------

k <- ols_step_both_p(linear_model,pent=0.10,prem=0.1,details=TRUE)
plot(k)

# Centralizing data
df2 <- data.frame(scale(x, center = TRUE, scale = FALSE))
df2
df1
pairs(df2)
#Scaling the data
df3 <- data.frame(scale(x))
df3 <-cbind(df1$Y, df3)
colnames(df3)<- c('Y', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7')
df3

pairs(df3)

df2 <-cbind(df1$Y, df2)
colnames(df2)<- c('Y', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7')
df2

df3 <-cbind(df1$Y, df3)
colnames(df3)<- c('Y', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7')
df3

pairs(df3)
#-----------------------------------------------------------------------------------------------
#--------------------------------VIF best model to avoid multicolinearity------------------------------------------------------------
vif(linear_model)

linear_model_x2 <- lm(Y ~ X1 + X3 + X4 + X5 + X6 + X7, data =df1)
vif(linear_model_x2)

linear_model_x3 <- lm(Y ~ X1 + X2 + X4 + X5 + X6 + X7, data =df1)
vif(linear_model_x3)

linear_model_x2_x3 <- lm(Y ~ X1 + X4 + X5 + X6 + X7, data =df1)
vif(linear_model_x2_x3)

linear_model_x4_x6 <- lm(Y ~ X1 + X2 +X3 + X5 + X7, data =df1)
vif(linear_model_x4_x6)

linear_model_x4_x6 <- lm(Y ~ X1 + X2 +X3 + X5 + X7, data =df3)
vif(linear_model_x4_x6)

linear_model_try <- lm(Y ~ X1 + X2 +X3 + X6, data =df1)
vif(linear_model_try)

summary(linear_model_try)

linear_model_2 <-lm(Y ~ X1 + X2+ X3 + X4 + X5 + X6 + X7, data =df2)
vif(linear_model_2)
cor(df2)

linear_model_3 <-lm(Y ~ X1 + X2+ X3 + X4 + X5 + X6 + X7, data =df3)
vif(linear_model_3)
cor(df3)


linear_model_temo <- lm(Y ~ X1 + X2 + X4 , data =df1)
vif(linear_model_temo)


#----------------------------------------------------------------------------------------------------

#-------------------------Removing Outlines------------------------------------------------------------
p <- ncol(X)
n <- nrow(df1)

outliers <- which(lev > 0.25*p/n)

leveragePlots(linear_model_x4_x6)
leveragePlots(linear_model_x2_x3)

hist(lev, main="Histogram of leverage points", xlab="leverage")
abline(v=0.25*p/n, lty=2, lwd=2, col="red")



p <- ncol(X)
n <- nrow(dataset_new)

outliers <- which(lev > 2*p/n)

outliers
dataset_new <- df1[-c(outliers),]



leveragePlots(linear_model_x4_x6)


linear_model_temo <- lm(Y ~ X1 + X4 + X5 + X6 +X7, data =dataset_new)
summary(linear_model_temo)

linear_model_temo <- lm(Y ~ X1 + X2 + X3 + X5 + X7 , data =dataset_new)
summary(linear_model_temo)



#-------------------------------------------------------------------------------------------------

#------------------------Stepwise regression------------------------------------------------------

linear_model_x4_x6 <- lm(Y~ X1+X2 +X3 +X5+X7 , data= dataset_new)
b <- ols_step_all_possible(linear_model_x4_x6)
plot(b)
b


#### Adjusted R2 ####

b.adjr = data.frame(n=b$n,predictors=b$predictors,adjr=b$adjr)
print(b.adjr)
print(b.adjr[c(1,6,16,26,31),])

bptest(linear_model_x4_x6)

#### Cp ####

b.cp = data.frame(n=b$n,predictors=b$predictors,cp=b$cp)
print(b.cp)
print(b.cp[c(1,6,16,26,31),])

# from results model with highest adjusted R-squared value is with 4 predictor variables, X1, X2, X3, X6 
#and adjusted R-squared value 0.278010647

#### AIC ####

b.aic = data.frame(n=b$n,predictors=b$predictors,aic=b$aic)
print(b.aic)
print(b.aic[c(1,6,16,26,31),])

# from results model with smallest AIC value is with 4 predictor variables, X1, X2, X3, X6 
#and AIC value 6175.200, supporting earlier selected model


#### BIC ####

b.bic = data.frame(n=b$n,predictors=b$predictors,bic=b$sbic)
print(b.bic)
print(b.bic[c(1,6,16,26,31),])

# from results model with smallest BIC value is with 4 predictor variables, X1, X2, X3, X6 
#and BIC value 5343.924, supporting earlier selected model

# from results model with smallest BIC value is with 4 predictor variables, X2, X3, X4, X7
#and BIC value 5403.399, supporting earlier selected model


#### PRESS ####

b.press = data.frame(n=b$n,predictors=b$predictors,press=b$msep)
print(b.press)
print(b.press[c(1,6,16,26,31),])


# from results model with smallest Press value is with 4 predictor variables, X1, X2, X3, X6 
#and Press value 23833448900, supporting earlier selected model

# from results model with smallest press value is with 4 predictor variables, X2, X3, X4, X7
#and press value 29396553667, supporting earlier selected model









install.packages("EnvStats")
library(EnvStats)

lmfit <- lm(Y~X1+X2+X3+X7,data=dataset_new)
summary(lmfit)


#----------------test Homoscedasticity------------------------

bp ( )
######################## Box-Cox Transformations ###############################

boxcox(lmfit,optimize = TRUE)
boxcox.summary <- boxcox(lmfit, optimize = FALSE)

plot(boxcox.summary$objective ~ boxcox.summary$lambda, xlab=bquote(lambda), ylab="Likelihood",
     main="Likelihood vs. Power (lambda)", type="l",lwd=2)
abline(v=boxcox(lmfit,optimize = TRUE)$lambda,lty=2,col="red")

boxcox.data <- cbind(boxcox.summary$lambda, boxcox.summary$objective)
colnames(boxcox.data) <- c("lambda","Likelihood")

boxcox.summary <- boxcox(lmfit, optimize = TRUE)
lambda <- boxcox.summary$lambda

trans.Y <- dataset_new$Y^(lambda)
dataset_new <- cbind(dataset_new, trans.Y)
colnames(dataset_new)[9] <- c("trans.Y")

############################# Transformation #####################################

boxcox.lmfit <- lm(trans.Y ~ X1+X2+X3+X7, data=dataset_new)


res.boxcox.lmfit <- resid(boxcox.lmfit)

bptest(boxcox.lmfit)



############## Lack of Fit Test (Linearity) ##########################

full.lmfit <- lm(trans.Y ~ factor(X1)+, data=plasma)
anova(full.lmfit)

reduced.lmfit <- lm(trans.Y ~ X, data=plasma)
anova(reduced.lmfit)

anova(reduced.lmfit, full.lmfit)


################ Durbin-Watson Test (Independence) ##########################

library(lmtest)
dwtest(boxcox.lmfit, alternative = "two.sided")

################ Breusch-Pagan Test (Constant Variance) ##########################

bptest(boxcox.lmfit)


####################### QQ Plot & Shapiro-Wilk Test (Normality) ##################

qqnorm(res.boxcox.lmfit);qqline(res.boxcox.lmfit)
shapiro.test(res.boxcox.lmfit)

####################### Outliers? ########################

std.res <- rstudent(boxcox.lmfit)
plot(std.res ~ fitted(boxcox.lmfit),ylab="residual", xlab="fitted", ylim=c(-2,2),
     main="residual plot against fitted values")
abline(h=0)



#------------------------additional transformations and testing----------------------------------------------


library(sandwich)
coeftest(lmfit, vcov. = vcovHC(lmfit, "HC1"))


residd <-lmfit$residuals
varfunc <- lm(log(residd^2)~ log(X1+X2+X3+X7), data= dataset_new)
summary((varfunc))
bptest(varfunc)

nextst <- exp(varfunc$fitted.values)
gls <- lm(Y ~ X1+X2+X3+X7, weights= 1/sqrt(nextst), data= dataset_new)

summary(gls)

bptest(gls)

gls_resid <- rstudent(gls)
plot(gls_resid, xlab='time', ylab='residuals', main='plot of residual against index')
abline(h=0, lwd=2)


plot(gls_resid~ fitted(gls), xlab='fitted values', ylab='residuals', main='plot of residual against fitted values')
abline(h=0, lwd=2)

dwtest(gls, alternative = "two.sided")


dataset_new
142+125

vif(lmfit)
s

ols_plot_dffits(boxcox.lmfit)

###################### 2. Cook's D ##########################################

ols_plot_cooksd_chart(boxcox.lmfit)

###################### 3. DFBETAS ##########################################

ols_plot_dfbetas(boxcox.lmfit)

summary(boxcox.lmfit)

library(MASS)

bc <- boxcox(lmfit)
lambda <- bc$x[which(bc$y==max(bc$y))]
str(bc)
bc.power <-bc$X[which.max(bc$y)]
bc.power
