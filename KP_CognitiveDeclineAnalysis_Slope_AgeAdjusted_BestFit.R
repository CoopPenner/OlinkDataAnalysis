# Load settings
options(stringsAsFactors = FALSE)

# Load packages
library(data.table)
library(dendextend)
library(stabs)
library(lars)
library(GGally)
library(mlbench) # for cross validation
library(caret)
library(ggplot2)
library(plyr)
library(dplyr) # data manipulation
library(ROCR) #for plotting ROC
library(e1071)
library(ggpubr) # for plotting
library(tidyverse) #for data manipulation
library(ROCR) #for ROC
library(pROC) #for generating ROC curves
library(lme4) # for LME model
library(survival) # for survival analysis
library(survminer) #for survival analysis
library(ggfortify) # for plotting 
library(gridExtra)
library(grid)
library(sjPlot)
library(sjmisc)
library(nlme)
library(lme4)
library(knitr)
library(kableExtra)
library(lubridate)
library(glmnet)
library(leaps)
options(warn = -1)

# Directory
getwd()
setwd("/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/")

# Load Data
#olink.dat <- read.csv("Olink data 20210120.csv") #original olink data set
olink.dat <- read.csv("20220517_olinkqualitycontrolvalues.csv")
original_olink.dat <- read.csv("20220517_olinkqualitycontrolvalues.csv") #olink data set with quality control
olink_log <- as.data.frame(sapply(original_olink.dat,as.numeric))
olink_log <- log(olink_log[12:length(olink_log)])
olink.dat <- cbind(original_olink.dat[1:11], olink_log)

olink.dat <- olink.dat[-which(olink.dat$INDDID == 111864),] #remove patient demented at baseline 
olink_pd.dat <- olink.dat[olink.dat$Group == "PD",]
olink_pd.dat <- olink_pd.dat[colSums(!is.na(olink_pd.dat)) > 2]
olink_pd.dat <- olink_pd.dat[!is.na(olink_pd.dat$AGE..YRS),] 

#################################################################################
# (1) DEFINE DRS SLOPE FROM BEST-FITTING LINE
# remove DRS taken before sample,
# compute DRS slope from best-fitting line,
# assign fast vs. slow 
#################################################################################

# Extract DRS data for PD patients 
drs.dat <- read.csv("01_pd_drs.csv") #PD DRS data
pd_drs.dat <- drs.dat %>% filter(drs.dat$INDDID %in% olink_pd.dat$INDDID)

pdnum <- nrow(olink_pd.dat) # number of PD patients (50)
data1 <- data.frame(INDDID = numeric(pdnum),
                    DRS.initial = numeric(pdnum),
                    date.initial = character(pdnum),
                    DRS.final = numeric(pdnum),
                    date.final = character(pdnum),
                    yearDiff = numeric(pdnum),
                    DRSnum = numeric(pdnum),
                    declineRate = numeric(pdnum),
                    declineSpeed = character(pdnum),
                    fastProgression = numeric(pdnum))
path_out2 <- "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/figures/"

proteinDrawDiff <- list()
DRSrowstoremove <- c()

#for-loop to remove DRS taken before sample, compute DRS slope, and assign fast vs. slow
for (i in 1:pdnum) {
  data1$INDDID[i] <- olink_pd.dat$INDDID[i]
  
  initialIndex <- match(olink_pd.dat$INDDID[i], pd_drs.dat$INDDID)

  allIndex <- which(pd_drs.dat$INDDID %in% olink_pd.dat$INDDID[i])
  allDRS <- pd_drs.dat$DRSTotalAge[allIndex]
  allYears <- (as.Date(pd_drs.dat$TestDate[allIndex], "%m/%d/%Y"))
  allYears <- decimal_date(allYears)
  
  if(sum(is.na(allDRS)) > 0) {
    allIndex <- allIndex[-length(allIndex)]
    allDRS <- allDRS[-length(allDRS)]
    allYears <- allYears[-length(allYears)]
  }
  
  #when olink sample was taken 
  sampleDate <- olink_pd.dat$SampleDate[i]
  sampleDate <- decimal_date(as.Date(sampleDate, "%m/%d/%Y"))
  
  #remove DRS taken before sample
  dateDiff <- allYears - sampleDate
  toRemoveIndex <- which(dateDiff < 0)
  if(length(toRemoveIndex) != 0) {
    allYears <- allYears[-toRemoveIndex]
    allDRS <- allDRS[-toRemoveIndex]
    DRSrowstoremove <- c(DRSrowstoremove, allIndex[toRemoveIndex])
  }
  data1$DRSnum[i] <- length(allDRS[!is.na(allDRS)])
  
  if(length(allDRS)<1){
    data1$declineRate[i] <- NaN 
    next 
  }

  data1$DRS.initial[i] <- allDRS[1]
  data1$date.initial[i] <- format(date_decimal(allYears[1]),"%m/%d/%Y")
  data1$DRS.final[i] <- allDRS[length(allDRS)]
  data1$date.final[i] <- format(date_decimal(allYears[length(allYears)]),"%m/%d/%Y")
  data1$yearDiff[i] <- allYears[length(allYears)]-allYears[1]
  
  proteinDrawDiff[i] <- min(abs(dateDiff))
  
  # find best fitting line
  bestLine <- lm(allDRS ~ allYears) 
  
  #set decline rate to slope of best fitting line
  slope <- coef(bestLine)[2]
  data1$declineRate[i] <- coef(bestLine)[2] 
  data1$bestFitIntercept[i] <- summary(bestLine)$coefficients[1,1]
  
  
  jpeg(file=paste(path_out2, data1$INDDID[i],"bestFitLineDRS.jpeg"))
  plot(allYears, allDRS,
       main = paste("INDDID", data1$INDDID[i], ":", slope),
       ylab = "DRS", xlab = "Year", cex.lab = 1.5,
       pch = 21, bg = "black")
  abline(bestLine, col = "blue", lwd = 3, lty = 2)
  dev.off()  
  
  #assign fast vs. slow
  if (data1$declineRate[i] < -0.5) {
    data1$declineSpeed[i] <- "Fast"
    data1$fastProgression[i] <- 1
  }
  else {
    data1$declineSpeed[i] <- "Slow"
    data1$fastProgression[i] <- 0
  }
  print(i)
}

#remove DRS scores before sample
pd_drs.dat <- pd_drs.dat[-DRSrowstoremove,]

#look at when DRS was taken relative to protein draw
hist(as.numeric(proteinDrawDiff),
     main = "Absolute Difference Between DRS and Protein Draw (Years) ")
summary(as.numeric(proteinDrawDiff))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.07178 0.09016 0.45205 
       
data1 <- data1 %>% arrange(declineRate)
slowdata <- subset(data1, declineSpeed == "Slow")
nrow(slowdata) #24 slow progressing patients 
fastdata <- subset(data1, declineSpeed == "Fast")
nrow(fastdata) #25 fast progressing patients 

#################################################################################
# (2) CHARACTERIZATION OF COGNITIVE PROGRESSION GROUPS
# LME model to look at fast vs. slow groups
#################################################################################

data2 <- merge(pd_drs.dat, data1, by = "INDDID")

# Calculate # of years since first visit for each test
drsnum <- nrow(data2)
for (i in 1:drsnum){
  data2$yearsToVisit[i] <- time_length(difftime(as.Date(data2$TestDate[i], "%m/%d/%Y"), 
                                                as.Date(data2$date.initial[i], "%m/%d/%Y")),
                                       "years")
}

# Add olink proteins
data2 <- merge(data2, olink_pd.dat, by = "INDDID")

# Code sex from female male into 0 (Female) and 1 (Male)
data2$Sex <- factor(data2$Sex, levels=c("Female","Male"), labels=c(0,1))
data2$fastProgression <- factor(data2$fastProgression, levels = c(0, 1), labels=c("Slow","Fast"))

# Train LME model to see how fast and slow progressors differ
DRS <- lmer(DRSTotalAge ~ fastProgression*yearsToVisit + AgeatTest + Sex +
              (1 | INDDID), data=data2, REML=F) #no disease duration 

data(efc)
theme_set(theme_sjplot())
p1 = plot_model(DRS, type = "pred", terms = c("yearsToVisit", "fastProgression"),
                title = "Cognitive Change",
                axis.title = c("Years Since Plasma Sample","DRS"),
                axis.labels = NULL, legend.title = "Progression Group",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set2", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.25,0.18), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 20),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 20, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
p1

  
#################################################################################
# (3) MODEL SELECTION
# identify slope-associated proteins,
# stability selection,
# exhaustive search for model selection,
# create model and predict DRS slope
#################################################################################

#IDENTIFY SLOPE-ASSOCIATED PROTEINS
protein_num <- ncol(olink_pd.dat) - 11 #olink panel size #1214

data3 <- merge(data1, olink_pd.dat, by = "INDDID")
data3$Sex <- factor(data3$Sex, levels=c("Female","Male"), labels=c(0,1))

slopemodel.df <- data.frame(Protein = character(protein_num),
                            Protein.coef = numeric(protein_num),
                            Protein.pval = numeric(protein_num),
                            Protein.FDR = numeric(protein_num),
                            Male.coef = numeric(protein_num),
                            Male.pval = numeric(protein_num),
                            Male.FDR = numeric(protein_num),
                            Age.coef = numeric(protein_num),
                            Age.pval = numeric(protein_num),
                            Age.FDR = numeric(protein_num))

# LINEAR REGRESSION 
# DRS slope as a function of protein + sex + age 
for (i in 1:protein_num) {
  protein_name <- names(olink_pd.dat)[i+11]
  slopemodel.df$Protein[i] <- protein_name
  formula_i <- paste("-declineRate ~", protein_name, "+Sex", "+AGE..YRS", sep = "")
  model_i <- lm(formula_i, data3)
  slopemodel.df$Protein.coef[i] <- summary(model_i)$coefficients[2,1]
  slopemodel.df$Protein.pval[i] <- summary(model_i)$coefficients[2,4]
  slopemodel.df$Male.coef[i] <- summary(model_i)$coefficients[3,1]
  slopemodel.df$Male.pval[i] <- summary(model_i)$coefficients[3,4]
  slopemodel.df$Age.coef[i] <- summary(model_i)$coefficients[4,1]
  slopemodel.df$Age.pval[i] <- summary(model_i)$coefficients[4,4]
  print(i)
}

# Benjamini-Hochberg Procedure (p-value adjustment)
slopemodel.df$Protein.FDR <- p.adjust(slopemodel.df$Protein.pval, method = 'BH')
slopemodel.df$Male.FDR <- p.adjust(slopemodel.df$Male.pval, method = 'BH')
slopemodel.df$Age.FDR <- p.adjust(slopemodel.df$Age.pval, method = 'BH')

# Summary
summary(slopemodel.df$Protein.FDR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4063  0.7283  0.8612  0.8111  0.9108  0.9994 
summary(slopemodel.df$Male.FDR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2249  0.2249  0.2249  0.3484  0.5065  0.9835 
summary(slopemodel.df$Age.FDR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9993  0.9993  0.9993  0.9993  0.9993  0.9993 

summary(slopemodel.df$Protein.pval)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0003341 0.1803452 0.4249561 0.4434535 0.6828746 0.99937287 
summary(slopemodel.df$Male.pval)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01176 0.09002 0.11350 0.23311 0.38000 0.98348 
summary(slopemodel.df$Age.pval)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.006634 0.571736 0.823004 0.710892 0.887284 0.999291 

sum(slopemodel.df$Protein.pval < 0.1) #211 #202
sum(slopemodel.df$Protein.pval < 0.05) #110 slope-associated proteins #95
sum(slopemodel.df$Protein.pval < 0.01) #15 #14

slopemodel.df <- slopemodel.df %>% arrange(Protein.pval)
path_out = "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/"
fileName = paste(path_out, '20220726_slope_linearregression.csv', sep = '')
write.csv(slopemodel.df,
          file = fileName,
          row.names=FALSE)

#STABILITY SELECTION: elastic net regression with 10-fold CV (100 proteins)
sig.proteins <- slopemodel.df$Protein[slopemodel.df$Protein.pval < 0.05]
idx.sig <- match(sig.proteins, names(data3))
Y <- data3$declineRate
X <- model.matrix.lm(declineRate~ ., data=data3[,c(8,17,19,idx.sig)], na.action="na.pass")[,-1]#proteins,sex,age

#First, using 100% samples
count.matrix <- data.frame(Param=colnames(X), Non.Zero=numeric(length(colnames(X))))
for (i in 1:1000) {
  set.seed(i)
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),] #get non-zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix$Param)
  count.matrix$Non.Zero[idx_i] <- count.matrix$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (100% of Samples)",
       x = "Times with Non-Zero Coefficient", y = "Parameters")
sum(count.matrix$Non.Zero >= 1000) #26 #22
sum(count.matrix$Non.Zero >= 900)  #26 #25

# Second, using 90% samples:
count.matrix_0.9 <- data.frame(Param = colnames(X), Non.Zero = numeric(length(colnames(X))))
for (i in 1:1000) {
  set.seed(i)
  rand_0.9 <- sample(1:49, 44, replace = FALSE) #44 total samples (90% of 49 is 44)
  Y <- data3$declineRate[rand_0.9]
  X <- model.matrix.lm(declineRate ~ ., data = data3[rand_0.9, c(8,17,19,idx.sig)], na.action = "na.pass")[, -1] #sex, proteins
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix_0.9$Param)
  count.matrix_0.9$Non.Zero[idx_i] <- count.matrix_0.9$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix_0.9) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (90% samples, DRS Slope-Associated Proteins + Sex + Age)",
       x = "Times with non-zero coefficient", y = "Parameters")
sum(count.matrix_0.9$Non.Zero == 1000) #0 #1
sum(count.matrix_0.9$Non.Zero >= 900)  #11 #14
sum(count.matrix_0.9$Non.Zero >= 500)  #25 #24

# Lastly, using 70% samples:
count.matrix_0.7 <- data.frame(Param = colnames(X), Non.Zero = numeric(length(colnames(X))))
for (i in 1:1000) {
  set.seed(i)
  rand_0.7 <- sample(1:49, 34, replace = FALSE) #70% of 49 is 34
  Y <- data3$declineRate[rand_0.7]
  X <- model.matrix.lm(declineRate ~ ., data = data3[rand_0.7, c(8,17,19,idx.sig)], na.action = "na.pass")[, -1] #sex, proteins
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix_0.7$Param)
  count.matrix_0.7$Non.Zero[idx_i] <- count.matrix_0.7$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix_0.7) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (70% of Samples)",
       x = "Times with Non-Zero Coefficient", y = "Parameters")

sum(count.matrix_0.7$Non.Zero == 1000) #0
sum(count.matrix_0.7$Non.Zero >= 900)  #2 #4
sum(count.matrix_0.7$Non.Zero >= 500)  #20 #20
sum(count.matrix_0.7$Param[count.matrix_0.7$Non.Zero >= 500] %in%
      count.matrix$Param[count.matrix$Non.Zero == 1000])
# 19/20 70% >= 500 in 100% == 1000
count.matrix$Param[count.matrix$Non.Zero == 1000][!(count.matrix$Param[
  count.matrix$Non.Zero == 1000] %in% 
    count.matrix_0.7$Param[count.matrix_0.7$Non.Zero > 500])]
#"LAT2"  "LTA4H" "MVK"   "EPO"   "NUB1"  "NELL1" "SDC4"   
#"ERP44" "PTN.1" "SDC4" 

fileName = paste(path_out, 
                 "20220526 stability selection (100-pct data).csv",
                 sep = '')
write.csv(count.matrix %>% arrange(1000-Non.Zero), 
          file = fileName,
          row.names = FALSE)
fileName = paste(path_out, 
                 "20220526 stability selection (90-pct data).csv",
                 sep = '')
write.csv(count.matrix_0.9 %>% arrange(1000-Non.Zero), 
          file = fileName,
          row.names = FALSE)
fileName = paste(path_out, 
                 "20220526 stability selection (70-pct data).csv",
                 sep = '')
write.csv(count.matrix_0.7 %>% arrange(1000-Non.Zero), 
          file = fileName,
          row.names = FALSE)

# MODEL SELECTION WITH EXHAUSTIVE SEARCH

top <- sum(count.matrix$Non.Zero >= 1000) 
# Now, let's perform exhaustive search for the top proteins (non-zero 1000)
ss_19 <- count.matrix$Param[count.matrix$Non.Zero == 1000] #top 22 proteins
idx_19 <- match(ss_19, names(data3))
olink_19 <- data3[,c(8,17,19,idx_19)]

#Exhaustive Search
fit.exh <- regsubsets(declineRate~., olink_19, nvmax = top,
                      method = 'exhaustive', really.big = TRUE)
summary.exh <- summary(fit.exh)

# Graph Rsq, Cp, BIC using training data
par(mfrow=c(3,1), mar=c(2.5,4,0.5,1), mgp=c(1.5,0.5,0))
plot(summary.exh$rsq, xlab="Number of predictors",
     ylab="R-squared", col="green", type="p", pch=16)
abline(v=12,lty = 2, lwd = 1, col = "black")
plot(summary.exh$cp, xlab="Number of predictors",
     ylab="Cp", col="red", type="p", pch=16)
abline(v=12,lty = 2, lwd = 1, col = "black")
plot(summary.exh$bic, xlab="Number of predictors",
     ylab="BIC", col="blue", type="p", pch=16)
abline(v=12,lty = 2, lwd = 1, col = "black")

# Results of exh selection
all.models <- data.frame(Var.Num = 1:top,
                         Formula = character(top),
                         R.sqrd = numeric(top),
                         Cp = numeric(top),
                         BIC = numeric(top),
                         AIC = numeric(top),
                         Full.Cor = numeric(top))
for (i in 1:top){
  #Save the essentials
  all.models$R.sqrd[i] <- summary.exh$rsq[i]
  all.models$Cp[i] <- summary.exh$cp[i]
  all.models$BIC[i] <- summary.exh$bic[i]
  #Get the names
  names_i <- names(olink_19[summary.exh$which[i,]])[-1]
  formula_i <- paste("declineRate~", paste(names_i, collapse = "+"), sep = "")
  all.models$Formula[i] <- formula_i
  #Fit model
  model_i <- lm(formula_i, data3)
  #AIC
  all.models$AIC[i] <- AIC(model_i)
  #Get correlation between actual and predicted
  all.models$Full.Cor[i] <- cor(data3$declineRate, model_i$fitted.values,
                                method = 'pearson')
}
all.models.exh <- all.models
View(all.models.exh)

# Let's generate the following loop

# Create empty model list
# Create empty data frames from RMSE, R.sqr, MAE
# For each model in fit.exh
# Perform 10-fold CV repeated 100x and save in list
# Save resample data to data.tables

model.list <- list()
RMSE.exh <- matrix(0, nrow = 1000, ncol = top)
R.sqr.exh <- matrix(0, nrow = 1000, ncol = top)
MAE.exh <- matrix(0, nrow = 1000, ncol = top)
for (i in 1:top){
  model_i <- train(formula(all.models.exh$Formula[i]), 
                   data = data3, method = "lm",
                   trControl = trainControl(method = "repeatedcv",
                                            number = 10,
                                            repeats = 100))
  model.list[[i]] <- model_i
  RMSE.exh[,i] <- model.list[[i]]$resample$RMSE
  R.sqr.exh[,i] <- model.list[[i]]$resample$Rsquared
  MAE.exh[,i] <- model.list[[i]]$resample$MAE
}

which.min(apply(RMSE.exh, 2, mean))  # 10 #11
which.max(apply(R.sqr.exh, 2, mean)) # 13 #11
which.min(apply(MAE.exh, 2, mean))   # 10 #11

boxplot(RMSE.exh, main = 'RMSE for 100x 10-fold CV', 
        ylab = 'RMSE', xlab = 'Parameter #')
boxplot(R.sqr.exh, main = 'Rsquared for 100x 10-fold CV', 
        ylab = 'Rsqrd', xlab = 'Parameter #')
boxplot(MAE.exh, main = 'MAE for 100x 10-fold CV',
        ylab = 'MAE', xlab = 'Parameter #')

# Save results
all.models.exh <- as.data.frame(all.models.exh)
fileName = paste(path_out, '20220526 all models exhaustive search.csv',sep = '')
write.csv(all.models.exh, 
          file = fileName,
          row.names = FALSE)

# My repeatCV
my.repeat.cv <- function(formula, data, kfold, repeats) {
  # Generate output data frame
  output.df <- data.frame(Iter = character(kfold*repeats), 
                          Rsquared.train = numeric(kfold*repeats),
                          RMSE.train = numeric(kfold*repeats),
                          Train.cor = numeric(kfold*repeats),
                          Test.cor = numeric(kfold*repeats))
  # Y variable
  yvar <- strsplit(formula, split = '~')[[1]][1]
  # Generate kfold.vector
  kfold.vector <- numeric(nrow(data))
  # Calculate AVERAGE size per fold (will vary if kfold is not a factor of nrow(data))
  size.per.fold <- nrow(data)/kfold
  # For k folds 1:kfold
  for (k in 1:kfold) {
    # Determine size of fold
    kfold.size <- length((round(size.per.fold*(k-1))+1):round(size.per.fold*k))
    # Replace vector entries with values indicating fold
    kfold.vector[(round(size.per.fold*(k-1))+1):round(size.per.fold*k)] <- rep(k, kfold.size)
    # Ex: kfold = 10, nrow = 147, size.per.fold = 14.7
    # k = 1, size.per.fold*k = 14.7, 1:15 , n = 15
    # k = 2, size.per.fold*k = 29.4, 16:29, n = 14
    # k = 3, size.per.fold*k = 44.1, 30:44, n = 15
  }
  for (i in 1:repeats){
    # Generate random fold vector
    set.seed(i)
    random.folds.i <- sample(kfold.vector, nrow(data), replace = FALSE)
    for (j in 1:kfold) {
      # Define Iter
      output.df$Iter[((i-1)*kfold)+j] <- paste('Repeat', i, 'Fold', j)
      # Define train/test sets 1:kfold
      test.idx <- c(which(random.folds.i == j))
      test.set <- data[test.idx,]
      train.set <- data[-test.idx,]
      # Train model
      model.i <- lm(formula, train.set)
      # Calculate rsquared for this fold
      output.df$Rsquared.train[((i-1)*kfold)+j] <- summary(model.i)$r.squared
      # Calculate RMSE for this fold
      output.df$RMSE.train[((i-1)*kfold)+j] <- RMSE(model.i$fitted.values, train.set[,yvar])
      # Train correlation
      output.df$Train.cor[((i-1)*kfold)+j] <- cor(train.set[,yvar], model.i$fitted.values)
      # Test correlation
      agepred <- predict(model.i, test.set)
      output.df$Test.cor[((i-1)*kfold)+j] <- cor(test.set[,yvar], agepred)
    }
    # if (i %in% seq(from = 0, to = 1000, by = 100)) {
    #   print(paste("Just completed ", i, "th iteration"))
    # }
  }
  return(output.df)
}

# Apply to 22 models
model.list <- list()
RMSE.exh <- matrix(0, nrow = 1000, ncol = top)
R.sqr.exh <- matrix(0, nrow = 1000, ncol = top)
Train.cor.exh <- matrix(0, nrow = 1000, ncol = top)
Test.cor.exh <- matrix(0, nrow = 1000, ncol = top)
for (i in 1:top){
  model_i <- my.repeat.cv(formula = all.models.exh$Formula[i],
                          data = data3,
                          kfold = 10, 
                          repeats = 100)
  model.list[[i]] <- model_i
  RMSE.exh[,i] <- model.list[[i]]$RMSE.train
  R.sqr.exh[,i] <- model.list[[i]]$Rsquared.train
  Train.cor.exh[,i] <- model.list[[i]]$Train.cor
  Test.cor.exh[,i] <- model.list[[i]]$Test.cor
  print(paste("Just completed ", i, "th iteration"))
}

which.min(apply(RMSE.exh, 2, mean))  # 26 #22
which.max(apply(R.sqr.exh, 2, mean)) # 26 #22
which.max(apply(Train.cor.exh, 2, mean)) # 26 #22
which.max(apply(Test.cor.exh, 2, mean)) # 12 #10
mean(apply(Test.cor.exh, 2, mean) ) #0.7580995 #0.7490471
final <- which.max(apply(Test.cor.exh, 2, mean)) #final # of parameters

fileName = paste(path_out, '20220526 RMSE for exh search models.csv',sep = '')
write.csv(RMSE.exh, 
          file = fileName,
          row.names = FALSE)
fileName = paste(path_out, '20220526 R.sqr for exh search models.csv',sep = '')
write.csv(R.sqr.exh, 
          file = fileName,
          row.names = FALSE)
fileName = paste(path_out, '20220526 Train set correlations.csv',sep = '')
write.csv(Train.cor.exh, 
          file = fileName,
          row.names = FALSE)
fileName = paste(path_out, '20220526 Test set correlations.csv',sep = '')
write.csv(Test.cor.exh, 
          file = fileName,
          row.names = FALSE)

# Plot
boxplot(Test.cor.exh, main = 'Test Correlation for 100x 10-fold CV',
        ylab = 'Pearson Correlation Coefficient', xlab = 'Parameter #')

#################################################################################
# (4) FINAL MODEL
#################################################################################
# Predict
#final_model <- lm(all.models.exh$Formula[final], data3)
final_model <- lm("declineRate~Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2", data3)
#data3 <- data3[-which(data3$INDDID == 111864),]
data3$SlopePred <- predict(final_model, data3)
data3$DeltaSlope <- data3$SlopePred - data3$declineRate


# Diagnostic Plots
plot(final_model,1)
plot(final_model,2)
plot(final_model,3)

# Generate residual plots
ggplot(data3, aes(x = Group, y = DeltaSlope)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Delta DRS Slope\n(11-parameter model)", 
       y = "Delta Slope", x = "Disease Group") +
  theme(legend.position = "right")

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(data3$SlopePred) #-2.227116 LOG -2.134599
max(data3$SlopePred) #1.944714 LOG 2.0283

cor(data3$declineRate,data3$SlopePred, method = "pearson")
#0.945008 LOG #0.9316612

dev.off()
par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = data3$SlopePred,
     xlim = c(-2.5,2.5), ylim = c(-2.5,2.5),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Predicted DRS Slope',
     cex.lab = 1.5)
     #main = 'Proteomic DRS Slope vs. Actual DRS Slope')
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.945")

#for paper figure
data3$color <- "black"
data3$color[which(data3$INDDID == 106049)] <- "red"
data3$color[which(data3$INDDID == 101287)] <- "orange"
data3$color[which(data3$INDDID == 118111)] <- "green"
data3$color[which(data3$INDDID == 108524)] <- "purple"
plot(x = data3$declineRate,
     y = data3$SlopePred,
     xlim = c(-2.5,2), ylim = c(-2.5,2),
     type = 'p', col = data3$color, pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Predicted DRS Slope',
     cex.lab = 1.5)
#main = 'Proteomic DRS Slope vs. Actual DRS Slope')
abline(a = 0, b = 1, col = 'black')
text(x = 0.5, y = -1.5, labels = "cor = 0.945")

all.models.exh$Formula[final]
#"declineRate~Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
#"declineRate~GALNT10+PDGFC+SELL+ICAM.2+SCF+ENAH+FCRL5+PTN.1+SDC4+VWC2" -- log transformation selection
final_model$coefficients[2:(final+1)]
#Dkk.4      PDGFC       SELL      PVALB       PPIB     ICAM.2       KYNU        SCF      FCRL5      PTN.1       VWC2         C2 
#-0.1866448 -0.6498774 -0.7193840  0.1093743  0.2256641  0.4492017  0.3863914 -0.4394916  0.4344448 -0.1207548 -0.4946895 -0.4055188 
#LOG
#Dkk.4      PDGFC       SELL      PVALB       PPIB     ICAM.2       KYNU        SCF      FCRL5      PTN.1       VWC2 
#-0.9266094 -2.3990706 -5.5184636  0.8505698  0.2763727  2.3712868  3.2610284 -4.0343502  1.1380155 -0.3669695 -2.8562718 
#C2 
#-2.1619666 
barplot(final_model$coefficients[2:(final+1)])


get_intercept <- function(x1, y1, slope) {
  # Point slope formula (y - y1) = slope(x - x1)
  y_intercept = slope * (- x1) + y1
  return(y_intercept)
}

#FOR FIGURE 2A 
#currinddid <- 106049 #(-2.3, 2013), red
#currinddid <- 101287  #(-1.3, 2016), orange
#currinddid <- 118111 #(-0.71, 2013), dark green
currinddid <- 108524 #(-0.36, 2014), purple
color <- "purple"
patientind <- which(pd_drs.dat$INDDID == currinddid)
plot(decimal_date(as.Date(pd_drs.dat$TestDate[patientind],"%m/%d/%Y")), pd_drs.dat$DRSTotalAge[patientind],
     main = paste("Patient ID", currinddid),
     ylab = "DRS", xlab = "Year", cex.lab = 1.3,
     pch = 21, bg = "black",
     ylim = c(2,18), xlim = c(2013, 2023))
abline(a=data3$bestFitIntercept[which(data3$INDDID == currinddid)],
       b=data3$declineRate[which(data3$INDDID == currinddid)],
       col = "black", lwd = 1.5, lty = 1)
text(x = 2020, y = 17.3, col = "black",labels = paste("actual slope =",round(data3$declineRate[which(data3$INDDID == currinddid)],2)))
pred_ypoint <-  data3$declineRate[which(data3$INDDID == currinddid)] * decimal_date(as.Date(pd_drs.dat$TestDate[patientind],"%m/%d/%Y"))[1] + data3$bestFitIntercept[which(data3$INDDID == currinddid)]
pred_intercept <- get_intercept(decimal_date(as.Date(pd_drs.dat$TestDate[patientind],"%m/%d/%Y"))[1],
                                pred_ypoint,
                                data3$SlopePred[which(data3$INDDID == currinddid)])
abline(#a=data3$bestFitIntercept[which(data3$INDDID == 106049)],
       a = pred_intercept,
       b=data3$SlopePred[which(data3$INDDID == currinddid)],
       col = color, lwd = 1.5, lty = 1)
text(x = 2020, y = 15.8, col = color,labels = paste("predicted slope =",round(data3$SlopePred[which(data3$INDDID == currinddid)],2)))

#################################################################################
# (5) FORCING IN AGE, SEX, AGE + SEX
#################################################################################
ss_final <- names(final_model$coefficients[2:(final+1)])
idx_final <- match(ss_final, names(data3))
olink_final <- data3[,c(8,17,19,idx_final)]

#FORCE AGE 
#Exhaustive Search
age_fit.exh <- regsubsets(declineRate~., olink_final, nvmax = final+1,
                          force.in = c(2),
                          method = 'exhaustive', really.big = TRUE)
age_summary.exh <- summary(age_fit.exh)

# Graph Rsq, Cp, BIC using training data
par(mfrow=c(3,1), mar=c(2.5,4,0.5,1), mgp=c(1.5,0.5,0))
plot(age_summary.exh$rsq, xlab="Number of predictors",
     ylab="R-squared", col="green", type="p", pch=16)
plot(age_summary.exh$cp, xlab="Number of predictors",
     ylab="Cp", col="red", type="p", pch=16)
plot(age_summary.exh$bic, xlab="Number of predictors",
     ylab="BIC", col="blue", type="p", pch=16)

# Results of exh selection
all.models_age <- data.frame(Var.Num = 1:final,
                             Formula = character(final),
                             R.sqrd = numeric(final),
                             Cp = numeric(final),
                             BIC = numeric(final),
                             AIC = numeric(final),
                             Full.Cor = numeric(final))
for (i in 1:final){
  #Save the essentials
  all.models_age$R.sqrd[i] <- age_summary.exh$rsq[i]
  all.models_age$Cp[i] <- age_summary.exh$cp[i]
  all.models_age$BIC[i] <- age_summary.exh$bic[i]
  #Get the names
  names_i <- names(olink_final[age_summary.exh$which[i,]])[-1]
  formula_i <- paste("declineRate~", paste(names_i, collapse = "+"), sep = "")
  all.models_age$Formula[i] <- formula_i
  #Fit model
  model_i <- lm(formula_i, data3)
  #AIC
  all.models_age$AIC[i] <- AIC(model_i)
  #Get correlation between actual and predicted
  all.models_age$Full.Cor[i] <- cor(data3$declineRate, model_i$fitted.values,
                                    method = 'pearson')
}
all.models_age.exh <- all.models_age
View(all.models_age.exh)

#replace protein with age
# Predict
age_model <- lm(all.models_age.exh$Formula[final-1], data3)
age_SlopePred <- predict(age_model, data3)
age_DeltaSlope <- age_SlopePred - data3$declineRate

# Diagnostic Plots
plot(age_model,1)
plot(age_model,2)
plot(age_model,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(age_SlopePred) #-2.207707 #log -2.11613
max(age_SlopePred) #2.036584 #log 1.876798

cor(data3$declineRate,age_SlopePred, method = "pearson")
#0.9408435 #log 0.9260067

dev.off()  
par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = age_SlopePred,
     xlim = c(-2.5,2.5), ylim = c(-2.5,2.5),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "Forced in Age (replaced protein)")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.926")

all.models_age.exh$Formula[final-1]
#"declineRate~AGE..YRS+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
#"declineRate~AGE..YRS+Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+VWC2+C2" #ptn.1 removed
age_model$coefficients[1:final+1]
#    AGE..YRS        PDGFC         SELL        PVALB         PPIB       ICAM.2         KYNU          SCF 
#0.006263518 -0.700417204 -0.647383588  0.107069271  0.235551830  0.399742869  0.404071822 -0.460450962 
#FCRL5        PTN.1         VWC2           C2 
#0.435850841 -0.156635329 -0.612223436 -0.447995877 
barplot(age_model$coefficients[1:final+1])

# add age in with all proteins
# Predict
age_model2 <- lm(all.models_age.exh$Formula[final], data3)
age_SlopePred2 <- predict(age_model2, data3)
age_DeltaSlope2 <- age_SlopePred2 - data3$declineRate

# Diagnostic Plots
plot(age_model2,1)
plot(age_model2,2)
plot(age_model2,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(age_SlopePred2) #-2.203646 #log -2.137961
max(age_SlopePred2) #2.015238 #log 2.044061

cor(data3$declineRate,age_SlopePred2, method = "pearson")
#0.9466908 #log 0.9317385

par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = age_SlopePred2,
     xlim = c(-2.5,2), ylim = c(-2.5,2),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "Forced in Age (with all 12 proteins)")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.931")

all.models_age.exh$Formula[final]
#"declineRate~AGE..YRS+Sex+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
age_model2$coefficients[2:(final+2)]
#AGE..YRS        Dkk.4        PDGFC         SELL        PVALB         PPIB       ICAM.2         KYNU          SCF 
#0.007985061 -0.199020512 -0.645444284 -0.622785869  0.107322273  0.239165388  0.448034821  0.390820270 -0.424977256 
#FCRL5        PTN.1         VWC2           C2 
#0.433840089 -0.148741082 -0.501691846 -0.414214449
barplot(age_model2$coefficients[2:(final+2)])

#FORCE SEX
sex_final <- data3[,c(8,19,idx_final)]
#Exhaustive Search
sex_fit.exh <- regsubsets(declineRate~., sex_final, nvmax = final+1,
                          force.in = c(2),
                          method = 'exhaustive', really.big = TRUE)
sex_summary.exh <- summary(sex_fit.exh)

# Graph Rsq, Cp, BIC using training data
par(mfrow=c(3,1), mar=c(2.5,4,0.5,1), mgp=c(1.5,0.5,0))
plot(sex_summary.exh$rsq, xlab="Number of predictors",
     ylab="R-squared", col="green", type="p", pch=16)
plot(sex_summary.exh$cp, xlab="Number of predictors",
     ylab="Cp", col="red", type="p", pch=16)
plot(sex_summary.exh$bic, xlab="Number of predictors",
     ylab="BIC", col="blue", type="p", pch=16)

# Results of exh selection
all.models_sex <- data.frame(Var.Num = 1:final,
                             Formula = character(final),
                             R.sqrd = numeric(final),
                             Cp = numeric(final),
                             BIC = numeric(final),
                             AIC = numeric(final),
                             Full.Cor = numeric(final))
for (i in 1:final){
  #Save the essentials
  all.models_sex$R.sqrd[i] <- sex_summary.exh$rsq[i]
  all.models_sex$Cp[i] <- sex_summary.exh$cp[i]
  all.models_sex$BIC[i] <- sex_summary.exh$bic[i]
  #Get the names
  names_i <- names(sex_final[sex_summary.exh$which[i,]])[-1]
  formula_i <- paste("declineRate~", paste(names_i, collapse = "+"), sep = "")
  all.models_sex$Formula[i] <- formula_i
  #Fit model
  model_i <- lm(formula_i, data3)
  #AIC
  all.models_sex$AIC[i] <- AIC(model_i)
  #Get correlation between actual and predicted
  all.models_sex$Full.Cor[i] <- cor(data3$declineRate, model_i$fitted.values,
                                    method = 'pearson')
}
all.models_sex.exh <- all.models_sex
View(all.models_sex.exh)

# REPLACE PROTEIN
# Predict
sex_model <- lm(all.models_sex.exh$Formula[final-1], data3)
sex_SlopePred <- predict(sex_model, data3)
sex_DeltaSlope <- sex_SlopePred - data3$declineRate

# Diagnostic Plots
plot(sex_model,1)
plot(sex_model,2)
plot(sex_model,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(sex_SlopePred) # -2.232237 #log -2.149706
max(sex_SlopePred) # 1.977699 #log 2.090542

cor(data3$declineRate,sex_SlopePred, method = "pearson")
#0.9400296 #log 0.9237836

dev.off()
par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = sex_SlopePred,
     xlim = c(-2.5,2.5), ylim = c(-2.5,2.5),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "Forced in Sex (replaced protein)")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.924")

all.models_sex.exh$Formula[final-1]
#declineRate~Sex+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
sex_model$coefficients[2:(final+1)]
#(Intercept)        Sex1     GALNT10        SELL       PVALB      ICAM.2        KYNU         SCF       FCRL5 
#10.3969674  -0.1097016  -0.8830590  -0.5993314   0.1079669   0.8973836   0.3148789  -0.4978064   0.3752804 
#C2        FLT1        VWC2 
#-0.2912860  -1.1163233  -0.4905270 
barplot(sex_model$coefficients[2:(final+1)])

# with all proteins
# Predict
sex_model2 <- lm(all.models_sex.exh$Formula[final], data3)
sex_SlopePred2 <- predict(sex_model2, data3)
sex_DeltaSlope2 <- sex_SlopePred2 - data3$declineRate2

# Diagnostic Plots
plot(sex_model2,1)
plot(sex_model2,2)
plot(sex_model2,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(sex_SlopePred2) # -2.224301 #log -2.130426
max(sex_SlopePred2) # 1.944467 #log 2.024771

cor(data3$declineRate,sex_SlopePred2, method = "pearson")
#0.9450559 #log 0.9319674

par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = sex_SlopePred2,
     xlim = c(-2.5,2.5), ylim = c(-2.5,2.5),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "Forced in Sex (with all 12 proteins)")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.932")


sex_model2$coefficients[2:(final+2)]
#Sex1       Dkk.4       PDGFC        SELL       PVALB        PPIB      ICAM.2        KYNU 
#0.01966575 -0.19239851 -0.65130512 -0.70968806  0.10815162  0.22748992  0.44762258  0.38971990 
#SCF       FCRL5       PTN.1        VWC2          C2 
#-0.43903469  0.42855310 -0.11814681 -0.48953161 -0.40650919 
barplot(sex_model2$coefficients[2:(final+2)])

#FORCE AGE AND SEX
#Exhaustive Search
agesex_fit.exh <- regsubsets(declineRate~., olink_final, nvmax = final+2,
                             force.in = c(2,3),
                             method = 'exhaustive', really.big = TRUE)
agesex_summary.exh <- summary(agesex_fit.exh)

# Graph Rsq, Cp, BIC using training data
par(mfrow=c(3,1), mar=c(2.5,4,0.5,1), mgp=c(1.5,0.5,0))
plot(agesex_summary.exh$rsq, xlab="Number of predictors",
     ylab="R-squared", col="green", type="p", pch=16)
plot(agesex_summary.exh$cp, xlab="Number of predictors",
     ylab="Cp", col="red", type="p", pch=16)
plot(agesex_summary.exh$bic, xlab="Number of predictors",
     ylab="BIC", col="blue", type="p", pch=16)

# Results of exh selection
all.models_agesex <- data.frame(Var.Num = 1:final,
                                Formula = character(final),
                                R.sqrd = numeric(final),
                                Cp = numeric(final),
                                BIC = numeric(final),
                                AIC = numeric(final),
                                Full.Cor = numeric(final))
for (i in 1:final){
  #Save the essentials
  all.models_agesex$R.sqrd[i] <- agesex_summary.exh$rsq[i]
  all.models_agesex$Cp[i] <- agesex_summary.exh$cp[i]
  all.models_agesex$BIC[i] <- agesex_summary.exh$bic[i]
  #Get the names
  names_i <- names(olink_final[agesex_summary.exh$which[i,]])[-1]
  formula_i <- paste("declineRate~", paste(names_i, collapse = "+"), sep = "")
  all.models_agesex$Formula[i] <- formula_i
  #Fit model
  model_i <- lm(formula_i, data3)
  #AIC
  all.models_agesex$AIC[i] <- AIC(model_i)
  #Get correlation between actual and predicted
  all.models_agesex$Full.Cor[i] <- cor(data3$declineRate, model_i$fitted.values,
                                       method = 'pearson')
}
all.models_agesex.exh <- all.models_agesex
View(all.models_agesex.exh)

# Predict
agesex_model <- lm(all.models_agesex.exh$Formula[final-2], data3)
agesex_SlopePred <- predict(agesex_model, data3)
agesex_DeltaSlope <- agesex_SlopePred - data3$declineRate

# Diagnostic Plots
plot(agesex_model,1)
plot(agesex_model,2)
plot(agesex_model,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(agesex_SlopePred) #-2.231992 #log 2.107579
max(agesex_SlopePred) #2.101176 #log 1.912365

cor(data3$declineRate,agesex_SlopePred, method = "pearson")
#0.9317609 #log 0.9169344

par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = agesex_SlopePred,
     xlim = c(-2.5,2), ylim = c(-2.5,2),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "Forced in Age + Sex (replaced proteins)")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.912")

all.models_agesex.exh$Formula[final-2]
#"declineRate~AGE..YRS+Sex+PDGFC+SELL+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
agesex_model$coefficients[2:(final+1)]
#AGE..YRS         Sex1        PDGFC         SELL         PPIB       ICAM.2         KYNU          SCF        FCRL5 
#0.007197442 -0.013359434 -0.704249178 -0.565896556  0.330249426  0.374242748  0.435722100 -0.555812203  0.561192378 
#PTN.1         VWC2           C2 
#-0.142253024 -0.636249254 -0.456300902 
barplot(agesex_model$coefficients[2:(final+1)])

#all proteins
# Predict
agesex_model2 <- lm(all.models_agesex.exh$Formula[final], data3)
agesex_SlopePred2 <- predict(agesex_model2, data3)
agesex_DeltaSlope2 <- agesex_SlopePred2 - data3$declineRate

# Diagnostic Plots
plot(agesex_model2,1)
plot(agesex_model2,2)
plot(agesex_model2,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(agesex_SlopePred2) #-2.20224 #-2.133356
max(agesex_SlopePred2) #2.014603 #2.037898

cor(data3$declineRate,agesex_SlopePred2, method = "pearson")
#0.9467056 #log0.9320194

par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = agesex_SlopePred2,
     xlim = c(-2.5,2), ylim = c(-2.5,2),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "Forced in Age + Sex (with all 12 proteins)")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.932")

agesex_model2$coefficients[2:(final+3)]
#AGE..YRS         Sex1        Dkk.4        PDGFC         SELL        PVALB         PPIB       ICAM.2         KYNU          SCF        FCRL5 
#0.007928736  0.010977378 -0.202144908 -0.646272478 -0.618054992  0.106654269  0.240089327  0.447161578  0.392646981 -0.424824588  0.430555632 
#PTN.1         VWC2           C2 
#-0.147087925 -0.498763352 -0.414705971 
barplot(agesex_model2$coefficients[2:(final+3)])

#no No Dkk.4
newForm = "declineRate~PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
# Predict
newForm_model <- lm(newForm, data3)
newForm_SlopePred <- predict(newForm_model, data3)
newForm_DeltaSlope <- newForm_SlopePred - data3$declineRate

# Diagnostic Plots
plot(newForm_model,1)
plot(newForm_model,2)
plot(newForm_model,3)

# Stats
min(data3$declineRate) #-2.363552
max(data3$declineRate) #1.792649
min(newForm_SlopePred) #-2.226167 #log -2.147772
max(newForm_SlopePred) #1.979458 #log 2.090779

cor(data3$declineRate,newForm_SlopePred, method = "pearson")
#0.9397879 #log0.9236696

par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = newForm_SlopePred,
     xlim = c(-2.5,2), ylim = c(-2.5,2),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = "No Dkk.4")
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -1.5, labels = "cor = 0.924")

newForm
#"declineRate~PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2"
newForm_model$coefficients[2:final]
#PDGFC       SELL      PVALB       PPIB     ICAM.2       KYNU        SCF      FCRL5      PTN.1       VWC2         C2 
#-0.7012240 -0.7229692  0.1087132  0.2249961  0.4030582  0.3998955 -0.4702370  0.4362322 -0.1339941 -0.6011915 -0.4394121 
barplot(newForm_model$coefficients[2:final])

#################################################################################
# (5) COMPARE PD WITH NC AND AD 
#################################################################################
setwd("/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project")

#PD MOCA 
# Load Data
moca_pd.dat <- read.csv("01_pd_moca.csv")

# Define Decline Rate 
pdPatients <- unique(moca_pd.dat$INDDID)
pdPatients <- pdPatients[-which(pdPatients == 111864)] 
pdnum <- length(pdPatients) # number of PD patients (49)
data1_pd <- data.frame(INDDID = numeric(pdnum),
                       sampleDate = character(pdnum),
                       numMOCA = numeric(pdnum), 
                       declineRate = numeric(pdnum),
                       declineSpeed = character(pdnum),
                       fastProgression = numeric(pdnum))
path_out3 <- "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/figures/1_pd moca/"

for (i in 1:pdnum) {
  data1_pd$INDDID[i] <- pdPatients[i]
  
  #when olink sample was taken 
  olinkInd <- match(pdPatients[i], olink_pd.dat$INDDID)
  data1_pd$sampleDate[i] <- olink_pd.dat$SampleDate[olinkInd]
  sampleDate <- decimal_date(as.Date(data1_pd$sampleDate[i], "%m/%d/%Y"))  
  
  #moca scores
  allIndex <- which(moca_pd.dat$INDDID %in% pdPatients[i])
  allMOCA <- moca_pd.dat$MoCATotal[allIndex]
  allYears <- (as.Date(moca_pd.dat$TestDate[allIndex], "%m/%d/%Y"))
  allYears <- decimal_date(allYears)
  
  #remove MOCA taken before sample
  dateDiff <- allYears - sampleDate
  toRemoveIndex <- which(dateDiff < 0)
  if(length(toRemoveIndex) != 0) {
    allYears <- allYears[-toRemoveIndex]
    allMOCA <- allMOCA[-toRemoveIndex]
  }
  
  data1_pd$numMOCA[i] <- length(allMOCA)
  
  if(length(allMOCA)<2){
    data1_pd$declineRate[i] <- NaN 
    next 
  }
  
  # find best fitting line
  bestLine <- lm(allMOCA ~ allYears) 
  
  #set decline rate to slope of best fitting line
  slope <- coef(bestLine)[2]
  data1_pd$declineRate[i] <- coef(bestLine)[2] 
  
  jpeg(file=paste(path_out3, data1_pd$INDDID[i],"bestFitLineDRS.jpeg"))
  plot(allYears, allMOCA,
       main = paste("INDDID", data1_pd$INDDID[i], ": Slope", slope),
       ylab = "MOCA Score", xlab = "Date")
  abline(bestLine)
  dev.off()  
  
  if (data1_pd$declineRate[i] < -0.5) {
    data1_pd$declineSpeed[i] <- "Fast"
    data1_pd$fastProgression[i] <- 1
  }
  else {
    data1_pd$declineSpeed[i] <- "Slow"
    data1_pd$fastProgression[i] <- 0
  }
}

data1_pd <- data1_pd[data1_pd$numMOCA>1,]
data1_pd <- data1_pd %>% arrange(declineRate)

slowdata_pd <- subset(data1_pd, declineSpeed == "Slow")
nrow(slowdata_pd) #26
fastdata_pd <- subset(data1_pd, declineSpeed == "Fast")
nrow(fastdata_pd) #23

# Predict
data2_pd <- merge(data1_pd, olink_pd.dat)
finalModel_pd <- lm(paste(all.models.exh$Formula[12], "+Sex+AGE..YRS"), data2_pd)
data2_pd$SlopePred <- predict(finalModel_pd, data2_pd)
data2_pd$DeltaSlope <- data2_pd$SlopePred - data2_pd$declineRate

cor(data2_pd$declineRate,data2_pd$SlopePred, method = "pearson")
#0.7372228

par(mfrow=c(1,1))
plot(x = data2_pd$declineRate,
     y = data2_pd$SlopePred,
     xlim = c(-5,4), ylim = c(-5,4),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual MOCA Slope',
     ylab = 'Proteomic MOCA Slope')
abline(a = 0, b = 1, col = 'blue')
text(x = 1, y = -4, labels = "cor = 0.74")

finalModel_pd$coefficients
#(Intercept)     TNFAIP8       Dkk.4       PDGFC     GALNT10      ICAM.2        KYNU    SERPINA7        ENAH         SCF 
#11.24201741  0.02557266  0.19506158 -1.00472833  0.35572286 -0.16240522  0.53045133 -0.42364925  0.67457404 -1.28839081 
#FCRL5          C2        FLT1        VWC2     SexMale    AGE..YRS 
#0.11057368  0.39172219 -0.42791226 -0.21228032  0.46490919 -0.04396749 
barplot(finalModel_pd$coefficients[2:14])

#Binary Classifier for full data set 
# Let's train the model using the 12 proteins
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
data2_pd$fastProgression <- factor(data2_pd$fastProgression, levels = c(0, 1), labels=c("Slow","Fast"))
data2_pd$Sex <- factor(data2_pd$Sex, levels=c("Female","Male"), labels=c(0,1))
model_all_pd <- train(as.factor(fastProgression) ~
                     AGE..YRS +
                     Sex + 
                     Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
                   method = "glm",
                   data = data2_pd,
                   trControl = train_control,
                   na.action=na.pass,
                   family = "binomial")
summary(model_all_pd)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control_all_pd <- train(as.factor(fastProgression) ~
                             AGE..YRS +
                             Sex,
                           method = "glm",  
                           data = data2_pd,
                           trControl = train_control,
                           na.action=na.pass,
                           family = "binomial")

# Predict using test set 
predict_pd <- predict(model_all_pd, newdata = data2_pd, type= "prob")
predictControl_pd <- predict(model.control_all_pd, newdata = data2_pd, type= "prob")

roc_pd <- roc(data2_pd$fastProgression, as.numeric(predict_pd$'Fast'))
rocControl_pd <- roc(data2_pd$fastProgression, as.numeric(predictControl_pd$'Fast'))

#95% confidence intervals for AUC and accuracy
ci.auc(roc_pd)
predict_pd <- predict(model_all_pd, newdata = data2_pd)
xtab <- table(predict_pd, data2_pd$fastProgression)
confusionMatrix(xtab)

ci.auc(rocControl_pd)
predictControl_pd <- predict(model.control_all_pd, newdata = data2_pd)
xtabcontrol <- table(predictControl_pd, data2_pd$fastProgression)
confusionMatrix(xtabcontrol)

par(mar = c(5, 5, 5, 5))
plot(roc_pd, col = "blue", legacy.axes = TRUE, print.auc=TRUE, print.auc.y=.5)
plot(rocControl_pd, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .4)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "ROC for Discovery Cohort (MOCA)", line = 3.5)

#NEGATIVE CONTROL
# Load Data
olink_nc.dat <- olink.dat[olink.dat$Group == "NC",]
olink_nc.dat <- olink_nc.dat[!is.na(olink_nc.dat$AGE..YRS),]
moca_nc.dat <- read.csv("02_nc_moca.csv")

ncToRemove <- which(is.na(moca_nc.dat$MOCATOTS)) #remove empty MOCA
moca_nc.dat <- moca_nc.dat[-ncToRemove,]

# Define Decline Rate 
ncPatients <- unique(moca_nc.dat$INDDID)
ncnum <- length(ncPatients) # number of NC patients (50)
data1_nc <- data.frame(INDDID = numeric(ncnum),
                       sampleDate = character(ncnum),
                      numMOCA = numeric(ncnum), 
                      declineRate = numeric(ncnum),
                      declineSpeed = character(ncnum),
                      fastProgression = numeric(ncnum))
path_out2 <- "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/figures/2_nc"

for (i in 1:ncnum) {
  data1_nc$INDDID[i] <- ncPatients[i]
  
  #when olink sample was taken 
  olinkInd <- match(ncPatients[i], olink_nc.dat$INDDID)
  data1_nc$sampleDate[i] <- olink_nc.dat$SampleDate[olinkInd]
  sampleDate <- decimal_date(as.Date(data1_nc$sampleDate[i], "%m/%d/%Y"))  
  
  #moca scores
  allIndex <- which(moca_nc.dat$INDDID %in% ncPatients[i])
  allMOCA <- moca_nc.dat$MOCATOTS[allIndex]
  allYears <- (as.Date(moca_nc.dat$TestDate[allIndex], "%m/%d/%Y"))
  allYears <- decimal_date(allYears)
  
  #remove DRS taken before sample
  dateDiff <- allYears - sampleDate
  toRemoveIndex <- which(dateDiff < 0)
  if(length(toRemoveIndex) != 0) {
    allYears <- allYears[-toRemoveIndex]
    allMOCA <- allMOCA[-toRemoveIndex]
  }
  
  data1_nc$numMOCA[i] <- length(allMOCA)
  
  if(length(allMOCA)<2){
    data1_nc$declineRate[i] <- NaN 
    next 
  }
  
  # find best fitting line
  bestLine <- lm(allMOCA ~ allYears) 
  
  #set decline rate to slope of best fitting line
  slope <- coef(bestLine)[2]
  data1_nc$declineRate[i] <- coef(bestLine)[2] 
  
  jpeg(file=paste(path_out2, data1_nc$INDDID[i],"bestFitLineDRS.jpeg"))
  plot(allYears, allMOCA,
       main = paste("INDDID", data1_nc$INDDID[i], ": Slope", slope),
       ylab = "DRS Score", xlab = "Date")
  abline(bestLine)
  dev.off()  
  
  if (data1_nc$declineRate[i] < -0.5) {
    data1_nc$declineSpeed[i] <- "Fast"
    data1_nc$fastProgression[i] <- 1
  }
  else {
    data1_nc$declineSpeed[i] <- "Slow"
    data1_nc$fastProgression[i] <- 0
  }
}

data1_nc <- data1_nc[data1_nc$numMOCA>1,]
data1_nc <- data1_nc %>% arrange(declineRate)

slowdata_nc <- subset(data1_nc, declineSpeed == "Slow")
nrow(slowdata_nc) #30 slow progressing patients #42
fastdata_nc <- subset(data1_nc, declineSpeed == "Fast")
nrow(fastdata_nc) #11 fast progressing patients #8

# Predict
data2_nc <- merge(data1_nc, olink_nc.dat)
finalModel_nc <- lm(paste(all.models.exh$Formula[13], "+Sex+AGE..YRS"), data2_nc)
data2_nc$SlopePred <- predict(finalModel_nc, data2_nc)
data2_nc$DeltaSlope <- data2_nc$SlopePred - data2_nc$declineRate

# Diagnostic Plots
plot(finalModel_nc,1)
plot(finalModel_nc,2)
plot(finalModel_nc,3)

# Generate residual plots
ggplot(data2_nc, aes(x = Group, y = DeltaSlope)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Delta DRS Slope\n(13-parameter model)", 
       y = "Delta Slope", x = "Disease Group") +
  theme(legend.position = "right")

# Stats
min(data2_nc$declineRate) #-4.905914
max(data2_nc$declineRate) #3.075843
min(data2_nc$SlopePred) #-2.237577
max(data2_nc$SlopePred) # 1.750873

cor(data2_nc$declineRate,data2_nc$SlopePred, method = "pearson")
#0.6779559

par(mfrow=c(1,1))
plot(x = data2_nc$declineRate,
     y = data2_nc$SlopePred,
     xlim = c(-5,3), ylim = c(-5,3),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope')
abline(a = 0, b = 1, col = 'blue')
text(x = 1, y = -4, labels = "cor = 0.678")

all.models.exh$Formula[13]
#"declineRate~TNFAIP8+Dkk.4+PDGFC+GALNT10+ICAM.2+KYNU+SERPINA7+ENAH+SCF+FCRL5+C2+FLT1+VWC2"
finalModel_nc$coefficients
#(Intercept)        ENAH     S100A16       EGFL7        PRCP       PVALB        RFNG       CRHBP       PDGFC VWC2 
#-2.89755584  0.59996561 -0.16992219  0.27550721  0.25566150  0.05573779 -0.27957598  0.72311097 -0.35408227 VWC2 
#-0.41178338 
barplot(finalModel_nc$coefficients)

#AD MOCA 
# Load Data
moca_ad.dat <- read.csv("03_ad_moca.csv")
olink_ad.dat <- olink.dat[olink.dat$Group == "AD",]
olink_ad.dat <- olink_ad.dat[!is.na(olink_ad.dat$AGE..YRS),]

adToRemove <- which(is.na(moca_ad.dat$MOCATOTS)) #remove empty MOCA
moca_ad.dat <- moca_ad.dat[-adToRemove,]

# Define Decline Rate 
adPatients <- unique(moca_ad.dat$INDDID)
adnum <- length(adPatients) # number of ad patients (47)
data1_ad <- data.frame(INDDID = numeric(adnum),
                       sampleDate = character(adnum),
                       numMOCA = numeric(adnum), 
                       declineRate = numeric(adnum),
                       declineSpeed = character(adnum),
                       fastProgression = numeric(adnum))
path_out4 <- "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/figures/3_ad moca/"

for (i in 1:adnum) {
  data1_ad$INDDID[i] <- adPatients[i]
  
  #when olink sample was taken 
  olinkInd <- match(adPatients[i], olink_ad.dat$INDDID)
  data1_ad$sampleDate[i] <- olink_ad.dat$SampleDate[olinkInd]
  sampleDate <- decimal_date(as.Date(data1_ad$sampleDate[i], "%m/%d/%Y"))  
  
  #moca scores
  allIndex <- which(moca_ad.dat$INDDID %in% adPatients[i])
  allMOCA <- moca_ad.dat$MOCATOTS[allIndex]
  allYears <- (as.Date(moca_ad.dat$TestDate[allIndex], "%m/%d/%Y"))
  allYears <- decimal_date(allYears)
  
  #remove DRS taken before sample
  dateDiff <- allYears - sampleDate
  toRemoveIndex <- which(dateDiff < 0)
  if(length(toRemoveIndex) != 0) {
    allYears <- allYears[-toRemoveIndex]
    allMOCA <- allMOCA[-toRemoveIndex]
  }
  
  if(i == 18){
    data1_ad$numMOCA[i] <- 1
    data1$declineRate[i] <- NaN 
    next     
  }
  data1_ad$numMOCA[i] <- length(allMOCA)
  
  if(length(allMOCA)<2){
    data1_ad$declineRate[i] <- NaN 
    next 
  }
  
  # find best fitting line
  bestLine <- lm(allMOCA ~ allYears) 
  
  #set decline rate to slope of best fitting line
  slope <- coef(bestLine)[2]
  data1_ad$declineRate[i] <- coef(bestLine)[2] 
  
  jpeg(file=paste(path_out4, data1_ad$INDDID[i],"bestFitLineDRS.jpeg"))
  plot(allYears, allMOCA,
       main = paste("INDDID", data1_ad$INDDID[i], ": Slope", slope),
       ylab = "DRS Score", xlab = "Date")
  abline(bestLine)
  dev.off()  
  
  if (data1_ad$declineRate[i] < -0.5) {
    data1_ad$declineSpeed[i] <- "Fast"
    data1_ad$fastProgression[i] <- 1
  }
  else {
    data1_ad$declineSpeed[i] <- "Slow"
    data1_ad$fastProgression[i] <- 0
  }
}

data1_ad <- data1_ad[data1_ad$numMOCA>1,]
data1_ad <- data1_ad %>% arrange(declineRate)

slowdata_ad <- subset(data1_ad, declineSpeed == "Slow")
nrow(slowdata_ad) #9
fastdata_ad <- subset(data1_ad, declineSpeed == "Fast")
nrow(fastdata_ad) #23

# Predict
data2_ad <- merge(data1_ad, olink_ad.dat)
data2_ad_updated <- data2_ad[-c(23,31),]
finalModel_ad <- lm(paste(all.models.exh$Formula[13]), data2_ad_updated)
data2_ad_updated$SlopePred <- predict(finalModel_ad, data2_ad_updated)
data2_ad_updated$DeltaSlope <- data2_ad_updated$SlopePred - data2_ad_updated$declineRate

# Diagnostic Plots
plot(finalModel_ad,1)
plot(finalModel_ad,2)
plot(finalModel_ad,3)

# Generate residual plots
ggplot(data2_ad_updated, aes(x = Group, y = DeltaSlope)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Delta DRS Slope\n(13-parameter model)", 
       y = "Delta Slope", x = "Disease Group") +
  theme(legend.position = "right")

# Stats
min(data2_ad_updated$declineRate) #-7.038567
max(data2_ad_updated$declineRate) #2.991803
min(data2_ad_updated$SlopePred) #-6.519446
max(data2_ad_updated$SlopePred) # 0.5619672

cor(data2_ad_updated$declineRate,data2_ad_updated$SlopePred, method = "pearson")
#0.7337011

par(mfrow=c(1,1))
plot(x = data2_ad_updated$declineRate,
     y = data2_ad_updated$SlopePred,
     xlim = c(-8,3), ylim = c(-8,3),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope')
abline(a = 0, b = 1, col = 'blue')
text(x = 1, y = -6, labels = "cor = 0.546")

all.models.exh$Formula[13]
#"declineRate~TNFAIP8+Dkk.4+PDGFC+GALNT10+ICAM.2+KYNU+SERPINA7+ENAH+SCF+FCRL5+C2+FLT1+VWC2"
finalModel_ad$coefficients
#(Intercept)     TNFAIP8       Dkk.4       PDGFC     GALNT10      ICAM.2        KYNU    SERPINA7        ENAH 
#8.37653566 -0.03317852  0.14434664 -0.90579087  0.35258433 -0.12529244  0.68290315 -0.58482546  0.54496129 
#SCF       FCRL5          C2        FLT1        VWC2 
#-1.10556385  0.18798625  0.40802434 -0.57354692 -0.57236985 
barplot(finalModel_ad$coefficients)



#################################################################################
# (6) FAST VS. SLOW PROGRESSION CLASSIFIER 
#################################################################################

#Binary Classifier for full data set 
# Let's train the model using the 12 proteins
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
model_all <- train(as.factor(fastProgression) ~
                 AGE..YRS +
                 Sex + 
                 Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
               method = "glm",
               data = data3,
               trControl = train_control,
               na.action=na.pass,
               family = "binomial")
summary(model_all)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control_all <- train(as.factor(fastProgression) ~
                         AGE..YRS +
                         Sex,
                       method = "glm",  
                       data = data3,
                       trControl = train_control,
                       na.action=na.pass,
                       family = "binomial")

# Predict using test set 
predict <- predict(model_all, newdata = data3, type= "prob")
predictControl <- predict(model.control_all, newdata = data3, type= "prob")

roc <- roc(data3$fastProgression, as.numeric(predict$'1'))
rocControl <- roc(data3$fastProgression, as.numeric(predictControl$'1'))

#95% confidence intervals for AUC and accuracy
ci.auc(roc) #auc 1-1 
predict <- predict(model_all, newdata = data3)
xtab <- table(predict, data3$fastProgression)
confusionMatrix(xtab) #accuracy 1 (0.9275-1)

ci.auc(rocControl) #auc (0.4007-0.7327)
predictControl <- predict(model.control_all, newdata = data3)
xtabcontrol <- table(predictControl, data3$fastProgression)
confusionMatrix(xtabcontrol) #accuracy0.5714 (0.4221, 0.7118)

par(mar = c(5, 5, 5, 5))
plot(roc, col = "blue", legacy.axes = TRUE, print.auc=TRUE, print.auc.y=.5)
plot(rocControl, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .4)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "ROC for Discovery Cohort", line = 3.5)
# AUC: 1.000, control AUC: 0.567


#Determine train/test divisions (75% train, 25% test)
set.seed(140)
trainIndex <- createDataPartition(data3$fastProgression,p=0.75,list=FALSE)
trainSet <- data3[trainIndex,] #training data (75% of data)
testSet <- data3[-trainIndex,] #testing data (25% of data)

# Let's train the model using the 12 proteins
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
model <- train(as.factor(fastProgression) ~
                 AGE..YRS +
                 Sex + 
                 Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
               method = "glm",
               data = trainSet,
               trControl = train_control,
               na.action=na.pass,
               family = "binomial")
summary(model)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control <- train(as.factor(fastProgression) ~
                         AGE..YRS +
                         Sex,
                       method = "glm",  
                       data = trainSet,
                       trControl = train_control,
                       na.action=na.pass,
                       family = "binomial")

# Predict using test set 
predict <- predict(model, newdata = testSet, type= "prob")
predictControl <- predict(model.control, newdata = testSet, type= "prob")

roc <- roc(testSet$fastProgression, as.numeric(predict$'1'))
rocControl <- roc(testSet$fastProgression, as.numeric(predictControl$'1'))

par(mar = c(5, 5, 5, 5))
plot(roc, col = "blue", legacy.axes = TRUE, print.auc=TRUE)
plot(rocControl, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .3)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "75% Train, 25% Test", line = 3.5)
# AUC: 0.814, control AUC: 0.543

#Determine train/test divisions (50% train, 50% test)
set.seed(1410)
trainIndex_50 <- createDataPartition(data3$fastProgression,p=0.50,list=FALSE)
trainSet_50 <- data3[trainIndex_50,] #training data (50% of data)
testSet_50 <- data3[-trainIndex_50,] #testing data (50% of data)

# Let's train the model using the 13 proteins
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
model_50 <- train(as.factor(fastProgression) ~
                 AGE..YRS +
                 Sex + 
                   Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
               method = "glm",
               data = trainSet_50,
               trControl = train_control,
               na.action=na.pass,
               family = "binomial")
summary(model)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control_50 <- train(as.factor(fastProgression) ~
                         AGE..YRS +
                         Sex,
                       method = "glm",  
                       data = trainSet_50,
                       trControl = train_control,
                       na.action=na.pass,
                       family = "binomial")

# Predict using test set 
predict_50 <- predict(model_50, newdata = testSet_50, type= "prob")
predictControl_50 <- predict(model.control_50, newdata = testSet_50, type= "prob")

roc_50 <- roc(testSet_50$fastProgression, as.numeric(predict_50$'1'))
rocControl_50 <- roc(testSet_50$fastProgression, as.numeric(predictControl_50$'1'))

par(mar = c(5, 5, 5, 5))
plot(roc_50, col = "blue", legacy.axes = TRUE, print.auc=TRUE)
plot(rocControl_50, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .3)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "50% Train, 50% Test", line = 3.5)
#AUC = 0.556, control AUC = 0.507


#Determine train/test divisions (80% train, 20% test)
set.seed(1405)
trainIndex_80 <- createDataPartition(data3$fastProgression,p=0.80,list=FALSE)
trainSet_80 <- data3[trainIndex_80,] #training data (80% of data)
testSet_80 <- data3[-trainIndex_80,] #testing data (20% of data)

# Let's train the model using the 12 proteins
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
model_80 <- train(as.factor(fastProgression) ~
                    AGE..YRS +
                    Sex + 
                    Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
                  method = "glm",
                  data = trainSet_80,
                  trControl = train_control,
                  na.action=na.pass,
                  family = "binomial")
summary(model)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control_80 <- train(as.factor(fastProgression) ~
                            AGE..YRS +
                            Sex,
                          method = "glm",  
                          data = trainSet_80,
                          trControl = train_control,
                          na.action=na.pass,
                          family = "binomial")

# Predict using test set 
predict_80 <- predict(model_80, newdata = testSet_80, type= "prob")
predictControl_80 <- predict(model.control_80, newdata = testSet_80, type= "prob")

roc_80 <- roc(testSet_80$fastProgression, as.numeric(predict_80$'1'))
rocControl_80 <- roc(testSet_80$fastProgression, as.numeric(predictControl_80$'1'))

par(mar = c(5, 5, 5, 5))
plot(roc_80, col = "blue", legacy.axes = TRUE, print.auc=TRUE)
plot(rocControl_80, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .3)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "80% Train, 20% Test", line = 3.5)
#auc 0.944, 0.667

#################################################################################
# (7) SURVIVAL ANALYSIS
#################################################################################

setwd("/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/3_survivalanalysis")
survival.dat <- read.csv("20220531_survival.csv")

for (i in 1:nrow(survival.dat)){
  currentINDDID <- survival.dat$INDDID[i]
  ind <- which(data3$INDDID == currentINDDID)
  survival.dat$Sex[i] <- data3$Sex[ind]
  survival.dat$AGE..YRS[i] <- data3$AGE..YRS[ind]
  survival.dat$declineRate[i] <- data3$declineRate[ind]
  survival.dat$declineSpeed[i] <- data3$declineSpeed[ind]
  survival.dat$fastProgression[i] <- data3$fastProgression[ind]
  survival.dat$slopePred[i] <- data3$SlopePred[ind]
  if (survival.dat$slopePred[i] < -0.5){
    survival.dat$predFastProgression[i] <- 1
  }
  else{
    survival.dat$predFastProgression[i] <- 0
  }
}

survival.dat$Sex <- factor(survival.dat$Sex, levels=c(1,2), labels=c(0,1)) #Male = 1, Female = 0
survival.dat$Sex <- factor(survival.dat$Sex, levels=c(0,1), labels=c("Female", "Male"))

# Create survival object
survive_object = Surv(survival.dat$YearsSinceSample, survival.dat$ClinicalConversion)

mykm1 <- survfit(Surv(YearsSinceSample, ClinicalConversion) ~ predFastProgression, data = survival.dat)
dev.off()
ggsurvplot(mykm1, data = survival.dat, pval = TRUE, 
           xlim = c(0, 11), break.y.by = 0.2,
           palette = c('#1CBDC3','#EF4438'),
           breaks = 0:11, 
           title = " Clinical Conversion",
           legend.title="Group",ylab = "%Cognitively Normal",xlab = "Years",
           censor.size=0.1, size = 1,ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
           legend.labs=c("Predicted Slow","Predicted Fast"),
           font.legend = list(size = 10, color = "black", face = "bold"),
           font.main = list(size= 20, color = "black", face="bold"),
           font.x = list(size= 13, color = "black", face="bold"),
           font.y = list(size= 13, color = "black", face="bold"),
           legend = "right"
)

#################################################################################
# (8) MENDELIAN RANDOMIZATION
#################################################################################
setwd("/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/2_mendelian randomization/")
MRdata <- read.csv("20220624_allMR.csv")

MR_DRS_data <- merge(data2, MRdata, by = "INDDID")

proteinlist <- c("Dkk.4","Dkk.4","PDGFC","PDGFC","SELL","PVALB","PPIB",
                 "ICAM.2","KYNU","SCF","SCF","FCRL5","PTN.1","VWC2")
snplist <- c("rs10677223","rs12678359","rs55947119","rs17035367","rs4987358","rs4821544",
             "rs73453022","rs635634","rs3768844","rs6065904","rs3782181","rs914724",
             "rs1431093","rs740083")
baseilist <- c("A","T","C","A","G","T","T","C","A","G","A","A","C","C")
baseflist <- c("G","C","T","G","T","C","C","T","G","A","C","G","A","T")

MR_mixedeffects <- function(protein, snp, basei, basef) {
  currcolname <- paste(snp,basef,sep="_")
  MR_DRS_data[currcolname]<- factor(MR_DRS_data[[currcolname]], 
                                     levels=c(0,1,2), 
                                     labels = c(paste(basei,basei,sep="/"), #A/A
                                                paste(basei,basef,sep="/"), #A/G
                                                paste(basef,basef,sep="/"))) #G/G
  #MRdata[currcolname] <- factor(MRdata[[currcolname]], 
  #%                              levels=c(0,1,2), 
  ##                              labels = c(paste(basei,basei,sep="/"), #A/A
  #                                         paste(basei,basef,sep="/"), #A/G
  #                                         paste(basef,basef,sep="/"))) #G/G
  myformula <- paste("DRSTotalAge~", paste(currcolname,"yearsToVisit + AgeatTest + Sex + DRS.initial + (1 | INDDID)",sep="*"))
  
  MR_DRS <- lmer(myformula, data = MR_DRS_data, REML = F) #no disease duration 

  p1 = plot_model(MR_DRS, type = "pred", terms = c("yearsToVisit", currcolname),
                  title = toString(protein),
                  axis.title = c("Years Since Plasma Sample","Cognition (DRS)"),
                  axis.labels = NULL, legend.title = snp,
                  wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                  ci.lvl = 0.95, se = NULL, colors = "Set1", order.terms = c(3,1,2)) + theme_classic() +
    theme(legend.position = c(0.15,0.25), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
          title = element_text(size = 18),plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
          axis.text = element_text(size = 16, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
  print(p1)
  return(MRdata)
}

for(i in 1:length(proteinlist)){
  MRdata <- MR_mixedeffects(protein = proteinlist[i], 
                                      snp = snplist[i], 
                                      basei = baseilist[i], 
                                      basef = baseflist[i])
  
}

#PDGFC: combine C/C and C/T
MR_DRS_data$rs55947119_T <- factor(MR_DRS_data$rs55947119_T, 
                                  levels=c(0,1,2), 
                                  labels = c("C/C + C/T", "C/C + C/T", "T/T"))
#MRdata$rs55947119_T <- factor(MRdata$rs55947119_T, 
#                              levels=c(0,1,2), 
#                              labels = c("C/C + C/T", "C/C + C/T", "T/T"))
MR_DRS <- lmer(DRSTotalAge~rs55947119_T*yearsToVisit +
               AgeatTest + Sex + DRS.initial + (1 | INDDID), data = MR_DRS_data, REML = F) #no disease duration 
p1 = plot_model(MR_DRS, type = "pred", terms = c("yearsToVisit", "rs55947119_T"),
                title = "PDGFC",
                axis.title = c("Years Since Plasma Sample","Cognition (DRS)"),
                axis.labels = NULL, legend.title = "rs55947119",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set1", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.18,0.20), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 18),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 16, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
print(p1)

#SELL: combine G/T and T/T 
MR_DRS_data$rs4987358_T <- factor(MR_DRS_data$rs4987358_T, 
                                   levels=c(0,1,2), 
                                   labels = c("G/G", "G/T + T/T", "G/T + T/T"))
MRdata$rs4987358_T <- MR_DRS_data$rs4987358_T
MR_DRS <- lmer(DRSTotalAge~rs4987358_T*yearsToVisit +
                 AgeatTest + Sex + DRS.initial + (1 | INDDID), data = MR_DRS_data, REML = F) #no disease duration 
p1 = plot_model(MR_DRS, type = "pred", terms = c("yearsToVisit", "rs4987358_T"),
                title = "SELL",
                axis.title = c("Years Since Plasma Sample","Cognition (DRS)"),
                axis.labels = NULL, legend.title = "rs4987358",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set1", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.18,0.20), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 18),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 16, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
print(p1)

#ICAM.2: combine C/C and C/T
MR_DRS_data$rs635634_T <- factor(MR_DRS_data$rs635634_T, 
                                  levels=c(0,1,2), 
                                  labels = c("C/C", "C/T + T/T", "C/T + T/T"))
MRdata$rs635634_T <- MR_DRS_data$rs635634_T
MR_DRS <- lmer(DRSTotalAge~rs635634_T*yearsToVisit +
                 AgeatTest + Sex + DRS.initial + (1 | INDDID), data = MR_DRS_data, REML = F) #no disease duration 
p1 = plot_model(MR_DRS, type = "pred", terms = c("yearsToVisit", "rs635634_T"),
                title = "ICAM.2",
                axis.title = c("Years Since Plasma Sample","Cognition (DRS)"),
                axis.labels = NULL, legend.title = "rs635634",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set1", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.18,0.20), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 18),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 16, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
print(p1)

#FCRL5: combine A/G and G/G
MR_DRS_data$rs914724_G <- factor(MR_DRS_data$rs914724_G, 
                                 levels=c(0,1,2), 
                                 labels = c("A/A", "A/G + G/G", "A/G + G/G"))
MRdata$rs914724_G <- MR_DRS_data$rs914724_G
MR_DRS <- lmer(DRSTotalAge~rs914724_G*yearsToVisit +
                 AgeatTest + Sex + DRS.initial + (1 | INDDID), data = MR_DRS_data, REML = F) #no disease duration 
p1 = plot_model(MR_DRS, type = "pred", terms = c("yearsToVisit", "rs914724_G"),
                title = "FCRL5",
                axis.title = c("Years Since Plasma Sample","Cognition (DRS)"),
                axis.labels = NULL, legend.title = "rs914724",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set1", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.18,0.20), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 18),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 16, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
print(p1)

#VWC: combine C/C and C/T
MR_DRS_data$rs740083_T <- factor(MR_DRS_data$rs740083_T, 
                                 levels=c(0,1,2), 
                                 labels = c("C/C + C/T","C/C + C/T", "T/T"))
MRdata$rs740083_T <- MR_DRS_data$rs740083_T
MR_DRS <- lmer(DRSTotalAge~rs740083_T*yearsToVisit +
                 AgeatTest + Sex + DRS.initial + (1 | INDDID), data = MR_DRS_data, REML = F) #no disease duration 
p1 = plot_model(MR_DRS, type = "pred", terms = c("yearsToVisit", "rs740083_T"),
                title = "VWC",
                axis.title = c("Years Since Plasma Sample","Cognition (DRS)"),
                axis.labels = NULL, legend.title = "rs740083",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set1", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.18,0.20), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 18),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 16, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
print(p1)

#Survival curves and cox models
MR_survival <- merge(survival.dat, MRdata, by = "INDDID")
MR_survival_obj <- Surv(MR_survival$YearsSinceSample, MR_survival$ClinicalConversion)
MR_survival$rs740083_T <- factor(MR_survival$rs740083_T, 
                            levels=c(0,1,2), 
                            labels = c("C/C + C/T","C/C + C/T", "T/T"))
MR.cox <- coxph(MR_survival_obj ~ AGE..YRS + Sex + rs740083_T, data = MR_survival)

ggforest(MR.cox, data=MR_survival, main = "Hazard Ratios", cpositions = c(0.02,0.23,0.4), 
         fontsize = 1, refLabel="reference", noDigits=2)


MR_cox <- function(protein, snp, basei, basef) {
  currcolname <- paste(snp,basef,sep="_")
  MR_survival[currcolname]<- factor(MR_survival[[currcolname]], 
                                    levels=c(0,1,2), 
                                    labels = c(paste(basei,basei,sep="/"), #A/A
                                               paste(basei,basef,sep="/"), #A/G
                                               paste(basef,basef,sep="/"))) #G/G
  myformula <- as.formula(paste("MR_survival_obj ~ AGE..YRS + Sex + ", currcolname, sep=""))
  MR.cox <- coxph(myformula, data = MR_survival)
  p1 <- ggforest(MR.cox, data=MR_survival, main = "Hazard Ratios", 
           fontsize = 1, refLabel="reference", noDigits=2)
  print(p1)
  
  survivalcurve <- ggadjustedcurves(MR.cox, palette = c('#e4191c','#4485bc', '#669999'),
                                    data=MR_survival,
                                    method="average",
                                    variable = currcolname,
                                    legend.title = snp)+
    labs(title=paste(protein, snp, sep=":"), x = "Time (Years)", y = "% Cognitively Normal") +
    theme_classic()
  print(survivalcurve)
  return(invisible())
}

proteinlist <- c("Dkk.4","Dkk.4","PDGFC","PDGFC","SELL","PVALB","PPIB",
                 "ICAM.2","KYNU","SCF","SCF","FCRL5","PTN.1","VWC2")
snplist <- c("rs10677223","rs12678359","rs55947119","rs17035367","rs4987358","rs4821544",
             "rs73453022","rs635634","rs3768844","rs6065904","rs3782181","rs914724",
             "rs1431093","rs740083")
baseilist <- c("A","T","C","A","G","T","T","C","A","G","A","A","C","C")
baseflist <- c("G","C","T","G","T","C","C","T","G","A","C","G","A","T")

MR_survival <- merge(survival.dat, MRdata, by = "INDDID")
MR_survival_obj <- Surv(MR_survival$YearsSinceSample, MR_survival$ClinicalConversion)
for(i in 3:length(proteinlist)){
  cox <- MR_cox(protein = proteinlist[i], 
                            snp = snplist[i], 
                            basei = baseilist[i], 
                            basef = baseflist[i])
}

i = 14
protein = proteinlist[i]
snp = snplist[i]
basei = baseilist[i]
basef = baseflist[i]
currcolname <- paste(snp,basef,sep="_")
MR_survival[currcolname]<- factor(MR_survival[[currcolname]], 
                                  levels=c(0,1,2), 
                                  labels = c(paste(basei,basei,sep="/"), #A/A
                                             paste(basei,basef,sep="/"), #A/G
                                             paste(basef,basef,sep="/"))) #G/G
myformula <- as.formula(paste("MR_survival_obj ~ AGE..YRS + Sex + ", currcolname, sep=""))
MR.cox <- coxph(myformula, data = MR_survival)
p1 <- ggforest(MR.cox, data=MR_survival, main = "Hazard Ratios", 
               fontsize = 1, refLabel="reference", noDigits=2)
print(p1)

survivalcurve <- ggadjustedcurves(MR.cox, palette = c('#e4191c','#4485bc', '#669999'),
                                  data=MR_survival,
                                  method="average",
                                  variable = currcolname,
                                  legend.title = snp)+
  labs(title=paste(protein, snp, sep=":"), x = "Time (Years)", y = "% Cognitively Normal") +
  theme_classic()
print(survivalcurve)
#################################################################################
# (9) COMPARE PROTEIN LEVELS ACROSS GROUPS (NC/AD/PD)
#################################################################################
olink.dat <- olink.dat[!is.na(olink.dat$AGE..YRS),]

#1) PDGFC 
ggplot(olink.dat, aes(x = Group, y = PDGFC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (PDGFC)", 
       y = "PDGFC Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(PDGFC ~ Group, olink.dat))
# Response: PDGFC
#Df  Sum Sq Mean Sq F value Pr(>F)
#Group       2  0.3014 0.15069  1.0132 0.3652
#Residuals 169 25.1332 0.14872  

#2) GALNT10
ggplot(olink.dat, aes(x = Group, y = GALNT10)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (GALNT10)", 
       y = "GALNT10 Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(GALNT10 ~ Group, olink.dat))
# Response: PDGFC
#Df  Sum Sq Mean Sq F value Pr(>F)
#Group       2  0.3014 0.15069  1.0132 0.3652
#Residuals 169 25.1332 0.14872  

#3) SELL
ggplot(olink.dat, aes(x = Group, y = SELL)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (SELL)", 
       y = "SELL Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(SELL ~ Group, olink.dat))
#Response: SELL
#Df  Sum Sq Mean Sq F value Pr(>F)
#Group       2  0.4695 0.23473  2.0386 0.1334
#Residuals 169 19.4584 0.11514   

#4) PVALB
ggplot(olink.dat, aes(x = Group, y = PVALB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (PVALB)", 
       y = "PVALB Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(PVALB ~ Group, olink.dat))
#Response: PVALB
#Df Sum Sq Mean Sq F value   Pr(>F)   
#Group       2  27.21 13.6069  5.3937 0.005361 **
#  Residuals 169 426.34  2.5227                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#5) ICAM.2
ggplot(olink.dat, aes(x = Group, y = ICAM.2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (ICAM.2)", 
       y = "ICAM.2 Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(ICAM.2 ~ Group, olink.dat))
#Response: ICAM.2
#Df  Sum Sq  Mean Sq F value Pr(>F)
#Group       2  0.1369 0.068436  0.5756 0.5635
#Residuals 169 20.0922 0.118889  

#6) KYNU
ggplot(olink.dat, aes(x = Group, y = KYNU)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (KYNU)", 
       y = "KYNU Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(KYNU ~ Group, olink.dat))
#Response: KYNU
#Df Sum Sq Mean Sq F value   Pr(>F)   
#Group       2  3.183 1.59130  6.7816 0.001469 **
#  Residuals 169 39.656 0.23465                    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#7) SCF
ggplot(olink.dat, aes(x = Group, y = SCF)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (SCF)", 
       y = "SCF Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(SCF ~ Group, olink.dat))
#Response: SCF
#Df  Sum Sq Mean Sq F value Pr(>F)
#Group       2  0.2345 0.11724  1.0783 0.3425
#Residuals 168 18.2672 0.10873 

#8) FCRL5
ggplot(olink.dat, aes(x = Group, y = FCRL5)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (FCRL5)", 
       y = "FCRL5 Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(FCRL5 ~ Group, olink.dat))
#Response: FCRL5
#Df Sum Sq Mean Sq F value Pr(>F)
#Group       2  0.863 0.43169  1.5113 0.2236
#Residuals 169 48.274 0.28564   

#9) C2
ggplot(olink.dat, aes(x = Group, y = C2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (C2)", 
       y = "C2 Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(C2 ~ Group, olink.dat))
#Response: C2
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Group       2  2.9555 1.47773  9.1597 0.0001672 ***
#  Residuals 169 27.2645 0.16133                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#10) FLT1
ggplot(olink.dat, aes(x = Group, y = FLT1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (FLT1)", 
       y = "FLT1 Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(FLT1 ~ Group, olink.dat))

#11) VWC2
ggplot(olink.dat, aes(x = Group, y = VWC2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Model Protein Levels By Disease Category (VWC2)", 
       y = "VWC2 Levels", x = "Disease Group") +
  theme(legend.position = "right")
# Stats
anova(lm(VWC2 ~ Group, olink.dat))


"declineRate~PDGFC+GALNT10+SELL+PVALB+ICAM.2+KYNU+SCF+FCRL5+C2+FLT1+VWC2"

#################################################################################
# (10) CHECK SOMASCAN DATASET
#################################################################################
setwd("/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/")
somalogic.dat <- read.csv("Soma compiled data with PC1,2 5-28-20.csv") #somalogic data set
somalogic_age.dat <- read.csv("Soma_Penn_AgeWithDecimal.csv")
somalogic_drs.dat <- read.csv("somalogic_DRS.csv")

somalogic_pd.dat <- somalogic.dat[somalogic.dat$Group == "PD",]
somalogic_pd.dat <- somalogic_pd.dat[somalogic_pd.dat$Site == "UPenn",] #96 somalogic PD Upenn patients
names(somalogic_pd.dat)[1] <- "INDDID"
soma_pd.dat <- merge(somalogic_pd.dat,somalogic_age.dat, by = "INDDID")
soma_pd.dat <- soma_pd.dat %>% relocate(SampleDate, .before = PlateId)
soma_pd.dat <- soma_pd.dat %>% relocate(DOB, .before = Sex)
soma_pd.dat <- soma_pd.dat %>% relocate(Age.y, .after = SampleDate)
names(soma_pd.dat)[8] <- "AGE..YRS"

#remove patients demented at baseline or DRS < 5
soma_pd.dat <- soma_pd.dat[-which(soma_pd.dat$INDDID == 101346),]
soma_pd.dat <- soma_pd.dat[-which(soma_pd.dat$INDDID == 107146),]
soma_pd.dat <- soma_pd.dat[-which(soma_pd.dat$INDDID == 116487),]
soma_pd.dat <- soma_pd.dat[-which(soma_pd.dat$INDDID == 107193),]
soma_pd.dat <- soma_pd.dat[-which(soma_pd.dat$INDDID == 111864),]

pdnum <- nrow(soma_pd.dat)
soma_data1 <- data.frame(INDDID = numeric(pdnum),
                    DRS.initial = numeric(pdnum),
                    date.initial = character(pdnum),
                    DRS.final = numeric(pdnum),
                    date.final = character(pdnum),
                    yearDiff = numeric(pdnum),
                    DRSnum = numeric(pdnum),
                    declineRate = numeric(pdnum),
                    declineSpeed = character(pdnum),
                    fastProgression = numeric(pdnum))
path_out2 <- "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/figures/"

proteinDrawDiff <- list()
patientstoremove <- c()
DRSrowstoremove <- c()
#for-loop to remove DRS taken before sample, compute DRS slope, and assign fast vs. slow
for (i in 1:pdnum) {
  soma_data1$INDDID[i] <- soma_pd.dat$INDDID[i]
  
  initialIndex <- match(soma_pd.dat$INDDID[i], somalogic_drs.dat$INDDID)
  
  allIndex <- which(somalogic_drs.dat$INDDID %in% soma_pd.dat$INDDID[i])
  allDRS <- somalogic_drs.dat$DRSTotalAge[allIndex]
  allYears <- (as.Date(somalogic_drs.dat$TestDate[allIndex], "%m/%d/%Y"))
  allYears <- decimal_date(allYears)
  
  #when olink sample was taken 
  sampleDate <- soma_pd.dat$SampleDate[i]
  sampleDate <- decimal_date(as.Date(sampleDate, "%m/%d/%Y"))
  
  #remove DRS taken before sample
  dateDiff <- allYears - sampleDate
  toRemoveIndex <- which(dateDiff < 0)
  if(length(toRemoveIndex) != 0) {
    allYears <- allYears[-toRemoveIndex]
    allDRS <- allDRS[-toRemoveIndex]
    DRSrowstoremove <- c(DRSrowstoremove, allIndex[toRemoveIndex])
  }
  soma_data1$DRSnum[i] <- length(allDRS)
  
  if(length(allDRS)<2){
    soma_data1$declineRate[i] <- NaN 
    proteinDrawDiff[i] <- min(abs(dateDiff))
    next 
  }
  
  soma_data1$DRS.initial[i] <- allDRS[1]
  soma_data1$date.initial[i] <- format(date_decimal(allYears[1]),"%m/%d/%Y")
  soma_data1$DRS.final[i] <- allDRS[length(allDRS)]
  soma_data1$date.final[i] <- format(date_decimal(allYears[length(allYears)]),"%m/%d/%Y")
  
  proteinDrawDiff[i] <- min(abs(dateDiff))
  if(proteinDrawDiff[i] > 0.5){
    patientstoremove <- c(patientstoremove,soma_data1$INDDID[i])
  }
  
  # find best fitting line
  bestLine <- lm(allDRS ~ allYears) 
  
  #set decline rate to slope of best fitting line
  slope <- coef(bestLine)[2]
  soma_data1$declineRate[i] <- coef(bestLine)[2] 
  
  
  jpeg(file=paste(path_out2, soma_data1$INDDID[i],"bestFitLineDRS.jpeg"))
  plot(allYears, allDRS,
       main = paste("INDDID", soma_data1$INDDID[i], ": Slope", slope),
       ylab = "DRS Score", xlab = "Date")
  abline(bestLine)
  dev.off()  
  
  #assign fast vs. slow
  if (soma_data1$declineRate[i] < -0.5) {
    soma_data1$declineSpeed[i] <- "Fast"
    soma_data1$fastProgression[i] <- 1
  }
  else {
    soma_data1$declineSpeed[i] <- "Slow"
    soma_data1$fastProgression[i] <- 0
  }
}

#remove DRS scores before sample
somalogic_drs.dat <- somalogic_drs.dat[-DRSrowstoremove,]

#remove patients with > 6 month difference between DRS and protein Draw
soma_data1 <- soma_data1[-match(patientstoremove,soma_data1$INDDID),]
soma_data1 <- soma_data1[-which(soma_data1$DRSnum == 1),]

#look at when DRS was taken relative to protein draw
hist(as.numeric(proteinDrawDiff),
     main = "Absolute Difference Between DRS and Protein Draw (Years) ")
summary(as.numeric(proteinDrawDiff))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.07178 0.09016 0.45205 

soma_data1 <- soma_data1 %>% arrange(declineRate)
slowdata <- subset(soma_data1, declineSpeed == "Slow")
nrow(slowdata) #31 slow progressing patients 
fastdata <- subset(soma_data1, declineSpeed == "Fast")
nrow(fastdata) #27 fast progressing patients 

#IDENTIFY SLOPE-ASSOCIATED PROTEINS
protein_num <- ncol(soma_pd.dat) - 9 - 1 #olink panel size

soma_pd.dat <- subset(soma_pd.dat, select=-SomaVersion)
soma_data3 <- merge(soma_data1, soma_pd.dat)
soma_data3$Sex <- factor(soma_data3$Sex, levels=c("Female","Male"), labels=c(0,1))

soma_slopemodel.df <- data.frame(Protein = character(protein_num),
                            Protein.coef = numeric(protein_num),
                            Protein.pval = numeric(protein_num),
                            Protein.FDR = numeric(protein_num),
                            Male.coef = numeric(protein_num),
                            Male.pval = numeric(protein_num),
                            Male.FDR = numeric(protein_num),
                            Age.coef = numeric(protein_num),
                            Age.pval = numeric(protein_num),
                            Age.FDR = numeric(protein_num))

# LINEAR REGRESSION 
# DRS slope as a function of protein + sex + age 
for (i in 1:protein_num) {
  protein_name <- names(soma_pd.dat)[i+9]
  soma_slopemodel.df$Protein[i] <- protein_name
  formula_i <- paste("-declineRate ~", protein_name, "+Sex", "+AGE..YRS", sep = "")
  model_i <- lm(formula_i, soma_data3)
  soma_slopemodel.df$Protein.coef[i] <- summary(model_i)$coefficients[2,1]
  soma_slopemodel.df$Protein.pval[i] <- summary(model_i)$coefficients[2,4]
  soma_slopemodel.df$Male.coef[i] <- summary(model_i)$coefficients[3,1]
  soma_slopemodel.df$Male.pval[i] <- summary(model_i)$coefficients[3,4]
  soma_slopemodel.df$Age.coef[i] <- summary(model_i)$coefficients[4,1]
  soma_slopemodel.df$Age.pval[i] <- summary(model_i)$coefficients[4,4]
}

# Benjamini-Hochberg Procedure (p-value adjustment)
soma_slopemodel.df$Protein.FDR <- p.adjust(soma_slopemodel.df$Protein.pval, method = 'BH')
soma_slopemodel.df$Male.FDR <- p.adjust(soma_slopemodel.df$Male.pval, method = 'BH')
soma_slopemodel.df$Age.FDR <- p.adjust(soma_slopemodel.df$Age.pval, method = 'BH')

sum(soma_slopemodel.df$Protein.pval < 0.1) #48
sum(soma_slopemodel.df$Protein.pval < 0.05) #16 slope-associated proteins
sum(soma_slopemodel.df$Protein.pval < 0.01) #2

soma_slopemodel.df <- soma_slopemodel.df %>% arrange(Protein.pval)
path_out = "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/output/"
fileName = paste(path_out, '20220605_soma_linearregression.csv', sep = '')
write.csv(soma_slopemodel.df,
          file = fileName,
          row.names=FALSE)

# Predict
soma_model <- lm("declineRate~Dkk.4+sL.Selectin+sICAM.2+KYNU+PTN+C2+AGE..YRS+Sex", soma_data3)
soma_data3$SlopePred <- predict(soma_model, soma_data3)
soma_data3$DeltaSlope <- soma_data3$SlopePred - soma_data3$declineRate

# Stats
min(soma_data3$declineRate) #-2.756962
max(soma_data3$declineRate) #0.9374868
min(soma_data3$SlopePred) # -1.76853
max(soma_data3$SlopePred) #0.07864312

cor(soma_data3$declineRate,soma_data3$SlopePred, method = "pearson")
#0.5391972

dev.off()
par(mfrow=c(1,1))
plot(x = soma_data3$declineRate,
     y = soma_data3$SlopePred,
     xlim = c(-3,1), ylim = c(-3,1),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = 'SOMALOGIC Proteomic vs. Actual DRS Slope')
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -2, labels = "cor = 0.539")


#OLINK COMPARISON
# Predict
olinksoma_model <- lm("declineRate~Dkk.4+SELL+ICAM.2+KYNU+PTN.1+C2+AGE..YRS+Sex", data3)
olinksoma_slopepred <- predict(olinksoma_model, data3)
olinksoma_DeltaSlope <- olinksoma_slopepred - data3$declineRate

# Stats
min(data3$declineRate) # -2.363552
max(data3$declineRate) # 1.792649
min(olinksoma_slopepred) #-2.203426
max(olinksoma_slopepred) #1.730783

cor(data3$declineRate,olinksoma_slopepred, method = "pearson")
#0.7685652

dev.off()
par(mfrow=c(1,1))
plot(x = data3$declineRate,
     y = olinksoma_slopepred,
     xlim = c(-2.5,2), ylim = c(-2.5,2),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Proteomic DRS Slope',
     main = 'Olink Subset Proteomic vs. Actual DRS Slope')
abline(a = 0, b = 1, col = 'blue')
text(x = 0.5, y = -2, labels = "cor = 0.769")

soma_model$coefficients[2:9]
#Dkk.4 sL.Selectin     sICAM.2        KYNU      SCF.sR         PTN          C2    AGE..YRS        Sex1 
#0.70971864 -1.55059793  0.99519478  1.00940047 -0.53101939 -0.96144925  1.72930880 -0.03328338  0.06343483 
barplot(soma_model$coefficients[2:9])
olinksoma_model$coefficients[2:9]
#Dkk.4         SELL       ICAM.2         KYNU          SCF        PTN.1           C2     AGE..YRS         Sex1 
#-0.536115555 -0.679534171  0.857069025  0.409882559 -0.590767131 -0.105246806 -0.331315998  0.004895992  0.253078649  
barplot(olinksoma_model$coefficients[2:9])


#################################################################################
# (11) OLINK REPLICATION 
#################################################################################
#get data
setwd("/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/1_olink_replication/replication/")
#olink_rep.dat <- read.csv("20220726_olinkreplication_log.csv") #log transformed data
olink_rep.dat <- read.csv("20220726_olinkreplication.csv") #non-log transformed
olink_rep_pd.dat <- olink_rep.dat[olink_rep.dat$Group == "PD",]
olink_rep_pd.dat <- olink_rep_pd.dat[!is.na(olink_rep_pd.dat$INDDID),]#remove NA inddid
drs_rep <- read.csv("3_drs.csv") 

## FIND DRS SLOPE
# Extract DRS data for PD patients 
pd_drs_rep.dat <- drs_rep %>% filter(drs_rep$INDDID %in% olink_rep_pd.dat$INDDID)

pdrepnum <- nrow(olink_rep_pd.dat) # number of PD patients (100)
data1_rep <- data.frame(INDDID = numeric(pdrepnum),
                    DRS.initial = numeric(pdrepnum),
                    date.initial = character(pdrepnum),
                    DRS.final = numeric(pdrepnum),
                    date.final = character(pdrepnum),
                    yearDiff = numeric(pdrepnum),
                    DRSnum = numeric(pdrepnum),
                    declineRate = numeric(pdrepnum),
                    declineSpeed = character(pdrepnum),
                    fastProgression = numeric(pdrepnum))
path_out2 <- "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/1_olink_replication/replication/output/"

proteinDrawDiff <- list()
DRSrowstoremove <- c()

#for-loop to remove DRS taken before sample, compute DRS slope, and assign fast vs. slow
for (i in 1:pdrepnum) {
  data1_rep$INDDID[i] <- olink_rep_pd.dat$INDDID[i]
  
  initialIndex <- match(olink_rep_pd.dat$INDDID[i], pd_drs_rep.dat$INDDID)
  
  allIndex <- which(pd_drs_rep.dat$INDDID %in% olink_rep_pd.dat$INDDID[i])
  allDRS <- pd_drs_rep.dat$DRSTotalAge[allIndex]
  allYears <- (as.Date(pd_drs_rep.dat$TestDate[allIndex], "%m/%d/%Y"))
  allYears <- decimal_date(allYears)
  
  #when olink sample was taken 
  sampleDate <- olink_rep_pd.dat$SampleDate[i]
  sampleDate <- decimal_date(as.Date(sampleDate, "%m/%d/%y"))
  
  #remove DRS taken before sample
  dateDiff <- allYears - sampleDate
  toRemoveIndex <- which(dateDiff < 0)
  if(length(toRemoveIndex) != 0) {
    allYears <- allYears[-toRemoveIndex]
    allDRS <- allDRS[-toRemoveIndex]
    DRSrowstoremove <- c(DRSrowstoremove, allIndex[toRemoveIndex])
  }
  data1_rep$DRSnum[i] <- length(allDRS)
  
  if(length(allDRS)<1){
    data1_rep$declineRate[i] <- NaN 
    next 
  }
  
  data1_rep$DRS.initial[i] <- allDRS[1]
  data1_rep$date.initial[i] <- format(date_decimal(allYears[1]),"%m/%d/%Y")
  data1_rep$DRS.final[i] <- allDRS[length(allDRS)]
  data1_rep$date.final[i] <- format(date_decimal(allYears[length(allYears)]),"%m/%d/%Y")
  data1_rep$yearDiff[i] <- allYears[length(allYears)]-allYears[1]
  
  proteinDrawDiff[i] <- min(abs(dateDiff))
  
  # find best fitting line
  bestLine <- lm(allDRS ~ allYears) 
  
  #set decline rate to slope of best fitting line
  slope <- coef(bestLine)[2]
  data1_rep$declineRate[i] <- coef(bestLine)[2] 
  data1_rep$bestFitIntercept[i] <- summary(bestLine)$coefficients[1,1]
  
  
  jpeg(file=paste(path_out2, data1_rep$INDDID[i],"bestFitLineDRS.jpeg"))
  plot(allYears, allDRS,
       main = paste("INDDID", data1_rep$INDDID[i], ":", slope),
       ylab = "DRS", xlab = "Year", cex.lab = 1.5,
       pch = 21, bg = "black")
  abline(bestLine, col = "blue", lwd = 3, lty = 2)
  dev.off()  
  
  #assign fast vs. slow
  if (data1_rep$declineRate[i] < -0.5) {
    data1_rep$declineSpeed[i] <- "Fast"
    data1_rep$fastProgression[i] <- 1
  }
  else {
    data1_rep$declineSpeed[i] <- "Slow"
    data1_rep$fastProgression[i] <- 0
  }
}

#remove DRS scores before sample
pd_drs_rep.dat <- pd_drs_rep.dat[-DRSrowstoremove,]

#look at when DRS was taken relative to protein draw
hist(as.numeric(proteinDrawDiff),
     main = "Absolute Difference Between DRS and Protein Draw (Years) ")
summary(as.numeric(proteinDrawDiff))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.06608 0.01982 0.49863 

data1_rep <- data1_rep %>% arrange(declineRate)
slowdata <- subset(data1_rep, declineSpeed == "Slow")
nrow(slowdata) #56 slow progressing patients 
fastdata <- subset(data1_rep, declineSpeed == "Fast")
nrow(fastdata) #44 fast progressing patients 

##MODEL 
protein_num_rep <- ncol(olink_rep_pd.dat) - 10 #olink panel size #552
data3_rep <- cbind(data1_rep, olink_rep_pd.dat)
data3_rep$Sex <- factor(data3_rep$Sex, levels=c("Female","Male"), labels=c(0,1))
data3_rep <- data3_rep[-which(data3_rep$TubeNumber == 122985),] #remove samples with quality control warnings
data3_rep <- data3_rep[-which(data3_rep$TubeNumber == 121670),] 
data3_rep <- data3_rep[-which(data3_rep$TubeNumber == 120834),] 
#data3_rep <- data3_rep[-which(data3_rep$INDDID == 104767),] #remove outlier

# Predict
final_model_rep <- final_model
final_model_rep <- lm("declineRate~Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2", data3_rep)
final_model_rep <- lm("declineRate~CRIM1+ GM.CSF.R.alpha+BANK1+TNC+CD69+MMP.3+LAT2+COMP+INPPL1+FSTL3",data3_rep)
#data3_rep$SlopePred <- predict(final_model, data3_rep) #discovery model coefficients
data3_rep$SlopePred <- predict(final_model_rep, data3_rep) #fit
data3_rep$DeltaSlope <- data3_rep$SlopePred - data3_rep$declineRate

# Stats
min(data3_rep$declineRate) #-3.805227
max(data3_rep$declineRate) #1.460725
min(data3_rep$SlopePred) #-1.102824 #-1.679344
max(data3_rep$SlopePred) #0.2089023 #2.217465

cor(data3_rep$declineRate,data3_rep$SlopePred, method = "pearson", use="complete.obs")
#not log 0.410872, 0.3777222
#log 0.3801344 
#log, removed qc warning samples 0.3804662

par(mfrow=c(1,1))
plot(x = data3_rep$declineRate,
     y = data3_rep$SlopePred,
     xlim = c(-4,3), ylim = c(-4,3),
     type = 'p', col = 'blue', pch = 16,
     xlab = 'Actual DRS Slope',
     ylab = 'Predicted DRS Slope',
     cex.lab = 1.5,
     main = 'Replication Cohort')
abline(a = 0, b = 1, col = 'blue')
text(x = 0, y = -3, labels = "cor = 0.58")

barplot(final_model_rep$coefficients[2:(13)])

colnames(data3_rep) <- make.unique(names(data3_rep))
ggplot(data3_rep, aes(x = Group, y = DeltaSlope)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Group)) +
  labs(title = "Delta DRS Slope\n(11-parameter model)", 
       y = "Delta Slope", x = "Disease Group") +
  theme(legend.position = "right")

#Binary Classifier for full data set 
# Let's train the model using the 12 proteins
#data3_rep <- data3_rep[-which(is.na(data3_rep$SlopePred)),] #remove samples without coverage
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
model_all_rep <- train(as.factor(fastProgression) ~
                        AGE..YRS +
                        Sex + 
                         Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
                      method = "glm",
                      data = data3_rep,
                      trControl = train_control,
                      na.action=na.pass,
                      family = "binomial")
summary(model_all_rep)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control_all_rep <- train(as.factor(fastProgression) ~
                                 AGE..YRS+
                                Sex,
                              method = "glm",  
                              data = data3_rep,
                              trControl = train_control,
                              na.action=na.pass,
                              family = "binomial")

# Predict using test set 
predict_rep <- predict(model_all_rep, newdata = data3, type= "prob")
predictControl_rep <- predict(model.control_all_rep, newdata = data3, type= "prob")

roc_rep <- roc(data3$fastProgression, as.numeric(predict_rep$'1'))
rocControl_rep <- roc(data3$fastProgression, as.numeric(predictControl_rep$'1'))

#95% confidence intervals for AUC and accuracy
ci.auc(roc_rep) #0.621-0.8348 #0.6199-0.8322 
predict_rep <- predict(model_all_rep, newdata = data3_rep)
xtab <- table(predict_rep, data3_rep$fastProgression)
confusionMatrix(xtab) #accuracy 0.6939 (0.5926, 0.783) #accuracy 0.7216 (0.6214, 0.8079)

ci.auc(rocControl_rep) #0.5381-0.757 #0.5403-0.7602 (DeLong)
predictControl_rep <- predict(model.control_all_rep, newdata = data3_rep)
xtabcontrol <- table(predictControl_rep, data3_rep$fastProgression)
confusionMatrix(xtabcontrol) #accuracy 0.6224 (0.5188, 0.7184) #0.6289 ((0.5248, 0.7248))  

par(mar = c(5, 5, 5, 5))
plot(roc_rep, col = "blue", legacy.axes = TRUE, print.auc=TRUE, print.auc.y=.5)
plot(rocControl_rep, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .4)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "Test", line = 3.5)

data2_rep <- merge(pd_drs_rep.dat, data1_rep, by = "INDDID")

# Calculate # of years since first visit for each test
drsnum <- nrow(data2_rep)
for (i in 1:drsnum){
  data2_rep$yearsToVisit[i] <- time_length(difftime(as.Date(data2_rep$TestDate[i], "%m/%d/%Y"), 
                                                as.Date(data2_rep$date.initial[i], "%m/%d/%Y")),
                                       "years")
}

# Add olink proteins
data2_rep <- merge(data2_rep, olink_rep_pd.dat, by = "INDDID")

# Code sex from female male into 0 (Female) and 1 (Male)
data2_rep$Sex <- factor(data2_rep$Sex, levels=c("Female","Male"), labels=c(0,1))
data2_rep$fastProgression <- factor(data2_rep$fastProgression, levels = c(0, 1), labels=c("Slow","Fast"))

# Train LME model to see how fast and slow progressors differ
DRS <- lmer(DRSTotalAge ~ fastProgression*yearsToVisit + Sex +
              (1 | INDDID), data=data2_rep, REML=F) #no disease duration 

data(efc)
theme_set(theme_sjplot())
p1 = plot_model(DRS, type = "pred", terms = c("yearsToVisit", "fastProgression"),
                title = "Cognitive Change",
                axis.title = c("Years Since Plasma Sample","DRS"),
                axis.labels = NULL, legend.title = "Progression Group",
                wrap.title = 50, wrap.labels = 25, axis.lim = NULL, grid.breaks = NULL,
                ci.lvl = 0.95, se = NULL, colors = "Set2", order.terms = c(3,1,2)) + theme_classic() +
  theme(legend.position = c(0.25,0.18), legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 15,),
        title = element_text(size = 20),plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text = element_text(size = 20, color = "black"), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1))
p1

#IDENTIFY SLOPE-ASSOCIATED PROTEINS
protein_num_rep <- ncol(olink_rep_pd.dat) - 10 #552

slopemodel_rep <- data.frame(Protein = character(protein_num_rep),
                            Protein.coef = numeric(protein_num_rep),
                            Protein.pval = numeric(protein_num_rep),
                            Protein.FDR = numeric(protein_num_rep),
                            Male.coef = numeric(protein_num_rep),
                            Male.pval = numeric(protein_num_rep),
                            Male.FDR = numeric(protein_num_rep),
                            Age.coef = numeric(protein_num_rep),
                            Age.pval = numeric(protein_num_rep),
                            Age.FDR = numeric(protein_num_rep))

# LINEAR REGRESSION 
# DRS slope as a function of protein + sex + age 
for (i in 1:protein_num_rep) {
  protein_name <- names(olink_rep_pd.dat)[i+10]
  if (sum(!is.na(data3_rep[[protein_name]])) < 4) {
    next
  }
  slopemodel_rep$Protein[i] <- protein_name
  formula_i <- paste("-declineRate ~", protein_name, "+Sex", "+AGE..YRS", sep = "")
  model_i <- lm(formula_i, data3_rep)
  slopemodel_rep$Protein.coef[i] <- summary(model_i)$coefficients[2,1]
  slopemodel_rep$Protein.pval[i] <- summary(model_i)$coefficients[2,4]
  slopemodel_rep$Male.coef[i] <- summary(model_i)$coefficients[3,1]
  slopemodel_rep$Male.pval[i] <- summary(model_i)$coefficients[3,4]
  slopemodel_rep$Age.coef[i] <- summary(model_i)$coefficients[4,1]
  slopemodel_rep$Age.pval[i] <- summary(model_i)$coefficients[4,4]
  print(i)
}

# Benjamini-Hochberg Procedure (p-value adjustment)
slopemodel_rep$Protein.FDR <- p.adjust(slopemodel_rep$Protein.pval, method = 'BH')
slopemodel_rep$Male.FDR <- p.adjust(slopemodel_rep$Male.pval, method = 'BH')
slopemodel_rep$Age.FDR <- p.adjust(slopemodel_rep$Age.pval, method = 'BH')

# Summary
summary(slopemodel_rep$Protein.FDR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.8912  0.8912  0.9096  0.9399  0.9999 
summary(slopemodel_rep$Male.FDR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1358  0.1358  0.1475  0.1358  0.8157
summary(slopemodel_rep$Age.FDR)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.3643  0.3643  0.3730  0.3643  0.7847 

summary(slopemodel_rep$Protein.pval)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.2564  0.4507  0.4809  0.7056  0.9999 
summary(slopemodel_rep$Male.pval)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.07588 0.08696 0.10108 0.10394 0.81573 
summary(slopemodel_rep$Age.pval)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.2058  0.2265  0.2481  0.2759  0.7847 

sum(slopemodel_rep$Protein.pval < 0.1, na.rm=TRUE ) #57
sum(slopemodel_rep$Protein.pval < 0.05, na.rm = TRUE) #20
sum(slopemodel_rep$Protein.pval < 0.01, na.rm = TRUE) #5

slopemodel_rep <- slopemodel_rep %>% arrange(Protein.pval)
path_out = "/Users/kristenpark/Desktop/Chen-Plotkin/2_Cognitive Decline Project/final model/1_olink_replication/replication/"
fileName = paste(path_out, '20220727_replication_slope_linearregression.csv', sep = '')
write.csv(slopemodel_rep,
          file = fileName,
          row.names=FALSE)

#################################################################################
# (12) OLINK REPLICATION 
#################################################################################

#baseline DRS, years of follow up, age, sex -- add in disease duration? 
ogdata <- select(data3, INDDID, DRS.initial, yearDiff, AGE..YRS, Sex, declineRate, fastProgression,
                 Dkk.4, PDGFC, SELL, PVALB, PPIB, ICAM.2, KYNU, SCF, 
                 FCRL5, PTN.1, VWC2, C2)
repdata <- select(data3_rep, INDDID, DRS.initial, yearDiff, AGE..YRS, Sex, declineRate, fastProgression,
                  Dkk.4, PDGFC, SELL, PVALB, PPIB, ICAM.2, KYNU, SCF, 
                  FCRL5, PTN.1, VWC2, C2)
fulldata <- rbind(ogdata, repdata)

#Determine train/test divisions (75% train, 25% test)
set.seed(140)
trainIndex <- createDataPartition(fulldata$fastProgression,p=0.8,list=FALSE)
trainSet <- fulldata[trainIndex,] #training data (75% of data)
testSet <- fulldata[-trainIndex,] #testing data (25% of data)

# Let's train the model using the 12 proteins
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 50)
model <- train(as.factor(fastProgression) ~
                 DRS.initial+yearDiff+AGE..YRS+Sex+
                 Dkk.4+PDGFC+SELL+PVALB+PPIB+ICAM.2+KYNU+SCF+FCRL5+PTN.1+VWC2+C2,
               method = "glm",
               data = trainSet,
               trControl = train_control,
               na.action=na.pass,
               family = "binomial")
summary(model)

# Let's also train a model with clinical variables (Sex + Age only as a control)
model.control <- train(as.factor(fastProgression) ~
                         DRS.initial+yearDiff+AGE..YRS+Sex,
                       method = "glm",  
                       data = trainSet,
                       trControl = train_control,
                       na.action=na.pass,
                       family = "binomial")

# Predict using test set 
predict <- predict(model, newdata = trainSet, type= "prob")
predictControl <- predict(model.control, newdata = trainSet, type= "prob")

roc <- roc(testSet$fastProgression, as.numeric(predict$'1'))
rocControl <- roc(testSet$fastProgression, as.numeric(predictControl$'1'))

par(mar = c(5, 5, 5, 5))
plot(roc, col = "blue", legacy.axes = TRUE, print.auc=TRUE)
plot(rocControl, add=TRUE, col = "grey", legacy.axes = TRUE, print.auc=TRUE,print.auc.y = .3)
legend("bottomright", legend = c("Model", "Clinical Only"), col = c("blue", "grey"), lty =1)
title(main = "75% Train, 25% Test", line = 3.5)
# AUC: 0.814, control AUC: 0.543

cor(data3$Dkk.4, data3$declineRate) #-0.3910111
cor(data3_rep$Dkk.4, data3_rep$declineRate) #0.1059275

cor(data3$PDGFC, data3$declineRate) #-0.3626149
cor(data3_rep$PDGFC, data3_rep$declineRate) #0.03475833

cor(data3$SELL, data3$declineRate) #-0.3837619
cor(data3_rep$SELL, data3_rep$declineRate) #0.008800445

cor(data3$PVALB, data3$declineRate) #0.3721087
cor(data3_rep$PVALB, data3_rep$declineRate) #-0.2174983

cor(data3$PPIB, data3$declineRate) #0.3073909
cor(data3_rep$PPIB, data3_rep$declineRate) #-0.1510173

cor(data3$ICAM.2, data3$declineRate) # 0.3079582
cor(data3_rep$ICAM.2, data3_rep$declineRate) #-0.1358944

cor(data3$KYNU, data3$declineRate) # 0.3051549
cor(data3_rep$KYNU, data3_rep$declineRate) #-0.06829105

cor(data3$SCF, data3$declineRate) # -0.3111904
cor(data3_rep$SCF, data3_rep$declineRate) # -0.03187457'

cor(data3$FCRL5, data3$declineRate) # 0.3567889
cor(data3_rep$FCRL5, data3_rep$declineRate) # -0.02381553

cor(data3$PTN.1, data3$declineRate) #-0.307803
cor(data3_rep$PTN.1, data3_rep$declineRate) #-0.03542907

cor(data3$VWC2, data3$declineRate) #-0.3103235
cor(data3_rep$VWC2, data3_rep$declineRate) #0.1533796

cor(data3$C2, data3$declineRate) #-0.2982819
cor(data3_rep$C2, data3_rep$declineRate) #-0.01250206


