# SLE-Meta-analysis

##Load metafor and import data
library(metafor)
library(ggplot2)
library(bayesmeta)
library(forestplot)


data11 = read.csv(file.choose())


###########################
#synthesis and meta analysis

data1 = data11[-c(13),]

##Calculate proportions to double check
data1$pi = data1$Pos/data1$Scr

##Calculate Freeman-Tukey Transformation
data1 <- escalc(measure="PFT", xi=Pos, ni=Scr, data=data1)
data1$check = transf.ipft(data1$yi, data1$Scr)

##Run Random-effects model using REML estimator
Model <- rma(yi, vi, method="REML", data=data1)
Model

#Funnel plot asymmetry
regtest(Model)

Predictor <- predict(Model, transf=transf.ipft.hm, targs=list(ni=data1$Scr))

##Putting this here for safekeeping
mlabfun <- function(text, res) { list(bquote(paste(.(text),
                                                   " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                                                   ", df = ", .(res$k - res$p),
                                                   ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                                                   I^2, " = ", .(formatC(res$I2, digits=1, format="f")), "%, ",
                                                   tau^2, " = ", .(formatC(res$tau2, digits=2, format="f")), ")")))}

forest(Model, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"), refline = 0.20,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,4), order = data1$pi, showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("All Studies, Study 11 dichotomized",Model))

funnel(Model,atransf=transf.ipft.hm, targs=list(ni=data1$Scr),
       las=1, ylim = c(0,0.1))

funnel(Model,atransf=transf.ipft.hm, targs=list(ni=data1$Scr),
       las=1)

funnel(Model, targs=list(ni=data1$Scr),
       las=1)

funnel(trimfill(Model), las = 1, ylim = c(0,0.1))

##Outlier diagnostics
par(mar=c(5,6,4,2))
plot(influence(Model),cex=0.8,las=1)
influence(Model)
rstudent(Model)
Model.tf = trimfill(Model)
Model.tf
baujat(Model)


##Remove Outliers
data2 = data1[-c(5, 7, 12),]
ModelOutlier <- rma(yi, vi, method="REML", data=data2)
ModelOutlier
forest(ModelOutlier, transf=transf.ipft.hm, targs=list(ni=data2$Scr),
       slab = paste(data2$Author,", ",data2$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"), refline = 0.22,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black', 
       efac = c(3,4), order = data2$pi, showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("Outliers Excluded",ModelOutlier))


##Meta regression model for year of study
mModel = rma(yi,vi,method = 'REML', mods=~Year, data=data1)
mModel
regplot(mModel, xlab = 'Year of Study', transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 16)

regplot(mModel, xlab = 'Year of Study, Estimate: 0.002, p value: 0.86', transf = transf.ipft.hm, 
        targs = list(ni = data1$Scr), pch = 21, col = "black", bg = "white")

##Meta regression model for FUP time
mModelFUP = rma(yi,vi,method = 'REML', mods=~FUP, data=data1)
mModelFUP
regplot(mModelFUP, xlab = 'Follow up time(months), Estimate: -0.0001, p value: 0.96', transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 21, col = "black", bg = "white")

##Meta regression model for CS Dose

mModelCSDOSE = rma(yi,vi,method = 'REML', mods=~CsDose.mg.day., data=data1)
mModelCSDOSE
regplot(mModelCSDOSE, xlab = 'CS Daily Dose (mg), Study 11 dichotomized, Estimate: 0.005, p value: 0.0052', transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 21, col = "black", bg = "white")

##Meta regression model for year of study

mModelPulse = rma(yi,vi,method = 'REML', mods=~Pulse , data=data1)
mModelPulse
regplot(mModelPulse, xlab = '% of participant received CS Pulse therapy, Estimate: 0.0006, p value: 0.81 ', transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 21, col = "black", bg = "white")

##Meta regression model for age of the participants
mModelAge = rma(yi,vi,method = 'REML', mods=~Age , data=data1)
mModelAge
regplot(mModelAge, xlab = 'mean Age All Patients, Estimate: 0.01, p value: 0.37', transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 21, col = "black", bg = "white")


##Meta regression model for age of the adult paricipants
mModelAgeadult <- rma(yi, vi, method = 'REML', mods = ~Age, data = subset(data1, Age >= 18))
mModelAgeadult
regplot(mModelAgeadult, xlab = 'mean Age Adult Patients (Age>=18), Estimate: -0.024, p value: 0.41', transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 21, col = "black", bg = "white")


##Meta regression model for APL antibody
mModelApl = rma(yi,vi,method = 'REML', mods=~aPL , data=data1)
mModelApl
regplot(mModelApl, xlab = "Proportion of APL Ab positive participants All Patients, Estimate: -0.007, p value: 0.006", transf=transf.ipft.hm, targs=list(ni=data1$Scr),
        pch = 21, col = "black", bg = "white")


##Sensitivity Analysis
ModelCs = rma(yi, vi, method = 'REML', subset=(Csuse == 1), data=data1)
ModelCsN = rma(yi, vi, method = 'REML', subset=(Csuse == 0), data=data1)

ModelMale = rma(yi,vi, method = 'REML', subset = (GenderM > 9), data = data1)
ModelMaleN = rma(yi,vi, method = 'REML', subset = (GenderM < 10), data = data1)

ModelFollowup = rma(yi,vi, method = 'REML', subset = (FUP > 12), data = data1)
ModelFollowupN = rma(yi,vi, method = 'REML', subset = (FUP <13), data = data1)

ModelAge30Below <- rma(yi, vi, method = 'REML', subset = (Age <= 30), data = data1)
ModelAgeAbove30 <- rma(yi, vi, method = 'REML', subset = (Age > 30), data = data1)

ModelTesla0.5 <- rma(yi, vi, method = 'REML', subset = (tesla == 0.5), data = data1)
ModelTesla1.5 <- rma(yi, vi, method = 'REML', subset = (tesla == 1.5), data = data1)

ModelMSK0 <- rma(yi, vi, method = 'REML', subset = (MSK == 0), data = data1)
ModelMSK1 <- rma(yi, vi, method = 'REML', subset = (MSK == 1), data = data1)

Modelinterpreter1 <- rma(yi, vi, method = 'REML', subset = (interpreter == 1), data = data1)
Modelinterpreter2 <- rma(yi, vi, method = 'REML', subset = (interpreter == 2), data = data1)

Modelbias <- rma(yi, vi, method = 'REML', subset = (bias == 1), data = data1)


#Function to add statistics to forest plots
forest(ModelCs, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.25,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("Corticosteroid Use", ModelCs))

forest(ModelCsN, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.25,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("No Corticosteroid Use", ModelCsN))

forest(ModelMale, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.22,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun(">10% Male Participants", ModelMale))

forest(ModelMaleN, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.22,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("<10% Male Participants", ModelMaleN))

forest(ModelOutlier, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data2$Author,", ",data2$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.26,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("Influential Study Removed", ModelOutlier))


forest(ModelFollowup, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun(">12 Months Followup", ModelFollowup))

forest(ModelFollowupN, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("<12 Months Followup", ModelFollowupN))

forest(ModelAge30Below, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("30 year old and below", ModelAge30Below))


forest(ModelAgeAbove30, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("above 30 year old", ModelAgeAbove30))

forest(ModelMSK0, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("MSK 0", ModelMSK0))

forest(ModelMSK1, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("MSK 1", ModelMSK1))

forest(Modelinterpreter1, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("interpreter 1", Modelinterpreter1))

forest(Modelinterpreter2, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("interpreter 2", Modelinterpreter2))

forest(ModelTesla0.5, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("Tesla 0.5", ModelTesla0.5))

forest(ModelTesla1.5, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("Tesla 1.5", ModelTesla1.5))


forest(Modelbias, transf=transf.ipft.hm, targs=list(ni=data1$Scr),
       slab = paste(data1$Author,", ",data1$Year, sep=''), 
       header=c("Study", "Weights Proportion [95% CI]"),refline = 0.29,
       col = 'azure', pch = 0, plim = c(0,2), colout = 'black',
       font = 'Arial', efac = c(3,3), showweights = T,
       addfit=T, addpred=T, mlab=mlabfun("Tesla 1.5", Modelbias))

##Comparison of groups
Influential <- data.frame(estimate = c(coef(Model), coef(ModelOutlier)), stderror = c(Model$se, ModelOutlier$se),
                          meta = c("ALL","Exclusion"), tau2 = round(c(Model$tau2, ModelOutlier$tau2),3))
Influential

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=Influential, digits=3)

#Comparison of Corticosteroid Use
Cort <- data.frame(estimate = c(coef(ModelCs), coef(ModelCsN)), stderror = c(ModelCs$se, ModelCsN$se),
                   meta = c("CSUSE","NOUSE"), tau2 = round(c(ModelCs$tau2, ModelCsN$tau2),3))
Cort

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=Cort, digits=3)

#Comparison of Male Participants

Male <- data.frame(estimate = c(coef(ModelMale), coef(ModelMaleN)), stderror = c(ModelMale$se, ModelMaleN$se),
                   meta = c(">10","<10"), tau2 = round(c(ModelMale$tau2, ModelMaleN$tau2),3))
Male

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=Male, digits=3)

#Comparison of Follow Up Time

FollowUp <- data.frame(estimate = c(coef(ModelFollowup), coef(ModelFollowupN)), stderror = c(ModelFollowup$se, ModelFollowupN$se),
                       meta = c(">12","<12"), tau2 = round(c(ModelFollowup$tau2, ModelFollowupN$tau2),3))
FollowUp

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=FollowUp, digits=3)


#Comparison of MSK interpreter
MSK <- data.frame(estimate = c(coef(ModelMSK0), coef(ModelMSK1)), stderror = c(ModelMSK0$se, ModelMSK1$se),
                  meta = c("0","1"), tau2 = round(c(ModelMSK0$tau2, ModelMSK1$tau2),3))
MSK

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=MSK, digits=3)



#Comparison of number of interpreters
interpreter <- data.frame(estimate = c(coef(Modelinterpreter1), coef(Modelinterpreter2)), stderror = c(Modelinterpreter1$se, Modelinterpreter2$se),
                          meta = c("0","1"), tau2 = round(c(Modelinterpreter1$tau2, Modelinterpreter2$tau2),3))
interpreter

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=interpreter, digits=3)


#Comparison of Age
age <- data.frame(estimate = c(coef(ModelAge30Below), coef(ModelAgeAbove30)), stderror = c(ModelAge30Below$se, ModelAgeAbove30$se),
                  meta = c("0","1"), tau2 = round(c(ModelAge30Below$tau2, ModelAgeAbove30$tau2),3))
age

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=age, digits=3)



#Comparison of Bia
bias <- data.frame(estimate = c(coef(Modelbias), coef(Model)), stderror = c(Modelbias$se, Model$se),
                   meta = c("0","1"), tau2 = round(c(Modelbias$tau2, Model$tau2),3))
bias

rma(estimate, sei=stderror, mods = ~ meta, method="FE", data=bias, digits=3)


######additional analyses
#mean and SD of daily GC in study 10 
# Mean GC dose and SD in subgroup 1 (12 people)
mean1 <- 9.5
sd1 <- 1.8

# Mean GC dose and SD in subgroup 2 (26 people)
mean2 <- 4.4
sd2 <- 1
#
n <- 38
#
totalGC1 <- 12 * mean1
totalGC2 <- 26 * mean2
#
totalGC <- totalGC1 + totalGC2
#mean GC for the whole group
mean <- totalGC / n
#SD for the whole group
sd <- sqrt(((11 * sd1^2) + (25 * sd2^2) + (12 * (mean1 - mean)^2) + (26 * (mean2 - mean)^2)) / (n - 1))

cat("Mean GC in the whole group:", mean, "\n") 
cat("SD of GC in the whole group:", sd)







####bayesian analysis
taupriordensity <-  function(t){dhalfcauchy(t, scale = 25)}
dev.off()

bayesiananalysis <- bayesmeta(y = data1[,"yi"],
                              sigma = data1[,"vi"],
                              labels = data1[,"Author"], 
                              mu.prior.mean=0, mu.prior.sd=50,
                              tau.prior=taupriordensity) 
summary(bayesiananalysis)
forestplot(bayesiananalysis)



