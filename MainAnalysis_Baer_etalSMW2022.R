

####################################
#### TEMPERATURE & SUICIDES ########
####################################

# Code reproducing the study by Barr et al. SMW 2022
# Data available in: https://doi.org/10.48620/38

# LOAD LIBRARIES
library(dplyr); library(lubridate) ;library(dlnm); library(gnm)
library(splines) ; library(mixmeta) ;library(mixmeta); library(rlist)

# SET WORKING DIRECTORY TO FILE DIRECTORY (SESSION >> set working directory >> to source file location)

# LOAD DATA (PREPARED WITH 01.PREP CODE)
load("SuicideData_canton_19952016.RData")

# CREATE LIST REGIONS
listcantons <- c("Lake Geneva", "Midland", "North-West", "Zurich", "East", "Central", "Ticino ")   

# FUNCTION FOR COMPUTING THE Q-AIC
QAIC <- function(model) {
  phi <- summary(model)$dispersion
	loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}

# not needed in the version stored in BORIS
#dta_tmean_sui$tmean <- dta_tmean_sui$Tabsd

#### 1. 1ST STAGE ANALYSIS
# ESTIMATION OF REGION/GROUP SPECIFIC ESTIMATES

res_rr99 <- matrix(NA, nrow=length(listcantons), ncol=8)

cp_list <- list()
coefall <- matrix(NA, nrow=length(listcantons), ncol=2)
vcovall <- list()

listvariables <- names(dta_tmean_sui)[8:15]

# LOOP ACROSS VARIABLES
for (j in seq(length(listvariables))){

  varname <- listvariables[j]

  # LOOP ACROSS CANTONS
  for (i in seq(length(unique(dta_tmean_sui$region)))){
  
    data <- dta_tmean_sui[dta_tmean_sui$region==listcantons[i],]
  
    # PREPARE YOUR TIME VARIABLES
    data$year <- as.numeric(format(data$date,"%Y"))
    data$month <- as.numeric(format(data$date,"%m"))
    data$dow <- factor(weekdays(data$date,abbr=T),
     levels=c("Sun","Mon","Tue","Wed","Thu","Fri","Sat"))

    # CREATE THE STRATUM VARIABLE
    data$stratum <- with(data,factor(year):factor(month):factor(dow):factor(cant))
    
    # TO CHECK THAT THE N COLUM CORRESPONDS TO THE VARNAME
    idcol <- which(names(data)== varname)
  
    # ELIMINATE EMPTY STRATA (FOR CORRECT COMPUTATION OF CI IN gnm)
    ind <- tapply(data[,idcol], data$stratum, sum)[data$stratum]
  
    # CREATE CROSSBASIS 
    cbtmean <- crossbasis(data$tmean, lag=2, 
                          argvar=list(fun="ns",
                                      knots=quantile(data$tmean,c(0.75)))
                          , arglag=list(fun="strata", breaks=c(1)))
    
    # MODEL
    mod <- gnm(as.formula(paste0(varname, " ~ cbtmean")), eliminate=stratum,
                family=quasipoisson(), data=data, na.action="na.exclude",
                subset=ind>0)
    
    # STORE RESULTS
    cp_list[[i]] <- crossreduce(cbtmean,mod, cen=quantile(data$tmean,0.1,na.rm=T),
                                from=round(min(data$tmean),1), to=round(max(data$tmean),1))
    coefall[i,] <-  cp_list[[i]]$coefficients
    vcovall[[i]] <-  cp_list[[i]]$vcov
    
    # Predict at 99
    cp <- crossreduce(cbtmean,mod,at=quantile(data$tmean,0.99,na.rm=T),
                    cen=quantile(data$tmean,0.1,na.rm=T))
    res_rr99[i,j] <- paste0(round(cp$RRfit,2), " [", 
                      round(cp$RRlow,2)," - ",
                      round(cp$RRhigh,2) ,"]")
  }  
  
assign(paste("cp_list",varname,sep = "_"), cp_list) 
assign(paste("coefall", varname,sep = "_"), coefall) 
assign(paste("vcovall", varname,sep = "_"), vcovall) 

}  
  
rm(coefall, vcovall)

rownames(res_rr99) <-  listcantons
colnames(res_rr99) <-  listvariables

#### 2. 2ND STAGE ANALYSIS
# ESTIMATE THE POOL CURVES AND BLUPS

predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))

per <- list.cbind(tapply(dta_tmean_sui$tmean,dta_tmean_sui$region,
                    quantile, probs=predper/100,na.rm=T))
tmeanper <- rowMeans(per)

avertmean <- tapply(dta_tmean_sui$tmean,dta_tmean_sui$region,
                    mean)
rangetmean <- tapply(dta_tmean_sui$tmean,dta_tmean_sui$region,
                    IQR)

# PREDICT POOLED CURVE
bvar <- onebasis(tmeanper, fun="ns", knots=tmeanper[c("75.0%")],
                Boundary.knots=tmeanper[c("0.0%","100.0%")])
cen <- tmeanper["10.0%"]
datanew <- data.frame(avertmean=mean(avertmean),rangetmean=mean(rangetmean))

res_rr99_pool <- matrix(NA, nrow=1, ncol=8)


# OVERALL
mxmodel <- mixmeta(coefall_suicides ~ avertmean + rangetmean ,vcovall_suicides,
                     control=list(showiter=T),method="reml")

datanew <- data.frame(cbind(mean(avertmean),mean(rangetmean)))

colnames(datanew) <- c("avertmean","rangetmean")

pred <- predict(mxmodel, datanew, vcov=T)

cp.pool.comb.overall <- crosspred(bvar, coef=pred$fit, vcov=pred$vcov,
                     model.link="log", at=tmeanper, cen=cen)

    cp <- crosspred(bvar,coef=pred$fit, vcov=pred$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,1] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")

blups.overall <- blup(mxmodel,vcov=T)


# SEX 
coef.sex <- rbind(coefall_sui_sex1,coefall_sui_sex2)
vcov.sex <- c(vcovall_sui_sex1,vcovall_sui_sex2)
indsex <- c(rep(1,7),rep(2,7))

datamxmeta <- data.frame(cbind(rep(avertmean,2),rep(rangetmean,2),indsex))
colnames(datamxmeta) <- c("avertmean","rangetmean","indsex")

mxmodel <- mixmeta(coef.sex ~ avertmean + rangetmean + indsex ,vcov.sex,data=datamxmeta,
                     control=list(showiter=T),method="reml")

datanew <- data.frame(rbind(cbind(mean(avertmean),mean(rangetmean),1),
                      cbind(mean(avertmean),mean(rangetmean),2)))

colnames(datanew) <- c("avertmean","rangetmean","indsex")

pred <- predict(mxmodel, datanew, vcov=T)

cp.pool.comb.sex1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)
cp.pool.comb.sex2 <- crosspred(bvar, coef=pred[[2]]$fit, vcov=pred[[2]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)

    cp <- crosspred(bvar,coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,2] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    cp <- crosspred(bvar,coef=pred[[2]]$fit, vcov=pred[[2]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,3] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")

blups.sex <- blup(mxmodel,vcov=T)


# AGE
coef.age <- rbind(coefall_sui_ageb35,coefall_sui_age3565,coefall_sui_ageab65)
vcov.age <- c(vcovall_sui_ageb35,vcovall_sui_age3565,vcovall_sui_ageab65)
indage <- c(rep(1,7),rep(2,7),rep(3,7))

datamxmeta <- data.frame(cbind(rep(avertmean,3),rep(rangetmean,3),indage))
colnames(datamxmeta) <- c("avertmean","rangetmean","indage")

mxmodel <- mixmeta(coef.age ~ avertmean + rangetmean + as.factor(indage),vcov.age,
                   data=datamxmeta,control=list(showiter=T),method="reml")

datanew <- data.frame(rbind(cbind(mean(avertmean),mean(rangetmean),1),
                      cbind(mean(avertmean),mean(rangetmean),2),
                      cbind(mean(avertmean),mean(rangetmean),3)))

colnames(datanew) <- c("avertmean","rangetmean","indage")

pred <- predict(mxmodel, datanew, vcov=T)

cp.pool.comb.age1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)
cp.pool.comb.age2 <- crosspred(bvar, coef=pred[[2]]$fit, vcov=pred[[2]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)
cp.pool.comb.age3 <- crosspred(bvar, coef=pred[[3]]$fit, vcov=pred[[3]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)

    cp <- crosspred(bvar,coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,4] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    cp <- crosspred(bvar,coef=pred[[2]]$fit, vcov=pred[[2]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,5] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    cp <- crosspred(bvar,coef=pred[[3]]$fit, vcov=pred[[3]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,6] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")

blups.age <- blup(mxmodel,vcov=T)

# METHOD

coef.meth <- rbind(coefall_sui_nonviolent,coefall_sui_violent)
vcov.meth <- c(vcovall_sui_nonviolent,vcovall_sui_violent)
indmeth <- c(rep(1,7),rep(2,7))

datamxmeta <- data.frame(cbind(rep(avertmean,2),rep(rangetmean,2),indmeth))
colnames(datamxmeta) <- c("avertmean","rangetmean","indmeth")

mxmodel <- mixmeta(coef.meth ~ avertmean + rangetmean + indmeth ,vcov.meth,data=datamxmeta,
                     control=list(showiter=T),method="reml")

datanew <- data.frame(rbind(cbind(mean(avertmean),mean(rangetmean),1),
                      cbind(mean(avertmean),mean(rangetmean),2)))

colnames(datanew) <- c("avertmean","rangetmean","indmeth")

pred <- predict(mxmodel, datanew, vcov=T)

cp.pool.comb.meth1 <- crosspred(bvar, coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)
cp.pool.comb.meth2 <- crosspred(bvar, coef=pred[[2]]$fit, vcov=pred[[2]]$vcov,
                     model.link="log", at=tmeanper, cen=cen)

    cp <- crosspred(bvar,coef=pred[[1]]$fit, vcov=pred[[1]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,7] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    cp <- crosspred(bvar,coef=pred[[2]]$fit, vcov=pred[[2]]$vcov,
                     model.link="log",cen=cen,
                     at=tmeanper["99.0%"])
    
    res_rr99_pool[1,8] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")

blups.meth <- blup(mxmodel,vcov=T)



## 3. ESTIMATE THE ER - PREDICT THE BLUP

res_rr99_blup_comb <- matrix(NA, nrow=length(listcantons), ncol=8)
pred_list_comb_overall <- list()
pred_list_comb_sex1 <- pred_list_comb_sex2 <- list()
pred_list_comb_age1 <- pred_list_comb_age2 <- pred_list_comb_age3 <- list()
pred_list_comb_meth1 <- pred_list_comb_meth2 <- list()

# LOOP ACROSS CANTONS
  for (i in seq(length(unique(dta_tmean_sui$region)))){
  
    data <- dta_tmean_sui[dta_tmean_sui$region==listcantons[i],]

    argvar <- list(x=data$tmean,fun="ns",
      knots=quantile(data$tmean,c(0.75),na.rm=T),
      Bound=range(data$tmean,na.rm=T))
    
    bvar <- do.call(onebasis,argvar)
    
    # OVERALL
    pred_list_comb_overall[[i]] <- crosspred(bvar,coef=blups.overall[[i]]$blup,
                     vcov=blups.overall[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)

    cp <- crosspred(bvar,coef=blups.overall[[i]]$blup,
                     vcov=blups.overall[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,1] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # SEX1
    pred_list_comb_sex1[[i]] <- crosspred(bvar,coef=blups.sex[[i]]$blup,
                     vcov=blups.sex[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)

    cp <- crosspred(bvar,coef=blups.sex[[i]]$blup,
                     vcov=blups.sex[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,2] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # SEX2
    pred_list_comb_sex2[[i]] <- crosspred(bvar,coef=blups.sex[[i+7]]$blup,
                     vcov=blups.sex[[i+7]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)
    
    cp <- crosspred(bvar,coef=blups.sex[[i+7]]$blup,
                     vcov=blups.sex[[i+7]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,3] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # AGE 1
    pred_list_comb_age1[[i]] <- crosspred(bvar,coef=blups.age[[i]]$blup,
                     vcov=blups.age[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)
    
    cp <- crosspred(bvar,coef=blups.age[[i]]$blup,
                     vcov=blups.age[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,4] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # AGE 2
    pred_list_comb_age2[[i]] <- crosspred(bvar,coef=blups.age[[i+7]]$blup,
                     vcov=blups.age[[i+7]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)
    
    cp <- crosspred(bvar,coef=blups.age[[i+7]]$blup,
                     vcov=blups.age[[i+7]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,5] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # AGE 3
    pred_list_comb_age3[[i]] <- crosspred(bvar,coef=blups.age[[i+14]]$blup,
                     vcov=blups.age[[i+14]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)
    
    cp <- crosspred(bvar,coef=blups.age[[i+14]]$blup,
                     vcov=blups.age[[i+14]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,6] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # METH 1
    pred_list_comb_meth1[[i]] <- crosspred(bvar,coef=blups.meth[[i]]$blup,
                     vcov=blups.meth[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)
    
    cp <- crosspred(bvar,coef=blups.meth[[i]]$blup,
                     vcov=blups.meth[[i]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,7] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
    
    # METH 2
    pred_list_comb_meth2[[i]] <- crosspred(bvar,coef=blups.meth[[i+7]]$blup,
                     vcov=blups.meth[[i+7]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1)
    
    cp <- crosspred(bvar,coef=blups.meth[[i+7]]$blup,
                     vcov=blups.meth[[i+7]]$vcov,
                     model.link="log",cen=quantile(data$tmean,c(0.1),na.rm=T),
                     by=0.1,at=quantile(data$tmean,0.99,na.rm=T))
    
    res_rr99_blup_comb[i,8] <- paste0(round(cp$allRRfit,2), " [", 
                      round(cp$allRRlow,2)," - ",
                      round(cp$allRRhigh,2) ,"]")
}



rownames(res_rr99_blup_comb) <- listcantons
colnames(res_rr99_blup_comb) <- listvariables

res_rr99_blup_comb_all <- rbind(res_rr99_pool, res_rr99_blup_comb)

rownames(res_rr99_blup_comb_all)[1] <- "pooled"



########## PLOT OF THE E-R FUNCTIONS ##########

collist <- c("black","blue","red","forestgreen","orange","purple",
             "brown",gray(0.4))

coltrans_list <- vector()
for(j in seq(length(collist))){
  coltrans_list[j] <- do.call(rgb,c(as.list(col2rgb(collist[j])),alpha=255/8,max=255))
}

indlab <- predper%in%c(0,1,10,50,90,99,100)

png(filename = "er_oneplot_canton9517_blupscomb.png", 
    width = 2000, height = 2400, units = "px", pointsize = 24,
    bg = "white",  res = NA)


m <- matrix(c(c(1:8),rep(seq(9,16),3),rep(seq(17,24),3),
            rep(seq(25,32),3), rep(seq(33,40),3)),nrow = 8, ncol = 13)
layout(m)
par(mar=c(4,4,2,0),oma=c(1,0,1,0))

# introduce title rotated
titles <- c("Pooled",listcantons)
for (i in seq(length(listcantons)+1)){
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', 
       type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, titles[i], 
     cex = 1.8, col = "black", srt=90, font=2)
} 

# 1) OVERALL
  plot(cp.pool.comb.overall, ylim=c(0.6,3), xlab="Temp. Percentile",col=collist[1],
       ylab="RR (95% CI)", cex.lab=1.4, lwd=2, main="Overall", xaxt="n", cex.main=1.8)
  axis(1,at=tmeanper[indlab],labels=predper[indlab],cex.axis=1)
  abline(v=tmeanper["99.0%"],lty = "dashed", lwd=2)

for (i in seq(length(listcantons)-1)){
  plot(pred_list_comb_overall[[i]], ylim=c(0.6,3), xlab=" ",col=collist[1],
       ylab="RR (95% CI)", cex.lab=1.4, lwd=2, "overall")
  abline(v = per["99.0%",i], lty = "dashed", lwd=2)
  } 
plot(pred_list_comb_overall[[length(listcantons)]], ylim=c(0.6,3), xlab="Temperature (C)",
       col=collist[1], ylab="RR (95% CI)", cex.lab=1.4, lwd=2)
  abline(v = per["99.0%",length(listcantons)], lty = "dashed", lwd=2)

# 2) SEX 
plot(cp.pool.comb.sex1, ylim=c(0.6,3), xlab="Temp. Percentile", ylab="", 
       col=collist[2], ci="area",lwd=2,ci.arg=list(col=coltrans_list[2]),
     cex.lab=1.4, main="By sex", xaxt="n" ,cex.main=1.8)
  axis(1,at=tmeanper[indlab],labels=predper[indlab],cex.axis=1)
  lines(cp.pool.comb.sex2, col=collist[3], lwd=2, ci="area",ci.arg=list(col=coltrans_list[3]))
  abline(v=tmeanper["99.0%"],lty = "dashed", lwd=2)
  
for (i in seq(length(listcantons)-1)){
  plot(pred_list_comb_sex1[[i]], ylim=c(0.6,3), xlab="", ylab="", 
       col=collist[2], ci="area",lwd=2,ci.arg=list(col=coltrans_list[2]),cex.lab=1.4)
  lines(pred_list_comb_sex2[[i]], col=collist[3], lwd=2, ci="area",ci.arg=list(col=coltrans_list[3]))
  abline(v = per["99.0%",i], lty = "dashed", lwd=2)

} 
plot(pred_list_comb_sex1[[length(listcantons)]], ylim=c(0.6,3), xlab="Temperature (C)", ylab="", 
       col=collist[2], ci="area",lwd=2,ci.arg=list(col=coltrans_list[2]),cex.lab=1.4)
  lines(pred_list_comb_sex2[[length(listcantons)]], col=collist[3], lwd=2, ci="area",ci.arg=list(col=coltrans_list[3]))
  abline(v = per["99.0%",length(listcantons)], lty = "dashed", lwd=2)
  legend("topleft",c("Male", "Female"),col=c(collist[2],collist[3]), 
         lwd=1.5, lty=c(1,1), cex=1)


# 3) AGE GROUP

plot(cp.pool.comb.age1,  ylim=c(0.6,3), xlab="Temp. Percentile", ylab="", 
       col=collist[4], ci="area",lwd=2,ci.arg=list(col=coltrans_list[4]),
     cex.lab=1.4, main="By age group", xaxt="n", cex.main=1.8)
  lines(cp.pool.comb.age2, col=collist[5], lwd=2, ci="area",ci.arg=list(col=coltrans_list[5]))
  lines(cp.pool.comb.age3, col=collist[6], lwd=2, ci="area",ci.arg=list(col=coltrans_list[6]))
  axis(1,at=tmeanper[indlab],labels=predper[indlab],cex.axis=1)
  abline(v=tmeanper["99.0%"],lty = "dashed", lwd=2)

      
for (i in seq(length(listcantons)-1)){
  plot(pred_list_comb_age1[[i]],  ylim=c(0.6,3), xlab="", ylab="", 
       col=collist[4], ci="area",lwd=2,ci.arg=list(col=coltrans_list[4]),cex.lab=1.4)
  lines(pred_list_comb_age2[[i]], col=collist[5], lwd=2, ci="area",ci.arg=list(col=coltrans_list[5]))
  lines(pred_list_comb_age3[[i]], col=collist[6], lwd=2, ci="area",ci.arg=list(col=coltrans_list[6]))
  abline(v = per["99.0%",i], lty = "dashed", lwd=2)
} 
plot(pred_list_comb_age1[[length(listcantons)]], ylim=c(0.6,3), xlab="Temperature (C)", ylab="", 
       col=collist[4], ci="area",lwd=1.5,ci.arg=list(col=coltrans_list[4]),cex.lab=1.4)
  lines(pred_list_comb_age2[[length(listcantons)]], col=collist[5], lwd=1.5, 
        ci="area",ci.arg=list(col=coltrans_list[5]))
  lines(pred_list_comb_age3[[length(listcantons)]], col=collist[6], lwd=1.5, 
        ci="area",ci.arg=list(col=coltrans_list[6]))
  abline(v = per["99.0%",length(listcantons)], lty = "dashed", lwd=2)

  legend("topleft",c("<35 years", "35-65 years", ">65 years"),col=c(collist[4],collist[5],collist[6]), 
         lwd=1.5, lty=c(1,1,1), cex=1)

# 4) BY METHOD

plot(cp.pool.comb.meth1, ylim=c(0.6,3), xlab="Temp. Percentile", ylab="",
       col=collist[7], ci="area",lwd=1.5,ci.arg=list(col=coltrans_list[7]), 
     cex.lab=1.4, main="By method",xaxt="n", cex.main=1.8)
  lines(cp.pool.comb.meth2, col=collist[8], lwd=1.5, ci="area",ci.arg=list(col=coltrans_list[8]))
  axis(1,at=tmeanper[indlab],labels=predper[indlab],cex.axis=1)
  abline(v=tmeanper["99.0%"],lty = "dashed", lwd=2)
  
for (i in seq(length(listcantons)-1)){
  plot(pred_list_comb_meth1[[i]], ylim=c(0.6,3), xlab="", ylab="",
       col=collist[7], ci="area",lwd=1.5,ci.arg=list(col=coltrans_list[7]), cex.lab=1.4)
  lines(pred_list_comb_meth2[[i]], col=collist[8], lwd=1.5, ci="area",ci.arg=list(col=coltrans_list[8]))
  abline(v = per["99.0%",i], lty = "dashed", lwd=2)
} 

plot(pred_list_comb_meth1[[length(listcantons)]],  ylim=c(0.6,3), xlab="Temperature (C)", ylab="",
       col=collist[7], ci="area",lwd=1.5,ci.arg=list(col=coltrans_list[7]), cex.lab=1.4)
  lines(pred_list_comb_meth2[[length(listcantons)]], col=collist[8], lwd=1.5, ci="area",ci.arg=list(col=coltrans_list[8]))
  abline(v = per["99.0%",length(listcantons)], lty = "dashed", lwd=2)
  legend("topleft",c("Nonviolent","Violent"),
         col=c(collist[7],collist[8]), lwd=1.5, lty=c(1,1), cex=1)

dev.off() 


