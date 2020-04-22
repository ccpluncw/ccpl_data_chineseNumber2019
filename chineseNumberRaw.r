library(chutils)
library(plyr)
library(nlme)
library(lme4)

### thresholds
minSNpCor <- 0.85
minRT <- .005
maxRT <- .995

#read in data
df.all <- read.table("chinese2.txt", header=T, sep="\t")

data.c <- na.omit(df.all[df.all$cond == "CH1_CHIN",])
data.a <- na.omit(df.all[df.all$cond == "CH1_Eng",])

# get subject numbers - remove extranious random numbers
data.c$sn <- as.numeric(as.character(substr(data.c$sn, nchar(as.character(data.c$sn))-3, nchar(as.character(data.c$sn)))))
data.a$sn <- as.numeric(as.character(substr(data.a$sn, nchar(as.character(data.a$sn))-3, nchar(as.character(data.a$sn)))))
totalChineseSN <- length(unique(data.c$sn))
totalArabicSN <- length(unique(data.a$sn))
data.c$sn <- as.factor(data.c$sn)
data.a$sn <- as.factor(data.a$sn)

allAN <- c(data.c$sn, data.a$sn)

#remove practice trials
data.c <- data.c[data.c$block_c!=1,]
data.a <- data.a[data.a$block_c!=1,]

#fix correct column (two different codings)
data.a$correct2 <- ifelse(data.a$correct == "t" | data.a$correct =="YES", TRUE, FALSE)
data.c$correct2 <- ifelse(data.c$correct == "t" | data.c$correct =="YES", TRUE, FALSE)

data.a$ProbCorrect <- ifelse(data.a$correct2== TRUE, 1, 0)
data.c$ProbCorrect <- ifelse(data.c$correct2== TRUE, 1, 0)

data.c$probe2 <- data.c$probe
#create a parallel probe2 variable in the arabic dataset
data.a$probe2 <- data.a$probe

#prepare the arabic distance and physical similarity dataframe
number <- c(1,2,3,4,6,7,8,9)
PSc <- scale(c(.25, .667, 1.5, .22, .5,.6,.17,.75))
simNoFive <- scale(c(.2,.75,2,1,5,.5,2.5,2))

physical <- data.frame(probe2 = number,PSa = simNoFive, PSc = PSc)
#### now calculate the Numerical Distance from 5 for the encoded digits
physical$dist <- abs(5 - physical$probe2)
physical$WelProb  <- scale(ifelse (physical$probe2 > 5,((physical$probe2)/((physical$probe2)-5)), (5/(5-(physical$probe2)))))
predictorCorrs <- with(physical, cor(data.frame(PSa, PSc, WelProb)))

#merge the physical similarity and distance into the subjects data
data.c.2 <- merge(data.c, physical)
data.a.2 <- merge(data.a, physical)

#log RT
data.a.2$rt <- (data.a.2$rt)
data.c.2$rt <- (data.c.2$rt)

data.a.2$pRT<- rank(data.a.2$rt)/length(data.a.2$rt)
data.c.2$pRT<- rank(data.c.2$rt)/length(data.c.2$rt)

### filter data ###
#remove SNs based on RT criteria
outList <- ch.filterGrpBtwn(data.c.2, "ProbCorrect", "sn", lowThresh = minSNpCor, FUN=mean)
data.c.2 <- outList$datKeptRaw
numSN.belowPCorThresh.c <- outList$numRemoved
sn.removed.belowPCorThresh.c	<- outList$datRemoved

sn.rates.c	<- outList$datKept
sn.rates.c$errC	<- 1-sn.rates.c$out
sn.mean.errorRate.c	<- mean(sn.rates.c$errC)
sn.SD.errorRate.c	<- sd(sn.rates.c$errC)


outList <- ch.filterGrpBtwn(data.a.2, "ProbCorrect", "sn", lowThresh = minSNpCor, FUN=mean)
data.a.2 <- outList$datKeptRaw
numSN.belowPCorThresh.a <- outList$numRemoved
sn.removed.belowPCorThresh.a	<- outList$datRemoved

sn.rates.a	<- outList$datKept
sn.rates.a$errA	<- 1-sn.rates.a$out
sn.mean.errorRate.a	<- mean(sn.rates.a$errA)
sn.SD.errorRate.a	<- sd(sn.rates.a$errA)

#ttest on error rate difference
sn.all.rates <- merge (sn.rates.a, sn.rates.c, by='sn')
sn.err.ttest <- with(sn.all.rates, t.test(errC, errA, paired=T))


#remove individual trials based on RT criteria above
outList <- ch.filterDataBetween(data.c.2, "pRT", minRT, maxRT)
data.c.2 <- outList$datKept
pRTremoved.c <- outList$pRemoved

outList <- ch.filterDataBetween(data.a.2, "pRT", minRT, maxRT)
data.a.2 <- outList$datKept
pRTremoved.a <- outList$pRemoved

#average RT analysis
sn.RT.c <- ddply(data.c.2, .(sn), summarize, mRTc = mean(rt))
sn.RT.a <- ddply(data.a.2, .(sn), summarize, mRTa = mean(rt))
sn.all.RT <- merge (sn.RT.a, sn.RT.c, by='sn')
sn.mean.RT.c <- mean(sn.RT.c$mRTc)
sn.SD.RT.c <- sd(sn.RT.c$mRTc)
sn.mean.RT.a <- mean(sn.RT.a$mRTa)
sn.SD.RT.a <- sd(sn.RT.a$mRTa)

sn.RT.ttest <- with(sn.all.RT, t.test(mRTc, mRTa, paired=T))

###Mixed Model analysis
  data.c.lme.wp <- lme(rt ~ WelProb + PSa + PSc,random = ~1|sn/(WelProb/PSa/PSc), data=data.c.2[data.c.2$correct2,])
  data.a.lme.wp <- lme(rt ~ WelProb + PSa + PSc,random = ~1|sn/(WelProb/PSa/PSc), data=data.a.2[data.a.2$correct2,])

  data.a.lme.dat <- data.a.2[data.a.2$correct2,]
  data.a.lme.dat$fitted <- fitted(data.a.lme.wp)
  data.a.lme.wp.r2 <- cor(data.a.lme.dat$fitted, data.a.lme.dat$rt)^2

  data.c.lme.dat <- data.c.2[data.c.2$correct2,]
  data.c.lme.dat$fitted <- fitted(data.c.lme.wp)
  data.c.lme.wp.r2 <- cor(data.c.lme.dat$fitted, data.c.lme.dat$rt)^2

  data.c.lmer.wp <- lmer(rt ~ WelProb + PSa  + PSc + (1|sn) + (1|sn:WelProb), data=data.c.2[data.c.2$correct2,])
  data.a.lmer.wp <- lmer(rt ~ WelProb + PSa  + PSc + (1|sn) + (1|sn:WelProb), data=data.a.2[data.a.2$correct2,])

  sink("chin_final.raw.out.txt", append=F)
    cat("\n\n ***** Total Chinese SN\n\n")
    print(totalChineseSN)
    cat("\n\n ***** Total Arabic SN\n\n")
    print(totalArabicSN)

    cat("\n\n ***** Min SN p(Correct)\n\n")
    print(minSNpCor)
    cat("\n\n ***** Num SN Removed Chinese\n\n")
    print(numSN.belowPCorThresh.c)
    cat("\n\n ***** SNs Removed Chinese\n\n")
    print(sn.removed.belowPCorThresh.c)

    cat("\n\n ***** Num SN Removed Arabic\n\n")
    print(numSN.belowPCorThresh.a)
    cat("\n\n ***** SNs Removed Arabic\n\n")
    print(sn.removed.belowPCorThresh.a)

    cat("\n\n ***** Min RT\n\n")
    print(minRT)
    cat("\n\n ***** Max RT\n\n")
    print(maxRT)
    cat("\n\n ***** p(Num Min RT Removed) Chinese\n\n")
    print(pRTremoved.c)
    cat("\n\n ***** p(Num Max RT Removed) Arabic\n\n")
    print(pRTremoved.a)

		cat("\n\n********************************************\n\n")
    cat("\n\n ***** SNs Mean Error Rate Chinese\n\n")
    print(sn.mean.errorRate.c)
    cat("\n\n ***** SNs SD Error Rate Chinese\n\n")
    print(sn.SD.errorRate.c)
    cat("\n\n ***** SNs Mean Error Rate Arabic\n\n")
    print(sn.mean.errorRate.a)
    cat("\n\n ***** SNs SD Error Rate Arabic\n\n")
    print(sn.SD.errorRate.a)

    cat("\n\n ***** ttest on error rate difference Chinese vs Arabic\n\n")
    print(sn.err.ttest)

    cat("\n\n********************************************\n\n")
    cat("\n\n ***** SNs Mean RT Chinese\n\n")
    print(sn.mean.RT.c)
    cat("\n\n ***** SNs SD RT Chinese\n\n")
    print(sn.SD.RT.c)
    cat("\n\n ***** SNs Mean RT Arabic\n\n")
    print(sn.mean.RT.a)
    cat("\n\n ***** SNs SD RT Arabic\n\n")
    print(sn.SD.RT.a)

    cat("\n\n ***** ttest on RT difference Chinese vs Arabic\n\n")
    print(sn.RT.ttest)

    cat("\n\n********************************************\n\n")
    cat("\n\n ***** Correlation between PSa and PSc and Welford\n\n")
    print(predictorCorrs)

    cat("\n\n********************************************\n\n")
    cat("\n\n ***** nls LME Parameters for fit to a\n\n")
    print(fixef(data.a.lme.wp))
    cat("\n\n ***** nls LME SEs for parameters for fit to a\n\n")
    print(sqrt(diag(vcov(data.a.lme.wp))))
    cat("\n\n ***** nls LME ANOVA for fit to a\n\n")
    print(anova.lme(data.a.lme.wp))
    cat("\n\n ***** nls LME R2 for fit to a\n\n")
    print(data.a.lme.wp.r2)


    cat("\n\n ***** nls LME Parameters for fit to c\n\n")
    print(fixef(data.c.lme.wp))
    cat("\n\n ***** nls LME SEs for parameters for fit to c\n\n")
    print(sqrt(diag(vcov(data.c.lme.wp))))
    cat("\n\n ***** nls LME ANOVA for fit to c\n\n")
    print(anova.lme(data.c.lme.wp))
    cat("\n\n ***** nls LME R2 for fit to a\n\n")
    print(data.c.lme.wp.r2)
    cat("\n\n********************************************\n\n")
    cat("\n\n ***** nls LMER Parameters for fit to a\n\n")
    print(fixef(data.a.lmer.wp))
    cat("\n\n ***** nls LMER SEs for parameters for fit to a\n\n")
    print(sqrt(diag(vcov(data.a.lmer.wp))))
    cat("\n\n ***** nls LMER ANOVA for fit to a\n\n")
    print(anova(data.a.lmer.wp))
    cat("\n\n ***** nls LMER Parameters for fit to c\n\n")
    print(fixef(data.c.lmer.wp))
    cat("\n\n ***** nls LMER SEs for parameters for fit to c\n\n")
    print(sqrt(diag(vcov(data.c.lmer.wp))))
    cat("\n\n ***** nls LMER ANOVA for fit to c\n\n")
    print(anova(data.c.lmer.wp))


  sink(NULL)

  ###Individual analysis by subject
  subs.c <- unique(data.c.2$sn)
  subs.a <- unique(data.a.2$sn)
  sink("chin_final.raw.out.txt", append=T)
    print("******************individual subjects chinese - chinese ********************")
  sink(NULL)
  outDF <- NULL

  for (j in 1:2) {
    if (j == 1) {
      subs <- subs.c
      set <- "chi-chi"
      data <- data.c.2
    } else {
      subs <- subs.a
      set <- "chi-arb"
      data <- data.a.2
    }
    for(i in subs) {
      data.tmp <- data[data$correct2 & data$sn == i,]
      data.lm.wp <- lm(rt ~ WelProb + PSa + PSc, data=data.tmp)
      tmp.lm.df <- data.frame(sn = i, set = set,
          inter = summary(data.lm.wp)$coefficients[1,1],
          Wel = summary(data.lm.wp)$coefficients[2,1],
          PSa = summary(data.lm.wp)$coefficients[3,1],
          PSc = summary(data.lm.wp)$coefficients[4,1],
          interSE = summary(data.lm.wp)$coefficients[1,2],
          WelSE = summary(data.lm.wp)$coefficients[2,2],
          PSaSE = summary(data.lm.wp)$coefficients[3,2],
          PScSE = summary(data.lm.wp)$coefficients[4,2],
          interT = summary(data.lm.wp)$coefficients[1,3],
          WelT = summary(data.lm.wp)$coefficients[2,3],
          PSaT = summary(data.lm.wp)$coefficients[3,3],
          PScT = summary(data.lm.wp)$coefficients[4,3],
          interP = summary(data.lm.wp)$coefficients[1,4],
          WelP = summary(data.lm.wp)$coefficients[2,4],
          PSaP = summary(data.lm.wp)$coefficients[3,4],
          PScP = summary(data.lm.wp)$coefficients[4,4])

       if(is.null(outDF)) {
         outDF <- tmp.lm.df
       } else {
         outDF <- rbind(outDF, tmp.lm.df)
       }

      sink("chin1.raw.out.txt", append=T)
        cat("\n\n******************subject: ", i, set, "********************\n\n")
        cat("\n\n ***** lm Parameters\n\n")
        print(summary(data.lm.wp))
      sink(NULL)
    }
  }

  # prepare dataset for correlations between predictor variables.
  chiArb <- outDF[outDF$set=="chi-arb",]
  chiChi <- outDF[outDF$set=="chi-chi",]
  data.wide <- merge(chiArb,chiChi, by="sn", suffixes = c(".chiArb",".chiChi"))
  data.wide <- data.wide[,-c(2,19)]

  sink("chin_final.raw.out.txt", append=T)
    cat("\n\n****************** Correlations between PS parameters ********************\n\n")
		cat("\n\n ***** parameters for: PS arabic in Chi-Chi with PS chinese in Chi-Chi\n\n")
    print(with(data.wide, cor.test(PSa.chiChi, PSc.chiChi)))
		cat("\n\n ***** parameters for: PS arabic in Chi-Arb with PS chinese in Chi-Arb\n\n")
    print(with(data.wide, cor.test(PSa.chiArb, PSc.chiArb)))
    cat("\n\n ***** parameters for: PS arabic in Chi-Arb with PS arabic in Chi-Chi\n\n")
    print(with(data.wide, cor.test(PSa.chiArb, PSa.chiChi)))
    cat("\n\n ***** parameters for: PS arabic in Chi-Arb with PS chinese in Chi-Chi\n\n")
    print(with(data.wide, cor.test(PSa.chiArb, PSc.chiChi)))
    cat("\n\n ***** parameters for: PS arabic in Chi-Arb with PS chinese in Chi-Chi\n\n")
    print(with(data.wide, cor.test(PSc.chiArb, PSa.chiChi)))
    cat("\n\n ***** parameters for: PS chinese in Chi-Arb with PS chinese in Chi-Chi\n\n")
    print(with(data.wide, cor.test(PSc.chiArb, PSc.chiChi)))
  sink(NULL)
