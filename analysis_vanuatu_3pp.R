# To run at the command line, replace this with a function that calls setwd()
# to the local directory where you have the data and R files.
setwd_vanuatu_3pp()

library(splines)
library(gmodels)
library(epitools)

d2 <- read.table("data_vanuatu_3pp.tsv", sep="\t", fill=T, head=T, na.strings="*")
d2$AntiPun <- as.factor(ifelse(d2$AntiGot == "Bad", "Y", "N"))
d2$AntiEmpty <- as.factor(ifelse(d2$AntiGot == "Empty", "Y", "N"))
d2$AntiGood <- as.factor(ifelse(d2$AntiGot == "Good", "Y", "N"))
d2$NeutPun <- as.factor(ifelse(d2$NeutGot == "Bad", "Y", "N"))
d2$NeutEmpty <- as.factor(ifelse(d2$NeutGot == "Empty", "Y", "N"))
d2$NeutGood <- as.factor(ifelse(d2$NeutGot == "Good", "Y", "N"))

d2$AntiPunNum <- ifelse(d2$AntiGot == "Bad", 1, 0)
d2$AntiEmptyNum <- ifelse(d2$AntiGot == "Empty", 1, 0)
d2$AntiGoodNum <- ifelse(d2$AntiGot == "Good", 1, 0)
d2$NeutPunNum <- ifelse(d2$NeutGot == "Bad", 1, 0)
d2$NeutEmptyNum <- ifelse(d2$NeutGot == "Empty", 1, 0)
d2$NeutGoodNum <- ifelse(d2$NeutGot == "Good", 1, 0)

behaviours <- c("AntiPun", "NeutPun", "AntiEmpty", "NeutEmpty", "AntiGood", "NeutGood")
behavsOrder2 <- c("AntiPun", "AntiEmpty", "AntiGood", "NeutPun", "NeutEmpty", "NeutGood")
niceGreen <- "#009600FF"
s <- seq(4.5,9.5,by=1)
critval <- qnorm(0.975)

########################################### Plotting functions
     
plotTwoGreenRed <- function( mod1, mod2, predDat1, predDat2 ) {    
  plot(0,0,xlim=c(4,10), ylim=c(0,1),xlab="",ylab="" )
  plotPropWithConf(mod1, predDat1, niceGreen )
  plotPropWithConf(mod2, predDat2, "Red", lineType = 'longdash' )
}

plotPropWithConf <- function(model, predData, colour, alphaF=.1,lineType='solid' ) {
  preds <- predict(model, newdata = predData, type = "link", se.fit = TRUE)
  upr <- preds$fit + (critval * preds$se.fit)
  lwr <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
  fit2 <- model$family$linkinv(fit)
  upr2 <- model$family$linkinv(upr)
  lwr2 <- model$family$linkinv(lwr)
  topOfPoly <- data.frame( x = predData$Age, y = upr2 )
  topOfPoly <- topOfPoly[order(topOfPoly$x),]
  bottomOfPoly <- data.frame( x = predData$Age, y = lwr2 )
  bottomOfPoly <- bottomOfPoly[order(-topOfPoly$x),]
  polygon(x=c(topOfPoly$x,bottomOfPoly$x), y=c(topOfPoly$y,bottomOfPoly$y),col=adjustcolor(colour,alpha.f=alphaF),border=NA)
  lines(x=predData$Age,y=fit2, col=colour, lwd=2, lty=lineType)
  lines(x=predData$Age,y=upr2, col=colour, lty="dotted")
  lines(x=predData$Age,y=lwr2, col=colour, lty="dotted")
}

plotTwoPointsGreenRed <- function( DVGreen, DVRed, whichDatGreen, whichDatRed, sampleSize = F ) {
  meansGreen <- sapply(s,function(x) {
    datForMean <- paste( 'd2$', DVGreen, '[', whichDatGreen, '& d2$Age-', x, '<=.5 & d2$Age-', x, '>-.5 ]', sep = "" )
    mean(eval(parse(text=datForMean)))
  } )
  meansRed <- sapply(s,function(x) {
    datForMean <- paste( 'd2$', DVRed, '[', whichDatRed, '& d2$Age-', x, '<=.5 & d2$Age-', x, '>-.5 ]', sep = "" )
    mean(eval(parse(text=datForMean)))
  } )
  nsGreen <- sapply(s,function(x) {
    datForMean <- paste( 'd2$', DVGreen, '[', whichDatGreen, '& d2$Age-', x, '<=.5 & d2$Age-', x, '>-.5 ]', sep = "" )
    length(eval(parse(text=datForMean)))
  } )
  nsRed <- sapply(s,function(x) {
    datForMean <- paste( 'd2$', DVRed, '[', whichDatRed, '& d2$Age-', x, '<=.5 & d2$Age-', x, '>-.5 ]', sep = "" )
    length(eval(parse(text=datForMean)))
  } )
  points(s,meansGreen,col=niceGreen,pch=1)
  points(s,meansRed,col="red",pch=4)
  if( sampleSize ) {
    text(s-.3,meansGreen,nsGreen,col=niceGreen)
    text(s+.3,meansRed,nsRed,col="red")
  }
}

########################################### Output participant age table for methods section

ages <- seq(4,10,by=1)
ur <- subset( d2, Location == "Urban")
tableAges <- data.frame(
  Age = ages,
  FreePub = sapply(ages,function(x) length(ur$AntiPunNum[ ur$Cost == "N" & ur$Condition == "Public" & floor(ur$Age)==x])),
  FreePri = sapply(ages,function(x) length(ur$AntiPunNum[ ur$Cost == "N" & ur$Condition == "Private" & floor(ur$Age)==x])),
  CostPub = sapply(ages,function(x) length(ur$AntiPunNum[ ur$Cost == "Y" & ur$Condition == "Public" & floor(ur$Age)==x])),
  CostPri = sapply(ages,function(x) length(ur$AntiPunNum[ ur$Cost == "Y" & ur$Condition == "Private" & floor(ur$Age)==x]))
)
tableAges
sum(tableAges[,2:5])
length(d2[,1])

rur <- subset( d2, Location == "Rural")
table(cut(rur$Age, breaks = c(ages,11), right = FALSE, include.lowest = TRUE))

########################################### Punish antisocial more?

punTab <- rbind( table(d2$NeutPun), table(d2$AntiPun) )
punTab
round( punTab[,2] / ( punTab[,1] + punTab[,2] ), 2 )
fisher.test(punTab)

d2Age4 <- d2[d2$Age>=4&d2$Age<5,]
punTabAge4 <- rbind( table(d2Age4$NeutPun), table(d2Age4$AntiPun) )
punTabAge4
round( punTabAge4[,2] / ( punTabAge4[,1] + punTabAge4[,2] ), 2 )
sum(punTabAge4)
fisher.test(punTabAge4)

d2Age5 <- d2[d2$Age>=5&d2$Age<6,]
punTabAge5 <- rbind( table(d2Age5$NeutPun), table(d2Age5$AntiPun) )
punTabAge5
round( punTabAge5[,2] / ( punTabAge5[,1] + punTabAge5[,2] ), 2 )
sum(punTabAge5)
fisher.test(punTabAge5)

d2Age6Plus <- d2[d2$Age>=6,]
punTabAge6Plus <- rbind( table(d2Age6Plus$NeutPun), table(d2Age6Plus$AntiPun) )
punTabAge6Plus
round( punTabAge6Plus[,2] / ( punTabAge6Plus[,1] + punTabAge6Plus[,2] ), 2 )
sum(punTabAge6Plus)
fisher.test(punTabAge6Plus)

########################################### Calculate allocation variable summary statistics with CIs

ageForTable <- median(d2$Age)
ageForTable
results <- lapply( behaviours, function(beh) {
  modForPreds <- glm(as.formula(paste(beh,"~Condition+Cost+Age")), data=d2, family="binomial")
  predVals <- data.frame( Age = rep(ageForTable,4), Condition = c("Private","Public","Private","Public"), Cost = c("N","N","Y","Y"))
  preds <- predict(modForPreds, newdata = predVals, type = "link", se.fit = TRUE)
  upr <- preds$fit + (critval * preds$se.fit)
  lwr <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
  fit2 <- modForPreds$family$linkinv(fit)
  upr2 <- modForPreds$family$linkinv(upr)
  lwr2 <- modForPreds$family$linkinv(lwr)
  data.frame( behav = beh, condition = predVals$Condition, cost = predVals$Cost, prop = fit2, ci95l = lwr2, ci95u = upr2 )
})
results <- Reduce(function(f1, f2) rbind(f1, f2),results)
results

# Note that the graph corresponding to this is made in Excel from the above table.
# R appears to have no convenient way to make nested-category x-axis labels of the type used.

######################################### Makes main model table

meanAge <- mean(d2$Age)
d2$AgeCentred <- d2$Age - meanAge

mods <- c()
mods2 <- c()
pVals <- c()
coefs <- c()
coefAndPValTab <- c()
nothing <- lapply( behavsOrder2, function(beh) {
  cat("\n\n",beh,"\n")
  mods[[beh]] <<- glm(as.formula(paste(beh,"~Condition+Cost+Sex+AgeCentred")), data=d2, family="binomial")
  printCoefmat(summary(mods[[beh]])$coefficients)
  mods2[[beh]] <<- update(mods[[beh]], . ~ . + Condition:AgeCentred + Sex:AgeCentred )
  #NB adding the three-way interaction changes little
  printCoefmat(summary(mods2[[beh]])$coefficients)
  print(anova(mods[[beh]],mods2[[beh]],test="Chisq")[5])
  
  pVals[[beh]] <<- c(coef(summary(mods[[beh]]))[1:5,4],
    anova(mods[[beh]],mods2[[beh]],test="Chisq")$"Pr(>Chi)"[2],
    coef(summary(mods2[[beh]]))[6:7,4]
  )
  coefs[[beh]] <<- c(coef(summary(mods[[beh]]))[1:5,1], NA, coef(summary(mods2[[beh]]))[6:7,1]) 
  coefAndPValTab[[beh]] <<- cbind( coefs[[beh]], pVals[[beh]] )
})
tab <- Reduce(function(f1, f2) cbind(f1, f2),coefAndPValTab)
colnames(tab) <- paste(rep(behavsOrder2,each=2),rep(c("C","P"),6),sep="")
round(tab,digits=3) # with the empty box model which we don't report except to visualise in the figure
round(tab[, -c(3,4,9,10)], digits=3)

######################## Spline model plot of public/private against age

png("splineFig2.png", width = 5, height = 5, units = "in", res = 300)

# Create knot locations based on centered age
knotLoc <- quantile(d2$AgeCentred, probs=c(1/3, 2/3))

# Fit model with centered age
modKnot <- glm(AntiPun ~ Condition + ns(AgeCentred, knots=knotLoc) + 
    Condition:ns(AgeCentred, knots=knotLoc) + Cost,
    data=d2, family="binomial"
)
summary(modKnot)

# Create prediction data with BOTH centered age (for model) and raw age (for plotting)
preddatPublic <- data.frame( 
    AgeCentred = seq(4 - meanAge, 10 - meanAge, length = 100),  # Centered values for model prediction
    Age = seq(4, 10, length = 100),  # Raw values for x-axis in plot
    Condition = rep("Public", 100), 
    Cost = rep("N", 100) 
)
preddatPrivate <- data.frame( 
    AgeCentred = seq(4 - meanAge, 10 - meanAge, length = 100),  # Centered values for model prediction
    Age = seq(4, 10, length = 100),  # Raw values for x-axis in plot
    Condition = rep("Private", 100), 
    Cost = rep("N", 100) 
)

plotTwoGreenRed(modKnot, modKnot, preddatPrivate, preddatPublic)
plotTwoPointsGreenRed('AntiPunNum', 'AntiPunNum', 
    'd2$Condition == "Private" & d2$Cost == "N"', 
    'd2$Condition == "Public" & d2$Cost == "N"')

mtext("Age (years)", side=1, line=2) 
mtext("Proportion participants punishing antisocial actor", side=2, line=2) 
legend(7.1, .25, legend=c("Anonymous", "In person"), 
    col=c(niceGreen, "Red"), 
    lty=c("solid", "longdash"), 
    lwd=c(2, 2))

dev.off()

############################## Investigating rewarding of antisocial actor

cat("Rewarded antisocial actor:",
    round(100 * mean(d2$AntiGood[d2$NeutPun == "Y"] == "Y"), 1), "% (punished neutral) vs",
    round(100 * mean(d2$AntiGood[d2$NeutPun == "N"] == "Y"), 1), "% (didn't punish neutral)\n")


AntiGoodMod3 <- update(mods2[["AntiGood"]], . ~ . + NeutPun)
summary(mods2[["AntiGood"]])
summary(AntiGoodMod3)
anova(mods2[["AntiGood"]],AntiGoodMod3,test="Chisq")



# AntiGood by Condition and Age
knotLoc <- quantile(d2$AgeCentred, probs=c(1/3, 2/3))  # Use centered age for knots
AntiGoodModCondAge <- glm(AntiGood ~ Condition + ns(AgeCentred, knots=knotLoc) + 
                          Condition:ns(AgeCentred, knots=knotLoc), 
                          data=d2, family="binomial")
summary(AntiGoodModCondAge)

preddat1 <- data.frame( 
    AgeCentred = seq(4 - meanAge, 10 - meanAge, length = 100),
    Age = seq(4, 10, length = 100),
    Condition = rep("Private", 100) 
)
preddat2 <- data.frame( 
    AgeCentred = seq(4 - meanAge, 10 - meanAge, length = 100),
    Age = seq(4, 10, length = 100),
    Condition = rep("Public", 100) 
)
png("splineFig3a.png", width = 5, height = 5, units = "in", res = 300)
plotTwoGreenRed(AntiGoodModCondAge, AntiGoodModCondAge, preddat1, preddat2)
plotTwoPointsGreenRed('AntiGoodNum', 'AntiGoodNum', 'd2$Condition=="Private"', 'd2$Condition=="Public"')
mtext("(a)",at=4,line=.5) 
mtext("Age (years)", side=1, line=2) 
mtext("Proportion participants rewarding antisocial actor", side=2, line=2) 
legend(6, 1, legend=c("In person", "Anonymous"), col=c("Red", niceGreen), lty=c("solid", "longdash"), lwd=c(2, 2))
dev.off()

  
AntiGoodModSexAge <- glm(AntiGood~Sex+AgeCentred+Sex:AgeCentred, data=d2, family="binomial")
summary(AntiGoodModSexAge)
preddat1 <- data.frame( 
    AgeCentred = seq(4 - meanAge, 10 - meanAge, length = 100),
    Age = seq(4, 10, length = 100),
    Sex = rep("F",100) 
)
preddat2 <- data.frame( 
    AgeCentred = seq(4 - meanAge, 10 - meanAge, length = 100),
    Age = seq(4, 10, length = 100),
    Sex = rep("M",100) 
)
png("splineFig3b.png", width = 5, height = 5, units = "in", res = 300)
plotTwoGreenRed(AntiGoodModSexAge,AntiGoodModSexAge,preddat1,preddat2)
plotTwoPointsGreenRed('AntiGoodNum','AntiGoodNum','d2$Sex=="F"','d2$Sex=="M"' )
mtext("(b)",at=4,line=.5) 
#mtext("Proportion participants rewarding antisocial actor", side=2, line=2) 
mtext("Age (years)",side=1,line=2) 
legend(6,1,legend=c("Boys","Girls"),col=c("Red",niceGreen), lty=c("solid","longdash"),lwd=c(2,2))
dev.off()

#######################################################
############################### Urban rural comparison
#######################################################

d2Rural <- subset(d2, Location == "Rural" )
d2UrbPrivN <- subset(d2, Location == "Urban" & Condition == "Private" & Cost == "N" ) 
length(subset(d2UrbPrivN,Age>=min(d2Rural$Age) & Age<=max(d2Rural$Age))$Age)

tagNearest <- function(x) {
  d2UrbPrivN$diff <<- abs(d2UrbPrivN$Age - x)
  d2UrbPrivN$tagged[which(d2UrbPrivN$diff==min(subset(d2UrbPrivN, tagged==FALSE)$diff) & d2UrbPrivN$tagged == FALSE)[1]]<<-TRUE
  diffList <<- append( diffList, min(subset(d2UrbPrivN, tagged==FALSE)$diff))
}

sampleAgain <- function() {
  forget <- lapply(d2Rural$Age,tagNearest)
  urbanComp<<-subset(d2UrbPrivN, tagged==TRUE)
  cat(i, length(d2Rural$Age), mean(d2Rural$Age),sd(d2Rural$Age), min(d2Rural$Age), max(d2Rural$Age), "\n")
  cat(i, length(urbanComp$Age), mean(urbanComp$Age),sd(urbanComp$Age), min(urbanComp$Age), max(urbanComp$Age), "\n")
  cat(mean(diffList),"\n\n")
}

d2UrbPrivN$tagged = FALSE
diffList <- c()
for(i in 1:4) {
  sampleAgain()
}

# Function to extract proportion and CI from prop.test
extract_prop_ci <- function(outcome_var, sample_name) {
  counts <- summary(outcome_var)
  successes <- counts[2]
  total <- counts[1] + counts[2]
  prop_result <- prop.test(successes, total)
  proportion <- as.numeric(prop_result$estimate)
  ci_lower <- prop_result$conf.int[1]
  ci_upper <- prop_result$conf.int[2]
  result <- c(proportion, ci_lower, ci_upper)
  names(result) <- paste0(sample_name, c("_prop", "_ci_lower", "_ci_upper"))
  return(result)
}

# Calculate proportions and CIs for all conditions
results_list <- list()

# Rural sample results
rural_antipun <- extract_prop_ci(d2Rural$AntiPun, "Rural")
rural_antigood <- extract_prop_ci(d2Rural$AntiGood, "Rural")
rural_neutpun <- extract_prop_ci(d2Rural$NeutPun, "Rural")
rural_neutgood <- extract_prop_ci(d2Rural$NeutGood, "Rural")

# Urban sample results (using final urbanComp from the matching procedure)
urban_antipun <- extract_prop_ci(urbanComp$AntiPun, "Urban")
urban_antigood <- extract_prop_ci(urbanComp$AntiGood, "Urban")
urban_neutpun <- extract_prop_ci(urbanComp$NeutPun, "Urban")
urban_neutgood <- extract_prop_ci(urbanComp$NeutGood, "Urban")

# Build the final results table
results_table <- data.frame(
  Actor = c("Antisocial", "Antisocial", "Neutral", "Neutral"),
  Response = c("Punishment", "Reward", "Punishment", "Reward"),
  Urban_Proportion = c(urban_antipun[1], urban_antigood[1], urban_neutpun[1], urban_neutgood[1]),
  Urban_CI_Lower = c(urban_antipun[2], urban_antigood[2], urban_neutpun[2], urban_neutgood[2]),
  Urban_CI_Upper = c(urban_antipun[3], urban_antigood[3], urban_neutpun[3], urban_neutgood[3]),
  Rural_Proportion = c(rural_antipun[1], rural_antigood[1], rural_neutpun[1], rural_neutgood[1]),
  Rural_CI_Lower = c(rural_antipun[2], rural_antigood[2], rural_neutpun[2], rural_neutgood[2]),
  Rural_CI_Upper = c(rural_antipun[3], rural_antigood[3], rural_neutpun[3], rural_neutgood[3])
)

results_table[, 3:8] <- round(results_table[, 3:8], 4)
cat("\n=== FINAL RESULTS TABLE ===\n")
print(results_table, row.names = FALSE)

cat("\n=== FISHER'S EXACT TESTS FOR GROUP COMPARISONS ===\n")
fisher_antipun <- fisher.test(matrix(c(summary(d2Rural$AntiPun), summary(urbanComp$AntiPun)), 2, 2))
fisher_antigood <- fisher.test(matrix(c(summary(d2Rural$AntiGood), summary(urbanComp$AntiGood)), 2, 2))
fisher_neutpun <- fisher.test(matrix(c(summary(d2Rural$NeutPun), summary(urbanComp$NeutPun)), 2, 2))
fisher_neutgood <- fisher.test(matrix(c(summary(d2Rural$NeutGood), summary(urbanComp$NeutGood)), 2, 2))
cat("Antisocial Punishment p-value:", round(fisher_antipun$p.value, 4), "\n")
cat("Antisocial Reward p-value:", round(fisher_antigood$p.value, 4), "\n") 
cat("Neutral Punishment p-value:", round(fisher_neutpun$p.value, 4), "\n")
cat("Neutral Reward p-value:", round(fisher_neutgood$p.value, 4), "\n")

# Age matching verification - just the numbers
cat("Rural (n =", length(d2Rural$Age), "): mean =", round(mean(d2Rural$Age), 1), 
    ", range =", round(min(d2Rural$Age), 1), "to", round(max(d2Rural$Age), 1), "\n")
cat("Urban (n =", length(urbanComp$Age), "): mean =", round(mean(urbanComp$Age), 1), 
    ", range =", round(min(urbanComp$Age), 1), "to", round(max(urbanComp$Age), 1), "\n")

#############################################
# Cross-cultural comparison: Sweden vs Vanuatu (Urban, No-Cost)
#############################################

sweden_anon_anti_data <- factor(c(rep("N", 10), rep("Y", 14)), levels = c("N", "Y"))
sweden_public_anti_data <- factor(c(rep("N", 19), rep("Y", 5)), levels = c("N", "Y"))
sweden_anon_neut_data <- factor(c(rep("N", 23), rep("Y", 1)), levels = c("N", "Y"))
sweden_public_neut_data <- factor(c(rep("N", 23), rep("Y", 1)), levels = c("N", "Y"))

# Calculate Swedish proportions and CIs
sweden_anon_anti <- extract_prop_ci(sweden_anon_anti_data, "Sweden_Anon_Anti")
sweden_public_anti <- extract_prop_ci(sweden_public_anti_data, "Sweden_Public_Anti")
sweden_anon_neut <- extract_prop_ci(sweden_anon_neut_data, "Sweden_Anon_Neut")
sweden_public_neut <- extract_prop_ci(sweden_public_neut_data, "Sweden_Public_Neut")

# Vanuatu model fitting
SwedishAge <- 5 + 2/12
AntiPunModFS <- glm(AntiPun ~ Condition + Cost + Age + Condition:Age, data = d2, family = "binomial")
NeutPunModFS <- glm(NeutPun ~ Condition + Cost + Age + Condition:Age, data = d2, family = "binomial")

# Prediction data for Vanuatu (at Swedish age, no-cost condition)
predData <- data.frame(Age = rep(SwedishAge, 2), Cost = rep("N", 2), Condition = c("Private", "Public"))

# Get predictions and SEs
predsAntiPun <- predict(AntiPunModFS, newdata = predData, type = "link", se.fit = TRUE)
predsNeutPun <- predict(NeutPunModFS, newdata = predData, type = "link", se.fit = TRUE)

# Calculate Vanuatu fitted values

# Antisocial punishment
upr_anti <- predsAntiPun$fit + (critval * predsAntiPun$se.fit)
lwr_anti <- predsAntiPun$fit - (critval * predsAntiPun$se.fit)
vanuatu_anti_fit <- AntiPunModFS$family$linkinv(predsAntiPun$fit)
vanuatu_anti_lwr <- AntiPunModFS$family$linkinv(lwr_anti)
vanuatu_anti_upr <- AntiPunModFS$family$linkinv(upr_anti)

# Neutral punishment
upr_neut <- predsNeutPun$fit + (critval * predsNeutPun$se.fit)
lwr_neut <- predsNeutPun$fit - (critval * predsNeutPun$se.fit)
vanuatu_neut_fit <- NeutPunModFS$family$linkinv(predsNeutPun$fit)
vanuatu_neut_lwr <- NeutPunModFS$family$linkinv(lwr_neut)
vanuatu_neut_upr <- NeutPunModFS$family$linkinv(upr_neut)

# Build the final cross-cultural results table
cross_cultural_table <- data.frame(
  Actor = c("Antisocial", "Antisocial", "Neutral", "Neutral"),
  Condition = c("Anonymous", "In person", "Anonymous", "In person"),
  Vanuatu_Proportion = c(vanuatu_anti_fit[1], vanuatu_anti_fit[2], vanuatu_neut_fit[1], vanuatu_neut_fit[2]),
  Vanuatu_CI_Lower = c(vanuatu_anti_lwr[1], vanuatu_anti_lwr[2], vanuatu_neut_lwr[1], vanuatu_neut_lwr[2]),
  Vanuatu_CI_Upper = c(vanuatu_anti_upr[1], vanuatu_anti_upr[2], vanuatu_neut_upr[1], vanuatu_neut_upr[2]),
  Sweden_Proportion = c(sweden_anon_anti[1], sweden_public_anti[1], sweden_anon_neut[1], sweden_public_neut[1]),
  Sweden_CI_Lower = c(sweden_anon_anti[2], sweden_public_anti[2], sweden_anon_neut[2], sweden_public_neut[2]),
  Sweden_CI_Upper = c(sweden_anon_anti[3], sweden_public_anti[3], sweden_anon_neut[3], sweden_public_neut[3])
)

cat("\n=== CROSS-CULTURAL COMPARISON TABLE (SWEDEN VS VANUATU) ===\n")
print(cross_cultural_table, row.names = FALSE)



