# To run at the command line, replace this with a function that calls setwd()
# to the local directory where you have the data and R files.
setwd_vanuatu_3pp()

library(splines)
library(sandwich)
library(lmtest)

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
  bottomOfPoly <- bottomOfPoly[order(-bottomOfPoly$x),]
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

########################################### Output participant descriptives and age table for methods section

ur <- subset( d2, Location == "Urban")
rur <- subset( d2, Location == "Rural")

describeSample <- function( dat ) data.frame(
  n = nrow(dat), meanAge = round(mean(dat$Age), 2), sdAge = round(sd(dat$Age), 2),
  minAge = round(min(dat$Age), 2), maxAge = round(max(dat$Age), 2),
  girls = sum(dat$Sex == "F"), boys = sum(dat$Sex == "M") )

rbind( All = describeSample(d2), Urban = describeSample(ur), Rural = describeSample(rur) )

ages <- seq(4,10,by=1)
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

table(cut(rur$Age, breaks = c(ages,11), right = FALSE, include.lowest = TRUE))

########################################### Punish antisocial more?

# Each participant judged both actors, so these proportions are paired. The p value is
# McNemar's exact test, i.e. an exact binomial test on the discordant participants. The
# odds ratio compares the two marginal proportions, so its confidence interval uses
# cluster-robust standard errors (sandwich::vcovCL) clustered on participant.
punMcNemar <- function(dat) {
  nP <- nrow(dat)
  long <- data.frame( id = rep(seq_len(nP), 2), anti = rep(1:0, each = nP),
                      pun = c(dat$AntiPunNum, dat$NeutPunNum) )
  mod <- glm( pun ~ anti, data = long, family = "binomial" )
  ci <- exp( coefci( mod, parm = "anti", vcov. = vcovCL, cluster = ~id ) )
  data.frame( n = nP,
     antiOnly = sum( dat$AntiPunNum > dat$NeutPunNum ),
     neutOnly = sum( dat$AntiPunNum < dat$NeutPunNum ),
     OR = exp( coef(mod)[["anti"]] ), ci95l = ci[1], ci95u = ci[2],
     p = binom.test( sum( dat$AntiPunNum > dat$NeutPunNum ),
                     sum( dat$AntiPunNum != dat$NeutPunNum ), 0.5 )$p.value )
}

punTab <- rbind( table(d2$NeutPun), table(d2$AntiPun) )
punTab
round( punTab[,2] / ( punTab[,1] + punTab[,2] ), 2 )
punMcNemar(d2)

d2Age4 <- d2[d2$Age>=4&d2$Age<5,]
punTabAge4 <- rbind( table(d2Age4$NeutPun), table(d2Age4$AntiPun) )
punTabAge4
round( punTabAge4[,2] / ( punTabAge4[,1] + punTabAge4[,2] ), 2 )
punMcNemar(d2Age4)

d2Age5 <- d2[d2$Age>=5&d2$Age<6,]
punTabAge5 <- rbind( table(d2Age5$NeutPun), table(d2Age5$AntiPun) )
punTabAge5
round( punTabAge5[,2] / ( punTabAge5[,1] + punTabAge5[,2] ), 2 )
punMcNemar(d2Age5)

d2Age6Plus <- d2[d2$Age>=6,]
punTabAge6Plus <- rbind( table(d2Age6Plus$NeutPun), table(d2Age6Plus$AntiPun) )
punTabAge6Plus
round( punTabAge6Plus[,2] / ( punTabAge6Plus[,1] + punTabAge6Plus[,2] ), 2 )
punMcNemar(d2Age6Plus)

########################################### Allocation figure (Figure 2, and its supplementary version)

# Each allocation is modelled from the two experimental factors plus age, and
# plotted at the median age. Each column of the figure varies one experimental
# factor and holds the other at a reference level: the anonymity column is free
# of economic cost, the economic cost column is anonymous. Every bar is
# therefore a prediction for a design cell children were tested in, and the
# anonymous free-of-cost cell is common to both columns.

ageForFig <- median(d2$Age)
cat( "\nMedian age, at which figure predictions are made:", round(ageForFig, 3), "\n" )

allocMods <- lapply( behaviours, function(beh)
  glm(as.formula(paste(beh,"~Condition+Cost+Age")), data=d2, family="binomial") )
names(allocMods) <- behaviours

cellPred <- function( model, splitVar, splitLevel, holdVar, holdLevel ) {
  grid <- data.frame( Age = ageForFig )
  grid[[splitVar]] <- splitLevel
  grid[[holdVar]] <- holdLevel
  pred <- predict( model, newdata = grid, type = "link", se.fit = TRUE )
  data.frame( prop  = model$family$linkinv( pred$fit ),
              ci95l = model$family$linkinv( pred$fit - critval * pred$se.fit ),
              ci95u = model$family$linkinv( pred$fit + critval * pred$se.fit ) )
}

# Observed proportions for the same design cells, shown alongside the fit in the
# supplementary version of the figure. Note these pool over age, whereas the
# fitted bars hold it at the median.
rawPred <- function( responses ) {
  ci <- prop.test( sum(responses == "Y"), length(responses) )$conf.int
  data.frame( prop = mean(responses == "Y"), ci95l = ci[1], ci95u = ci[2] )
}

# One column per experimental factor, one row per allocation. Bars within a
# row run antisocial agent then neutral, each split by the column's factor.
axisLabelCex <- .8
allocColumns <- list(
  list( var = "Condition", levels = c("Private","Public"), title = "Anonymity",
        tickLabels = c("Anonymous","In person"), tickCex = axisLabelCex,
        holdVar = "Cost", holdLevel = "N" ),
  list( var = "Cost", levels = c("N","Y"), title = "Economic cost",
        tickLabels = c("Free","Econ. cost"), tickCex = axisLabelCex,
        holdVar = "Condition", holdLevel = "Private" )
)

allocRows <- list(
  list( label = "Bad sweets",  behavs = c("AntiPun","NeutPun"),     fill = "red",     rawFill = "#FF9999" ),
  list( label = "Empty box",   behavs = c("AntiEmpty","NeutEmpty"), fill = "grey60",  rawFill = "grey90" ),
  list( label = "Good sweets", behavs = c("AntiGood","NeutGood"),   fill = niceGreen, rawFill = "#99D699" )
)

stackRows <- function( pieces ) Reduce(function(f1, f2) rbind(f1, f2), pieces)

predTable <- function( column, predFn ) stackRows( lapply( behaviours, function(beh)
  stackRows( lapply( column$levels, function(lv)
    data.frame( behav = beh, level = lv, predFn(beh, lv) ) ) ) ) )

allocFitted <- lapply( allocColumns, function(column) predTable( column, function(beh, lv)
  cellPred( allocMods[[beh]], column$var, lv, column$holdVar, column$holdLevel ) ) )

allocRaw <- lapply( allocColumns, function(column) predTable( column, function(beh, lv)
  rawPred( d2[[beh]][ d2[[column$var]] == lv & d2[[column$holdVar]] == column$holdLevel ] ) ) )

for( j in seq_along(allocColumns) ) {
  cat( "\nPlotted values, split by", allocColumns[[j]]$title, "- holding",
       allocColumns[[j]]$holdVar, "at", allocColumns[[j]]$holdLevel, "\n" )
  print( data.frame( allocFitted[[j]][,1:2], round(allocFitted[[j]][,-(1:2)], 3) ), row.names = FALSE )
}

########################################### Drawing the allocation figure

barXs <- 1:4
barXLim <- c(0.5, 4.5)

panelBars <- function( tab, behavs ) stackRows( lapply( behavs, function(beh) tab[tab$behav == beh,] ) )

drawBars <- function( bars, centres, fill, halfWidth ) {
  rect( centres - halfWidth, 0, centres + halfWidth, bars$prop, col = fill, border = "black" )
  arrows( centres, bars$ci95l, centres, bars$ci95u, angle = 90, code = 3, length = .025 )
}

drawAllocPanel <- function( row, j, showRaw, showLegend, showYAxis ) {
  plot.new()
  plot.window( xlim = barXLim, ylim = c(0,1), xaxs = "i", yaxs = "i" )
  if( showRaw ) {
    drawBars( panelBars(allocFitted[[j]], row$behavs), barXs - .17, row$fill,    .16 )
    drawBars( panelBars(allocRaw[[j]],    row$behavs), barXs + .17, row$rawFill, .16 )
  } else {
    drawBars( panelBars(allocFitted[[j]], row$behavs), barXs, row$fill, .32 )
  }
  # Every panel shares the same 0-1 probability scale, so the axis is only
  # drawn on the leftmost column to avoid repeating it across the figure.
  if( showYAxis ) {
    axis( 2, at = seq(0,1,by=.1), labels = c("0", sprintf("%.1f", seq(.1,.9,by=.1)), "1"),
          las = 1, cex.axis = .7, tcl = -.3, mgp = c(3,.5,0) )
    mtext( "Probability", side = 2, line = 1.9, cex = .8 )
  }
  box()
  text( barXLim[1] + .25, .93, row$label, adj = c(0,1), font = 2, cex = panelTitleCex )
  if( showLegend )
    legend( "topright", inset = c(.015,.03), legend = c("Model fit","Raw data"),
            fill = c(row$fill, row$rawFill), bty = "n", cex = .8 )
}

# Base R has no nested-category x-axis, so it is drawn as two rows of labelled
# cells below the panels: the column's factor within agent.
drawNestedAxis <- function( column ) {
  plot.new()
  plot.window( xlim = barXLim, ylim = c(0,2), xaxs = "i", yaxs = "i" )
  axisRows <- list(
    list( yTop = 2, cellWidth = 1, labels = rep(column$tickLabels, 2), cex = column$tickCex ),
    list( yTop = 1, cellWidth = 2, labels = c("Antisocial agent","Neutral agent"), cex = column$tickCex )
  )
  for( axisRow in axisRows ) {
    edges <- seq( barXLim[1], barXLim[2], by = axisRow$cellWidth )
    segments( edges, axisRow$yTop, edges, axisRow$yTop - 1 )
    text( edges[-length(edges)] + axisRow$cellWidth/2,
          axisRow$yTop - .5 - strheight("A", cex = axisRow$cex)/2,
          axisRow$labels, cex = axisRow$cex, adj = c(.5, 0) )
  }
  segments( barXLim[1], 1, barXLim[2], 1 )
}

panelTitleCex <- 1
colMargins <- list( c(3, .6), c(.6, .6) )

drawAllocFigure <- function( showRaw ) {
  oldPar <- par( no.readonly = TRUE )
  on.exit( par(oldPar) )
  par( oma = c(0,0,1.8,0) )
  marIn <- sapply(colMargins, sum) * par("csi")
  plotWidthIn <- ( par("din")[1] - sum(marIn) ) / 2
  layout( matrix(1:8, ncol = 2), widths = plotWidthIn + marIn, heights = c(1, 1, 1, .36) )
  for( j in seq_along(allocColumns) ) {
    marLeftRight <- colMargins[[j]]
    for( i in seq_along(allocRows) ) {
      par( mar = c(.6, marLeftRight[1], .6, marLeftRight[2]), cex = 1 )
      # The legend sits in the economic-cost column: its bars run lower there
      # than in the anonymity column, leaving the legend more room to clear
      # the tallest (Good sweets) bars.
      drawAllocPanel( allocRows[[i]], j, showRaw, showLegend = showRaw && j == 2,
                      showYAxis = j == 1 )
      # Column titles sit in the outer margin so that all rows stay equal height.
      if( i == 1 ) mtext( allocColumns[[j]]$title, side = 3, outer = TRUE, line = .3,
                          font = 2, cex = panelTitleCex, at = mean(par("fig")[1:2]) )
    }
    par( mar = c(.2, marLeftRight[1], .2, marLeftRight[2]) )
    drawNestedAxis( allocColumns[[j]] )
  }
}

png("figure2.png", width = 7, height = 7.5, units = "in", res = 300)
drawAllocFigure( showRaw = FALSE )
dev.off()

png("supplementaryFigure.png", width = 7, height = 7.5, units = "in", res = 300)
drawAllocFigure( showRaw = TRUE )
dev.off()

######################################### Makes main model table

# Center age variable
meanAge <- mean(d2$Age)
d2$AgeCentred <- d2$Age - meanAge

# Fit all models and store in lists
mods <- list()
mods2 <- list()
for(beh in behaviours) {
  cat("\n\n", beh, "\n")
  
  # Fit simple model
  mods[[beh]] <- glm(as.formula(paste(beh, "~Condition+Cost+Sex+AgeCentred")), 
                     data=d2, family="binomial")
  printCoefmat(summary(mods[[beh]])$coefficients)
  
  # Fit interaction model
  mods2[[beh]] <- update(mods[[beh]], . ~ . + Condition:AgeCentred + Sex:AgeCentred)
  printCoefmat(summary(mods2[[beh]])$coefficients)
  
  # Test interaction terms
  print(anova(mods[[beh]], mods2[[beh]], test="Chisq")[5])
}

# Function to create OR table
create_OR_table <- function(model_names, main_label) {
  
  results_list <- list()
  
  for(beh in model_names) {
    model_simple <- mods[[beh]]  # Simple model for main effects
    model_interact <- mods2[[beh]]  # Interaction model for interactions
    
    # Get coefficients from SIMPLE model (rows 1-5)
    coef_summary_simple <- coef(summary(model_simple))
    coefs_log_simple <- coef_summary_simple[, "Estimate"]
    pvals_simple <- coef_summary_simple[, "Pr(>|z|)"]
    zvals_simple <- coef_summary_simple[, "z value"]
    ci_log_simple <- confint(model_simple)
    
    # Get coefficients from INTERACTION model (rows 7-8)
    coef_summary_interact <- coef(summary(model_interact))
    coefs_log_interact <- coef_summary_interact[, "Estimate"]
    pvals_interact <- coef_summary_interact[, "Pr(>|z|)"]
    zvals_interact <- coef_summary_interact[, "z value"]
    ci_log_interact <- confint(model_interact)
    
    # Convert to OR scale
    OR_simple <- exp(coefs_log_simple[1:5])
    ci_lower_OR_simple <- exp(ci_log_simple[1:5, 1])
    ci_upper_OR_simple <- exp(ci_log_simple[1:5, 2])
    
    OR_interact <- exp(coefs_log_interact[6:7])
    ci_lower_OR_interact <- exp(ci_log_interact[6:7, 1])
    ci_upper_OR_interact <- exp(ci_log_interact[6:7, 2])
    
    # Get interaction test p-value
    interaction_test_p <- anova(model_simple, model_interact, test="Chisq")$"Pr(>Chi)"[2]
    
    # Store results in table structure
    results_list[[beh]] <- data.frame(
      OR = c(OR_simple, NA, OR_interact),
      CI_lower = c(ci_lower_OR_simple, NA, ci_lower_OR_interact),
      CI_upper = c(ci_upper_OR_simple, NA, ci_upper_OR_interact),
      Z = c(zvals_simple[1:5], NA, zvals_interact[6:7]),
      P = c(pvals_simple, interaction_test_p, pvals_interact[6:7])
    )
  }
  
  # Combine into table
  combined <- cbind(results_list[[model_names[1]]], results_list[[model_names[2]]])
  
  # Set row names
  rownames(combined) <- c(
    "(Intercept)", "ConditionPublic", "CostY", "SexM", "AgeCentred",
    "Interaction Test", "Condition:Age", "Sex:Age"
  )
  
  base_names <- c("OR", "CI_lower", "CI_upper", "Z", "P")   # ADD Z HERE
  colnames(combined) <- c(
    paste(model_names[1], base_names, sep="_"),
    paste(model_names[2], base_names, sep="_")
  )
  
  return(combined)
}

# Format helper function
format_OR_with_CI <- function(OR, ci_lower, ci_upper, z_val, p_val) {
  if(is.na(OR)) {
    return(list(OR_CI = "---", Z = "---", P = sprintf("%.3f", p_val)))
  }
  OR_CI_string <- sprintf("%.2f [%.2f, %.2f]", OR, ci_lower, ci_upper)
  Z_string <- sprintf("%.2f", z_val)
  P_string <- sprintf("%.3f", p_val)
  return(list(OR_CI = OR_CI_string, Z = Z_string, P = P_string))
}

# === TABLE 1: PUNISHMENT ===

pun_table <- create_OR_table(c("AntiPun", "NeutPun"), "Punishment")
pun_formatted <- lapply(1:8, function(i) {
  list(
    Anti = format_OR_with_CI(pun_table[i, "AntiPun_OR"], pun_table[i, "AntiPun_CI_lower"],
                             pun_table[i, "AntiPun_CI_upper"], pun_table[i, "AntiPun_Z"],
                             pun_table[i, "AntiPun_P"]),
    Neut = format_OR_with_CI(pun_table[i, "NeutPun_OR"], pun_table[i, "NeutPun_CI_lower"],
                             pun_table[i, "NeutPun_CI_upper"], pun_table[i, "NeutPun_Z"],
                             pun_table[i, "NeutPun_P"])
  )
})

compact_pun <- data.frame(
  Variable    = rownames(pun_table),
  AntiPun_OR_CI = sapply(pun_formatted, function(x) x$Anti$OR_CI),
  AntiPun_Z     = sapply(pun_formatted, function(x) x$Anti$Z),
  AntiPun_P     = sapply(pun_formatted, function(x) x$Anti$P),
  NeutPun_OR_CI = sapply(pun_formatted, function(x) x$Neut$OR_CI),
  NeutPun_Z     = sapply(pun_formatted, function(x) x$Neut$Z),
  NeutPun_P     = sapply(pun_formatted, function(x) x$Neut$P)
)
print(compact_pun, row.names = FALSE)

# === TABLE 2: REWARD ===

good_table <- create_OR_table(c("AntiGood", "NeutGood"), "Reward")
good_formatted <- lapply(1:8, function(i) {
  list(
    Anti = format_OR_with_CI(good_table[i, "AntiGood_OR"], good_table[i, "AntiGood_CI_lower"],
                             good_table[i, "AntiGood_CI_upper"], good_table[i, "AntiGood_Z"],
                             good_table[i, "AntiGood_P"]),
    Neut = format_OR_with_CI(good_table[i, "NeutGood_OR"], good_table[i, "NeutGood_CI_lower"],
                             good_table[i, "NeutGood_CI_upper"], good_table[i, "NeutGood_Z"],
                             good_table[i, "NeutGood_P"])
  )
})

compact_good <- data.frame(
  Variable      = rownames(good_table),
  AntiGood_OR_CI = sapply(good_formatted, function(x) x$Anti$OR_CI),
  AntiGood_Z     = sapply(good_formatted, function(x) x$Anti$Z),
  AntiGood_P     = sapply(good_formatted, function(x) x$Anti$P),
  NeutGood_OR_CI = sapply(good_formatted, function(x) x$Neut$OR_CI),
  NeutGood_Z     = sapply(good_formatted, function(x) x$Neut$Z),
  NeutGood_P     = sapply(good_formatted, function(x) x$Neut$P)
)
print(compact_good, row.names = FALSE)

# Write to CSV
write.csv(compact_pun, "OR_table_punishment.csv", row.names = FALSE)
write.csv(compact_good, "OR_table_reward.csv", row.names = FALSE)

######################## Spline model plot of public/private against age

# Age-specific contrasts between the anonymity conditions, from the spline models
# below: is there an anonymity effect at each year of age?

# Design-matrix row for one age, differenced across the two anonymity
# conditions. Cost is set only for models that contain it.
anonGap <- function( mod, ag, cost = "N" ) {
  nd <- data.frame( AgeCentred = ag - meanAge,
                    Condition = factor( c("Private","Public"), levels = mod$xlevels$Condition ) )
  if( !is.null(mod$xlevels$Cost) ) nd$Cost <- factor( cost, levels = mod$xlevels$Cost )
  X <- model.matrix( delete.response(terms(mod)), nd )
  X["2",] - X["1",]
}

contrastTest <- function( mod, v, label ) {
  est <- sum( v * coef(mod) )
  se  <- sqrt( drop( t(v) %*% vcov(mod) %*% v ) )
  data.frame( contrast = label, OR = exp(est),
              ci95l = exp(est - critval*se), ci95u = exp(est + critval*se),
              z = est/se, p = 2*pnorm(-abs(est/se)) )
}

gapTable <- function( mod ) {
  out <- do.call( rbind,
    lapply( 4:10, function(ag) contrastTest( mod, anonGap(mod,ag), paste("age", ag) ) ) )
  out[,2:6] <- round( out[,2:6], 3 )
  out
}

png("figure3.png", width = 5, height = 5, units = "in", res = 300)

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

modKnotNoInt <- glm(AntiPun ~ Condition + ns(AgeCentred, knots=knotLoc) + Cost,
    data=d2, family="binomial")
cat("\nAntiPun: anonymity x age interaction block\n")
print(anova(modKnotNoInt, modKnot, test="Chisq"))

cat("\nAntiPun: anonymity contrast (OR in person vs anonymous), free of economic cost\n")
print(gapTable(modKnot), row.names = FALSE)

############################## Investigating rewarding of antisocial actor

cat("Rewarded antisocial actor:",
    round(100 * mean(d2$AntiGood[d2$NeutPun == "Y"] == "Y"), 1), "% (punished neutral) vs",
    round(100 * mean(d2$AntiGood[d2$NeutPun == "N"] == "Y"), 1), "% (didn't punish neutral)\n")


AntiGoodMod3 <- update(mods2[["AntiGood"]], . ~ . + NeutPun)
summary(mods2[["AntiGood"]])
summary(AntiGoodMod3)
coef_NeutPun <- coef(AntiGoodMod3)["NeutPunY"]
ci_log <- confint(AntiGoodMod3, parm = "NeutPunY")
OR_NeutPun <- exp(coef_NeutPun)
ci_OR <- exp(ci_log)
cat("NeutPunY: OR =", round(OR_NeutPun, 2), 
    "[95% CI:", round(ci_OR[1], 2), ",", round(ci_OR[2], 2), "]\n")
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
png("figure4a.png", width = 5, height = 5, units = "in", res = 300)
plotTwoGreenRed(AntiGoodModCondAge, AntiGoodModCondAge, preddat1, preddat2)
plotTwoPointsGreenRed('AntiGoodNum', 'AntiGoodNum', 'd2$Condition=="Private"', 'd2$Condition=="Public"')
mtext("(a)",at=4,line=.5) 
mtext("Age (years)", side=1, line=2) 
mtext("Proportion participants rewarding antisocial actor", side=2, line=2) 
legend(6, 1, legend=c("In person", "Anonymous"), col=c("Red", niceGreen), lty=c("solid", "longdash"), lwd=c(2, 2))
dev.off()

AntiGoodModCondAgeNoInt <- glm(AntiGood ~ Condition + ns(AgeCentred, knots=knotLoc),
                               data=d2, family="binomial")
cat("\nAntiGood: anonymity x age interaction block, spline age\n")
print(anova(AntiGoodModCondAgeNoInt, AntiGoodModCondAge, test="Chisq"))
cat("\nAntiGood: anonymity contrast (OR in person vs anonymous), spline model\n")
print(gapTable(AntiGoodModCondAge), row.names = FALSE)

  
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
png("figure4b.png", width = 5, height = 5, units = "in", res = 300)
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
  # Taken before tagging, so it is the distance to the child actually matched.
  nearest <- min(subset(d2UrbPrivN, tagged==FALSE)$diff)
  d2UrbPrivN$tagged[which(d2UrbPrivN$diff==nearest & d2UrbPrivN$tagged == FALSE)[1]]<<-TRUE
  diffList <<- append( diffList, nearest )
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

# Sensitivity power analysis for binomial GLMs — main effects only
# Uses pwr package. Assumes p0 = 0.5 (zero intercept / balanced baseline).
# Binary predictors: pwr.2p2n.test() with actual group sizes
# Age (continuous): pwr.r.test(), OR per year via r-to-logistic conversion

library(pwr)

# Group sizes for binary predictors
n_Condition1 <- sum(d2$Condition == "Private")
n_Condition2 <- sum(d2$Condition == "Public")
n_Cost1      <- sum(d2$Cost == "Y")
n_Cost2      <- sum(d2$Cost == "N")
n_Sex1       <- sum(d2$Sex == "M")
n_Sex2       <- sum(d2$Sex == "F")
n_total      <- nrow(d2)
age_sd       <- sd(d2$Age)

# ── Binary predictors ─────────────────────────────────────────────────────────
# pwr.2p2n.test returns Cohen's h; convert to OR assuming p0 = 0.5

binary_mde <- function(n1, n2) {
  h  <- pwr.2p2n.test(n1 = n1, n2 = n2, sig.level = 0.05, power = 0.80)$h
  p1 <- sin(asin(sqrt(0.5)) + h / 2)^2
  OR <- p1 / (1 - p1)   # since p0 = 0.5, OR = (p1/(1-p1)) / 1
  c(h = round(h, 3), p1 = round(p1, 3), OR = round(OR, 2))
}

cat("=== Minimum detectable effect at 80% power, alpha = 0.05, p0 = 0.5 ===\n\n")

cat("Condition (n =", n_Condition1, "/", n_Condition2, "):\n")
print(binary_mde(n_Condition1, n_Condition2))

cat("\nCost (n =", n_Cost1, "/", n_Cost2, "):\n")
print(binary_mde(n_Cost1, n_Cost2))

cat("\nSex (n =", n_Sex1, "/", n_Sex2, "):\n")
print(binary_mde(n_Sex1, n_Sex2))

# ── Age (continuous) ──────────────────────────────────────────────────────────
# pwr.r.test returns the minimum detectable correlation r, read here as the
# point-biserial correlation between age and the binary behaviour. Convert to
# Cohen's d, then to OR per SD via the logistic approximation, then to OR per year.

r      <- pwr.r.test(n = n_total, sig.level = 0.05, power = 0.80)$r
d      <- 2 * r / sqrt(1 - r^2)          # point-biserial r to Cohen's d
OR_sd  <- exp(d * pi / sqrt(3))          # OR per SD of age
OR_yr  <- OR_sd ^ (1 / age_sd)           # OR per year

cat("\nAge (n =", n_total, "):\n")
cat("r =", round(r, 3), " d =", round(d, 3), " OR per SD =", round(OR_sd, 3), "\n")
cat("OR per year =", round(OR_yr, 3), "\n")

#####################################
# To produce R Markdown for this file

# Create custom CSS file
# writeLines(".main-container {
#  max-width: 1130px !important;
# }", "path to custom.css")

# Render with custom CSS
# rmarkdown::render("path to analysis_vanuatu_3pp.R", 
#                  output_format = rmarkdown::html_document(
#                    pandoc_args = c("--metadata", "author=Anon For Peer Review"),
#                    css = "path to custom.css"
#                  ))