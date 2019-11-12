plot.effect <- function(var, frame = df)
{
  require(nlme)
  require(effects)
  require(multcomp)
  
  assign(".frame", frame, env = .GlobalEnv)
  f <- as.formula(paste0(var, " ~ site + year"))
  assign(".form", f, env=.GlobalEnv)
  fit <- lme(fixed = .form, random = ~ 1 | id, data = .frame, method = "REML", na.action = "na.omit")
  
  an <- anova(fit)
  ef <- effect("site", fit)
  let <- cld(summary(glht(fit, linfct = mcp(site = "Tukey"))))$mcletters$Letters
  
  bp <- barplot(as.vector(ef$fit), col = c("tomato", "steelblue", "skyblue", "limegreen"), 
                ylim = c(min(ef$lower), max(ef$upper) + (max(ef$upper) - min(ef$lower)) * 0.2 ),
                ylab = var)
  arrows(x0 = bp, y0 = ef$lower,
         y1 = ef$upper, code = 3, 	angle = 90)
  text(x = bp, y = ef$upper,
       labels = as.vector(let),
       pos = 3)
  
  legend("top", legend = "", bty = "n", title = paste0("ANOVA p-values: site - ", round(an[2, 4], 3),
                                                       ", year - ", round(an[3, 4], 3)))
  #remove(.frame, env=.GlobalEnv)
}

boxplot.effect <- function(var, frame = dfa)
{
  cld <- function(res)
  {
    require(multcompView)
    require(reshape)
    
    a <- melt(res$p.value)
    a.cc <- na.omit(a)
    a.pvals <- a.cc[, 3]
    names(a.pvals) <- paste(a.cc[, 2], a.cc[, 1], sep="-")
    multcompLetters(a.pvals)$Letters
  }
  
  res <- pairwise.t.test(frame[, var], frame$site)
  abc <- cld(res)
  
  pp <- tapply(frame[, var], frame$site, function(x) t.test(x)$p.value)
  
  abc[pp < 0.05] <- paste0(abc[pp < 0.05], "(*)")
  
  ylim <- c(min(frame[, var]), max(frame[, var]) +(max(frame[, var]) - min(frame[, var]))/10)
  col4 <- c("tomato", "steelblue", "skyblue", "limegreen")
  boxplot(frame[, var] ~ frame$site, col = col4, ylab = var, ylim = ylim, axes = F, cex.lab = 1)
  axis(1, at = 1:4, cex.axis = 1, labels = c("agro", "fallow 1", "fallow 2", "meadow"))
  axis(2, cex.axis = 1)
  box()
  
  text(x = 1:4, y = max(frame[, var]) +(max(frame[, var]) - min(frame[, var]))/12, labels = abc, cex = 1)
}