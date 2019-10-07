plot.effect <- function(var)
{
  require(nlme)
  require(effects)
  
  f <- as.formula(paste0(var, "~ site + year"))
  fit <- lme(f, random = ~ 1 | id, df, method = "REML")
  
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
}
