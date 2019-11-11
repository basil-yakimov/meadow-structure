# 
# 
# plot.effect.2 <- function(frame, ylab = "cwm.sla")
# {
#   require(nlme)
#   require(effects)
#   require(multcomp)
#   
#   #df <- data.frame(var = var, year = year, site = site, id = id)
#   
#   #f <- as.formula(paste0(var, " ~ site + year"))
#   assign(".frame", frame, env = .GlobalEnv)
#   fit <- lme(fixed = var ~ site + year, random = ~ 1 | id, .frame, method = "REML")
#   
#   an <- anova(fit)
#   ef <- effect("site", fit)
#   let <- cld(summary(glht(fit, linfct = mcp(site = "Tukey"))))$mcletters$Letters
#   
#   bp <- barplot(as.vector(ef$fit), col = c("tomato", "steelblue", "skyblue", "limegreen"), 
#                 ylim = c(min(ef$lower), max(ef$upper) + (max(ef$upper) - min(ef$lower)) * 0.2 ),
#                 ylab = ylab)
#   arrows(x0 = bp, y0 = ef$lower,
#          y1 = ef$upper, code = 3, 	angle = 90)
#   text(x = bp, y = ef$upper,
#        labels = as.vector(let),
#        pos = 3)
#   
#   legend("top", legend = "", bty = "n", title = paste0("ANOVA p-values: site - ", round(an[2, 4], 3),
#                                                        ", year - ", round(an[3, 4], 3)))
# }
# 
# library(effects)
# fc <- function(dta, formula, terms) 
# {
#   print(m1 <- lm(formula, .dta))
#   Effect(terms, m1)
# }
# form <- prestige ~ income*type + education
# terms <- c("income", "type")
# fc(Duncan, form, terms)
# 
# fc.working <- function(dta, var, terms)
# {
#   assign(".dta", dta, env=.GlobalEnv)
#   
#   f <- as.formula(paste0(var, " ~ site + year"))
#   assign(".f", f, env=.GlobalEnv)
#   print(m1 <- lme(fixed = .f, random = ~ 1 | id, data = .frame, method = "REML"))
#   Effect(terms, m1)
#   plot(allEffects(m1))
#   remove(".dta", envir=.GlobalEnv)
# }
# fc.working(Duncan, form, terms)
# 
# 
# fc.working(df, "cwm.sla", "site")



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

# plot.effect(frame = df, "pc1")
# 
# 
# plot.int <- function(var = "cwm.sla")
# {
#   require(nlme)
#   require(effects)
# 
#   f <- as.formula(paste0(var, " ~ site * year"))
#   fit <- lme(fixed = f, random = ~ 1 | id, data = df, method = "REML")
#   plot(allEffects(fit))
#   }

# df <- data.frame(y = rnorm(90), x = gl(3, 30), b = factor(rep(1:30, 3)))
# 
# fit <- lme(fixed = y ~ x, random = ~ 1 | b, data = df, method = "REML")
# ef <- effect("x", fit)
# bp <- barplot(as.vector(ef$fit), col = c("tomato", "skyblue", "limegreen"),
#         ylim = c(min(ef$lower), max(ef$upper) + (max(ef$upper) - min(ef$lower)) * 0.2 ))
# arrows(x0 = bp, y0 = ef$lower, y1 = ef$upper, code = 3, angle = 90)
# 
# test1 <- function(y, x, b)
# {
#   
#   assign(".y", y, env = .GlobalEnv)
#   assign(".x", x, env = .GlobalEnv)
#   fit <- lme(fixed = .y ~ .x, random = ~ 1 | b, method = "REML")
#   print(summary(fit))
#   ef <- effect(".x", fit)
#   bp <- barplot(as.vector(ef$fit), col = c("tomato", "skyblue", "limegreen"),
#                 ylim = c(min(ef$lower), max(ef$upper) + (max(ef$upper) - min(ef$lower)) * 0.2 ))
#   arrows(x0 = bp, y0 = ef$lower, y1 = ef$upper, code = 3, angle = 90)
#   #remove(".y", env = .GlobalEnv)
#   #remove(".x", env = .GlobalEnv)
# }
# 
# test1(df$y, df$x, df$b)
# 
# test2 <- function(y, x, b)
# {
#   frame <- data.frame(y, x, b)
# 
#   fit <- lme(fixed = y ~ x, random = ~ 1 | b, frame, method = "REML")
#   ef <- effect("x", fit)
#   bp <- barplot(as.vector(ef$fit), col = c("tomato", "skyblue", "limegreen"),
#                 ylim = c(min(ef$lower), max(ef$upper) + (max(ef$upper) - min(ef$lower)) * 0.2 ))
#   arrows(x0 = bp, y0 = ef$lower, y1 = ef$upper, code = 3, angle = 90)
# }
# 
# test2(df$y, df$x, df$b)
# 
# 
# test3 <- function(df) {
#   fit <- lme(fixed = y ~ x, random = ~ 1 | b, data = df, method = "REML")
#   ef <- effect("x", fit)
#   bp <- barplot(as.vector(ef$fit), col = c("tomato", "skyblue", "limegreen"),
#                 ylim = c(min(ef$lower),
#                          max(ef$upper) + (max(ef$upper) - min(ef$lower)) * 0.2 ))
#   arrows(x0 = bp, y0 = ef$lower, y1 = ef$upper, code = 3, angle = 90)
# }
# 
# test3(data.frame(df$y, df$x, df$b))
