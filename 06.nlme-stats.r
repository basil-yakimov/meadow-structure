load("clean.data/meadows-2018.rda")
load("clean.data/traits-2018-61.rda")

df <- data.frame(year = meta$year, 
                 site = sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dat <- dat[, row.names(tr)]
pp <- dat / rowSums(dat)


df$cwm.hgt <- as.matrix(pp) %*% tr$hgt
df$cwm.sla <- as.matrix(pp) %*% tr$SLA
df$cwm.ldmc <- as.matrix(pp) %*% tr$LDMC
df$cwm.seed <- as.matrix(pp[, !is.na(tr$seed_mass)]) %*% tr$seed_mass[!is.na(tr$seed_mass)]

source("R/plot.effect.r")

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("cwm.hgt")
plot.effect("cwm.sla")
plot.effect("cwm.ldmc")
plot.effect("cwm.seed")
par(op)


#---#

load("clean.data/meadows-single-trait-ses.rda")

df$mpd.z.hgt <- meta.hgt$mpd.z
df$mpd.a.z.hgt <- meta.hgt$mpd.a.z
df$mpd.z.sla <- meta.SLA$mpd.z
df$mpd.a.z.sla <- meta.SLA$mpd.a.z
df$mpd.z.ldmc <- meta.LDMC$mpd.z
df$mpd.a.z.ldmc <- meta.LDMC$mpd.a.z
df$mpd.z.seed <- meta.seedmass$mpd.z
df$mpd.a.z.seed <- meta.seedmass$mpd.a.z

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("mpd.z.hgt")
plot.effect("mpd.z.sla")
plot.effect("mpd.z.ldmc")
plot.effect("mpd.z.seed")
par(op)


op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("mpd.a.z.hgt")
plot.effect("mpd.a.z.sla")
plot.effect("mpd.a.z.ldmc")
plot.effect("mpd.a.z.seed")
par(op)

#---#

library(FD)

vec <- tr$hgt
names(vec) <- row.names(tr)
res <- dbFD(vec, pp, corr = "cailliez")

df$fric.hgt <- res$FRic
df$feve.hgt <- res$FEve
df$mpd.a.hgt <- res$RaoQ

vec <- tr$SLA
names(vec) <- row.names(tr)
res <- dbFD(vec, pp, corr = "cailliez")

df$fric.sla <- res$FRic
df$feve.sla <- res$FEve
df$mpd.a.sla <- res$RaoQ

vec <- tr$LDMC
names(vec) <- row.names(tr)
res <- dbFD(vec, pp, corr = "cailliez")

df$fric.ldmc <- res$FRic
df$feve.ldmc <- res$FEve
df$mpd.a.ldmc <- res$RaoQ

vec <- tr$seed_mass
names(vec) <- row.names(tr)
res <- dbFD(vec, pp, corr = "cailliez")

df$fric.seed <- res$FRic
df$feve.seed <- res$FEve
df$mpd.a.seed <- res$RaoQ

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("fric.hgt")
plot.effect("fric.sla")
plot.effect("fric.ldmc")
plot.effect("fric.seed")
par(op)


op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("feve.hgt")
plot.effect("feve.sla")
plot.effect("feve.ldmc")
plot.effect("feve.seed")
par(op)

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("mpd.a.hgt")
plot.effect("mpd.a.sla")
plot.effect("mpd.a.ldmc")
plot.effect("mpd.a.seed")
par(op)

#---#

res <- dbFD(tr[, 1:5], pp, corr = "cailliez")

df$fric.5 <- res$FRic
df$feve.5 <- res$FEve
df$fdiv.5 <- res$FDiv
df$mpd.a.5 <- res$RaoQ

tr$age_first_flowering <- factor(tr$age_first_flowering)
tr$branching <- factor(tr$branching)
tr$leaf_distribution <- factor(tr$leaf_distribution)
tr$growth_form <- factor(tr$growth_form)
tr$dispersal <- factor(tr$dispersal)
tr$life_span <- factor(tr$life_span)
library(cluster)
dist <- daisy(tr[, 6:11], metric = "gower")
res <- dbFD(x = dist, a = pp, corr = "cailliez")

df$fric.6 <- res$FRic
df$feve.6 <- res$FEve
df$fdiv.6 <- res$FDiv
df$mpd.a.6 <- res$RaoQ

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("fric.5")
plot.effect("feve.5")
plot.effect("fdiv.5")
plot.effect("mpd.a.5")
par(op)

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("fric.6")
plot.effect("feve.6")
plot.effect("fdiv.6")
plot.effect("mpd.a.6")
par(op)

#---#

