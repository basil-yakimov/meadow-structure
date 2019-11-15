source("R/assembly-processes.r")

load("clean.data/meadows-2018.rda")
load("clean.data/traits-2018-61.rda")

tr$age_first_flowering <- factor(tr$age_first_flowering)
tr$branching <- factor(tr$branching)
tr$leaf_distribution <- factor(tr$leaf_distribution)
tr$growth_form <- factor(tr$growth_form)
tr$dispersal <- factor(tr$dispersal)
tr$life_span <- factor(tr$life_span)


dd <- dat[, rownames(tr)]

library(cluster)

dist <- as.matrix(daisy(tr, metric = "gower"))
dist15 <- as.matrix(daisy(tr[, 1:5], metric = "gower"))
dist611 <- as.matrix(daisy(tr[, 6:11], metric = "gower"))


tp <- TotalPool.FDtest(1000, SpeciesMatrix = dd, TraitMatrix = tr[, 1:5], dist="gower")
rp <- RestrictedPool.FDtest(1000, SpeciesMatrix = dd, TraitMatrix = tr[, 1:5], dist="gower")
fc <- freqClasses.FDtest(1000, SpeciesMatrix = dd, TraitMatrix = tr[, 1:5], dist="gower")
fw <- freqWeights.FDtest(1000, SpeciesMatrix = dd, TraitMatrix = tr[, 1:5], dist="gower")

save(tp, rp, fc, fw, file = "clean.data/meadows-trait-r.rda")

load("clean.data/meadows-trait-r.rda")

tpq <- as.data.frame(tp$QAlphaFD)
rpq <- as.data.frame(rp$QAlphaFD)

plot(tpq$Q.SES, rpq$Q.SES, pch = 16, xlab = expression("SES of "*alpha*"-Rao"),
     ylab = expression("SES of "*alpha*"-Rao (restricted trait range)"))
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

