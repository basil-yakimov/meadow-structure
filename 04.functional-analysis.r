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

#__________________________________________________________________________________

library(picante)

meta.d <- meta

ses.mpd <- ses.mpd(dd, dist, null.model = "independentswap", runs = 999, iterations = 1000)
meta.d$mpd.z <- ses.mpd$mpd.obs.z
meta.d$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.d$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.d$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist, null.model = "independentswap", runs = 999, iterations = 1000)
meta.d$mntd.z <- ses.mntd$mntd.obs.z
meta.d$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.d$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.d$mntd.a.p <- ses.mntd.a$mntd.obs.p



meta.d15 <- meta

ses.mpd <- ses.mpd(dd, dist15, null.model = "independentswap", runs = 999, iterations = 1000)
meta.d15$mpd.z <- ses.mpd$mpd.obs.z
meta.d15$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist15, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.d15$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.d15$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist15, null.model = "independentswap", runs = 999, iterations = 1000)
meta.d15$mntd.z <- ses.mntd$mntd.obs.z
meta.d15$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist15, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.d15$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.d15$mntd.a.p <- ses.mntd.a$mntd.obs.p



meta.d611 <- meta

ses.mpd <- ses.mpd(dd, dist611, null.model = "independentswap", runs = 999, iterations = 1000)
meta.d611$mpd.z <- ses.mpd$mpd.obs.z
meta.d611$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist611, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.d611$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.d611$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist611, null.model = "independentswap", runs = 999, iterations = 1000)
meta.d611$mntd.z <- ses.mntd$mntd.obs.z
meta.d611$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist611, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.d611$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.d611$mntd.a.p <- ses.mntd.a$mntd.obs.p

save(meta.d, meta.d15, meta.d611, file = "clean.data/meadows-trait-ses.rda")

