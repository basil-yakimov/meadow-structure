load("clean.data/meadows-2018.rda")
load("clean.data/traits-2018-61.rda")

tr$age_first_flowering <- factor(tr$age_first_flowering)
tr$branching <- factor(tr$branching)
tr$leaf_distribution <- factor(tr$leaf_distribution)
tr$growth_form <- factor(tr$growth_form)
tr$dispersal <- factor(tr$dispersal)
tr$life_span <- factor(tr$life_span)

tr$SLA <- tr$SLA/tr$LDMC
tr["Vicia_sativa", "SLA"] <- 11600.923 / 0.37 / 1000
tr["Vicia_cracca", "SLA"] <- 19261.982 / 0.84 / 1000
tr["Vicia_tetrasperma", "SLA"] <- 6015.802 / 0.26 / 1000

trsm <- tr[!is.na(tr$seed_mass), ]
dd <- dat[, rownames(tr)]
ddsm <- dat[, rownames(trsm)]

library(cluster)

dist.hgt <- as.matrix(dist(tr$hgt))
colnames(dist.hgt) <- rownames(dist.hgt) <- rownames(tr)

dist.LA <- as.matrix(dist(tr$LA))
colnames(dist.LA) <- rownames(dist.LA) <- rownames(tr)

dist.SLA <- as.matrix(dist(tr$SLA))
colnames(dist.SLA) <- rownames(dist.SLA) <- rownames(tr)

dist.LDMC <- as.matrix(dist(tr$LDMC))
colnames(dist.LDMC) <- rownames(dist.LDMC) <- rownames(tr)

dist.seedmass <- as.matrix(dist(trsm$seed_mass))
colnames(dist.seedmass) <- rownames(dist.seedmass) <- rownames(trsm)

#__________________________________________________________________________________

library(picante)

meta.hgt <- meta

ses.mpd <- ses.mpd(dd, dist.hgt, null.model = "independentswap", runs = 999, iterations = 1000)
meta.hgt$mpd.z <- ses.mpd$mpd.obs.z
meta.hgt$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist.hgt, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.hgt$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.hgt$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist.hgt, null.model = "independentswap", runs = 999, iterations = 1000)
meta.hgt$mntd.z <- ses.mntd$mntd.obs.z
meta.hgt$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist.hgt, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.hgt$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.hgt$mntd.a.p <- ses.mntd.a$mntd.obs.p


meta.LA <- meta

ses.mpd <- ses.mpd(dd, dist.LA, null.model = "independentswap", runs = 999, iterations = 1000)
meta.LA$mpd.z <- ses.mpd$mpd.obs.z
meta.LA$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist.LA, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.LA$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.LA$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist.LA, null.model = "independentswap", runs = 999, iterations = 1000)
meta.LA$mntd.z <- ses.mntd$mntd.obs.z
meta.LA$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist.LA, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.LA$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.LA$mntd.a.p <- ses.mntd.a$mntd.obs.p


meta.SLA <- meta

ses.mpd <- ses.mpd(dd, dist.SLA, null.model = "independentswap", runs = 999, iterations = 1000)
meta.SLA$mpd.z <- ses.mpd$mpd.obs.z
meta.SLA$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist.SLA, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.SLA$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.SLA$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist.SLA, null.model = "independentswap", runs = 999, iterations = 1000)
meta.SLA$mntd.z <- ses.mntd$mntd.obs.z
meta.SLA$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist.SLA, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.SLA$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.SLA$mntd.a.p <- ses.mntd.a$mntd.obs.p


meta.LDMC <- meta

ses.mpd <- ses.mpd(dd, dist.LDMC, null.model = "independentswap", runs = 999, iterations = 1000)
meta.LDMC$mpd.z <- ses.mpd$mpd.obs.z
meta.LDMC$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dd, dist.LDMC, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.LDMC$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.LDMC$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dd, dist.LDMC, null.model = "independentswap", runs = 999, iterations = 1000)
meta.LDMC$mntd.z <- ses.mntd$mntd.obs.z
meta.LDMC$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dd, dist.LDMC, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.LDMC$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.LDMC$mntd.a.p <- ses.mntd.a$mntd.obs.p


meta.seedmass <- meta

ses.mpd <- ses.mpd(ddsm, dist.seedmass, null.model = "independentswap", runs = 999, iterations = 1000)
meta.seedmass$mpd.z <- ses.mpd$mpd.obs.z
meta.seedmass$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ddsm, dist.seedmass, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.seedmass$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta.seedmass$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(ddsm, dist.seedmass, null.model = "independentswap", runs = 999, iterations = 1000)
meta.seedmass$mntd.z <- ses.mntd$mntd.obs.z
meta.seedmass$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(ddsm, dist.seedmass, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta.seedmass$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta.seedmass$mntd.a.p <- ses.mntd.a$mntd.obs.p

save(meta.hgt, meta.LA, meta.SLA, meta.LDMC, meta.seedmass, file = "clean.data/meadows-single-trait-ses.rda")
