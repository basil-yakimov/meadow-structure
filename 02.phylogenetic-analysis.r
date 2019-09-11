library(ape)
library(picante)

load("clean.data/meadows-phylo.rda")
load("clean.data/meadows-2018.rda")

#all species

mat <- cophenetic(md.tree)
colnames(mat)[53] <- "Silene_flos.cuculi"
rownames(mat)[53] <- "Silene_flos.cuculi"

ses.mpd <- ses.mpd(dat, mat, null.model = "independentswap", runs = 999, iterations = 1000)
meta$mpd.z <- ses.mpd$mpd.obs.z
meta$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dat, mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dat, mat, null.model = "independentswap", runs = 999, iterations = 1000)
meta$mntd.z <- ses.mntd$mntd.obs.z
meta$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dat, mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)

meta$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta$mntd.a.p <- ses.mntd.a$mntd.obs.p

save(meta, file = "clean.data/meadows-phylo-ses.rda")

#_________________________________________________________________________________________

#angiosperm only

load("clean.data/meadows-phylo-angiosp-only.rda")
load("clean.data/meadows-2018-angiosp-only.rda")

meta2 <- meta

mat2 <- cophenetic(md.tree2)
colnames(mat2)[53] <- "Silene_flos.cuculi"
rownames(mat2)[53] <- "Silene_flos.cuculi"

ses.mpd <- ses.mpd(dat2, mat2, null.model = "independentswap", runs = 999, iterations = 1000)
meta2$mpd.z <- ses.mpd$mpd.obs.z
meta2$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(dat2, mat2, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta2$mpd.a.z <- ses.mpd.a$mpd.obs.z
meta2$mpd.a.p <- ses.mpd.a$mpd.obs.p

ses.mntd <- ses.mntd(dat2, mat2, null.model = "independentswap", runs = 999, iterations = 1000)
meta2$mntd.z <- ses.mntd$mntd.obs.z
meta2$mntd.p <- ses.mntd$mntd.obs.p

ses.mntd.a <- ses.mntd(dat2, mat2, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)

meta2$mntd.a.z <- ses.mntd.a$mntd.obs.z
meta2$mntd.a.p <- ses.mntd.a$mntd.obs.p

save(meta2, file = "clean.data/meadows-phylo-ses-ao.rda")
