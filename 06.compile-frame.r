load("clean.data/meadows-2018.rda")
load("clean.data/traits-2018-61.rda")

df <- data.frame(year = meta$year, 
                 site = sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dat <- dat[, row.names(tr)]
pp <- dat / rowSums(dat)
library(vegan)
dh <- decostand(dat, method = "hell")
pca <- rda(dh)

stage3 <- factor(substr(as.character(df$site), 1, 1), labels = c("a", "f", "m"))
col3 <- c("pink", "skyblue", "limegreen")
lab3 <- c("field", "fallow", "meadow")

ordiplot(pca, type = "n", scaling = 2)
ordispider(pca, groups = stage3, show.groups = levels(stage3), col = col3, scaling = 2)
legend("topleft", lab3, fill = col3)

wa <- scores(pca, display = "species", scaling = 2)
top10 <- which(rank(sqrt(rowSums(wa^2))) > 51)
arrows(0, 0, wa[top10, 1], wa[top10, 2])
text(wa[top10, ], labels = row.names(wa)[top10])

pca$CA$eig[1:5]/sum(pca$CA$eig)*100

df$pc1 <- scores(pca, display = "sites", scaling = 2)[, 1]
df$pc2 <- scores(pca, display = "sites", scaling = 2)[, 2]

#---qualitative trait cwm---#

df$cwm.hgt <- as.vector(as.matrix(pp) %*% tr$hgt)
df$cwm.sla <- as.vector(as.matrix(pp) %*% tr$SLA)
df$cwm.ldmc <- as.vector(as.matrix(pp) %*% tr$LDMC)
df$cwm.seed <- as.vector(as.matrix(pp[, !is.na(tr$seed_mass)]) %*% tr$seed_mass[!is.na(tr$seed_mass)])

#---single trait diversity---#

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

#---single-trait-ses---#

load("clean.data/meadows-single-trait-ses.rda")

df$mpd.z.hgt <- meta.hgt$mpd.z
df$mpd.a.z.hgt <- meta.hgt$mpd.a.z
df$mpd.z.sla <- meta.SLA$mpd.z
df$mpd.a.z.sla <- meta.SLA$mpd.a.z
df$mpd.z.ldmc <- meta.LDMC$mpd.z
df$mpd.a.z.ldmc <- meta.LDMC$mpd.a.z
df$mpd.z.seed <- meta.seedmass$mpd.z
df$mpd.a.z.seed <- meta.seedmass$mpd.a.z

#---qualitative (5) and quantitative (6) traits diversity---#

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

#---qualitative (5) and quantitative (6) traits ses---#

load("clean.data/meadows-trait-ses.rda")

df$mpd.z.5 <- meta.d15$mpd.z
df$mpd.a.z.5 <- meta.d15$mpd.a.z

df$mpd.z.6 <- meta.d611$mpd.z
df$mpd.a.z.6 <- meta.d611$mpd.a.z

#---all species and angionsperm-only phylogenetic ses---#

load("clean.data/meadows-phylo-ses.rda")

df$mpd.z.phylo <- meta$mpd.z
df$mpd.a.z.phylo <- meta$mpd.a.z
df$mntd.z.phylo <- meta$mntd.z
df$mntd.a.z.phylo <- meta$mntd.a.z

load("clean.data/meadows-phylo-ses-ao.rda")

df$mpd.z.phylo.ao <- meta$mpd.z
df$mpd.a.z.phylo.ao <- meta$mpd.a.z
df$mntd.z.phylo.ao <- meta$mntd.z
df$mntd.a.z.phylo.ao <- meta$mntd.a.z

#---#

save(df, file = "clean.data/meadows-df.rda")
