load("clean.data/meadows-2018.rda")
load("clean.data/traits-2018-61.rda")
load("clean.data/meadows-single-trait-ses.rda")


library(multcompView)
library(reshape)

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

boxplotm <- function(x, y, i){
  res <- pairwise.t.test(x, y)
  abc <- cld(res)
  
  pp <- tapply(x, y, function(x) t.test(x)$p.value)
  
  abc[pp < 0.05] <- paste0(abc[pp < 0.05], "(*)")
  
  col4 <- c("tomato", "steelblue", "skyblue", "limegreen")
  ylab <- c("SES MPD", expression("SES "*MPD[a]), "SES MNTD", expression("SES "*MNTD[a]))
  boxplot(x ~ y, col = col4, ylab = ylab[i], ylim = c(min(x), max(x)+0.21), axes = F, cex.lab = 1)
  axis(1, at = 1:4, cex.axis = 1, labels = c("поле", "залежь 1", "залежь 2", "луг"))
  axis(2, cex.axis = 1)
  box()
  
  text(x = 1:4, y = max(x) + 0.15, labels = abc, cex = 1)
}


#---#

df <- data.frame(mpd.z = meta.hgt$mpd.z, mpd.a.z = meta.hgt$mpd.a.z, mntd.z = meta.hgt$mntd.z, mntd.a.z = meta.hgt$mntd.a.z,
                 year = meta.hgt$year, 
                 site = sapply(strsplit(as.character(meta.hgt$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.hgt$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-height.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()


#---------------------------------------------------------#

df <- data.frame(mpd.z = meta.LA$mpd.z, mpd.a.z = meta.LA$mpd.a.z, mntd.z = meta.LA$mntd.z, mntd.a.z = meta.LA$mntd.a.z,
                 year = meta.LA$year, 
                 site = sapply(strsplit(as.character(meta.LA$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.LA$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-LA.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()


#---------------------------------------------------------#

df <- data.frame(mpd.z = meta.SLA$mpd.z, mpd.a.z = meta.SLA$mpd.a.z, mntd.z = meta.SLA$mntd.z, mntd.a.z = meta.SLA$mntd.a.z,
                 year = meta.SLA$year, 
                 site = sapply(strsplit(as.character(meta.SLA$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.SLA$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-SLA.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()


#---------------------------------------------------------#

df <- data.frame(mpd.z = meta.LDMC$mpd.z, mpd.a.z = meta.LDMC$mpd.a.z, mntd.z = meta.LDMC$mntd.z, mntd.a.z = meta.LDMC$mntd.a.z,
                 year = meta.LDMC$year, 
                 site = sapply(strsplit(as.character(meta.LDMC$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.LDMC$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-LDMC.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()


#---------------------------------------------------------#

df <- data.frame(mpd.z = meta.seedmass$mpd.z, mpd.a.z = meta.seedmass$mpd.a.z, mntd.z = meta.seedmass$mntd.z, mntd.a.z = meta.seedmass$mntd.a.z,
                 year = meta.seedmass$year, 
                 site = sapply(strsplit(as.character(meta.seedmass$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.seedmass$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-seed-mass.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()
