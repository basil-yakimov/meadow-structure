load("clean.data/meadows-2018.rda")
load("clean.data/traits-2018-61.rda")
load("clean.data/meadows-trait-ses.rda")



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

dd <- dat[, rownames(tr)]

df <- data.frame(mpd.z = meta.d$mpd.z, mpd.a.z = meta.d$mpd.a.z, mntd.z = meta.d$mntd.z, mntd.a.z = meta.d$mntd.a.z,
                 year = meta.d$year, 
                 site = sapply(strsplit(as.character(meta.d$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.d$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-all-traits.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()


#---------------------------------------------------------#

dd <- dat[, rownames(tr[1:5])]

df <- data.frame(mpd.z = meta.d15$mpd.z, mpd.a.z = meta.d15$mpd.a.z, mntd.z = meta.d15$mntd.z, mntd.a.z = meta.d15$mntd.a.z,
                 year = meta.d15$year, 
                 site = sapply(strsplit(as.character(meta.d15$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.d15$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-quantitative-traits.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()


#---------------------------------------------------------#

dd <- dat[, rownames(tr[6:11])]

df <- data.frame(mpd.z = meta.d611$mpd.z, mpd.a.z = meta.d611$mpd.a.z, mntd.z = meta.d611$mntd.z, mntd.a.z = meta.d611$mntd.a.z,
                 year = meta.d611$year, 
                 site = sapply(strsplit(as.character(meta.d611$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta.d611$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-qualitative-traits.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()
