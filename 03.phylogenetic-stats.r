load("clean.data/meadows-2018.rda")
load("clean.data/meadows-phylo-ses-ao.rda")

meta <- meta2

plot.year <- function(y , ind )
{
  dd <- meta[[ind]][meta$year == y]
  ss <- factor(meta$site[meta$year == y], levels = paste0(c("a", "f1", "f2", "m"), "-", substr(y, 3, 4)))
  boxplot(dd ~ ss, col = c("tomato", "steelblue", "skyblue", "limegreen"), ylab = ind)
}

plot.year(2014, "mpd.z")


op <- par(mfrow = c(2, 5),  mar = c(2, 4, 0.5, 0.5))

plot.year(2014, "mpd.z")
plot.year(2015, "mpd.z")
plot.year(2016, "mpd.z")
plot.year(2017, "mpd.z")
plot.year(2018, "mpd.z")

plot.year(2014, "mpd.a.z")
plot.year(2015, "mpd.a.z")
plot.year(2016, "mpd.a.z")
plot.year(2017, "mpd.a.z")
plot.year(2018, "mpd.a.z")

par(op)

#---#


df <- data.frame(mpd.z = meta$mpd.z, mpd.a.z = meta$mpd.a.z, mntd.z = meta$mntd.z, mntd.a.z = meta$mntd.a.z,
                 year = meta$year, 
                 site = sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)

#---------------------------------------------------------#

fit <- lm(mpd.z ~ year + id, data = df)
( an <- anova(fit) )

eta <- an[, 2] / sum(an[, 2])
names(eta) <- row.names(an)
eta * 100


fit <- lm(mpd.a.z ~ year + id, data = df)
( an <- anova(fit) )

eta <- an[, 2] / sum(an[, 2])
names(eta) <- row.names(an)
eta * 100



fit <- lm(mntd.z ~ year + id, data = df)
( an <- anova(fit) )

eta <- an[, 2] / sum(an[, 2])
names(eta) <- row.names(an)
eta * 100


fit <- lm(mntd.a.z ~ year + id, data = df)
( an <- anova(fit) )

eta <- an[, 2] / sum(an[, 2])
names(eta) <- row.names(an)
eta * 100


#---------------------------------------------------------#


col4 <- c("tomato", "steelblue", "skyblue", "limegreen")

op <- par(mfcol = c(2, 2))

boxplot(mpd.z ~ site, dfa, col = col4, ylab = "SES MPD")
boxplot(mpd.a.z ~ site, dfa, col = col4, ylab = expression("SES "*MPD[a]))

boxplot(mntd.z ~ site, dfa, col = col4, ylab = "SES MNTD")
boxplot(mntd.a.z ~ site, dfa, col = col4, ylab = expression("SES "*MNTD[a]))

par(op)

#---------------------------------------------------------#

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


fit <- lm(mpd.z ~ site, dfa)
anova(fit)

res <- pairwise.t.test(dfa$mpd.z, dfa$site)
abc <- cld(res)

pp <- tapply(dfa$mpd.z, dfa$site, function(x) t.test(x)$p.value)

abc[pp < 0.05] <- paste0(abc[pp < 0.05], "(*)")


boxplot(mpd.z ~ site, dfa, col = col4, ylab = "SES MPD", ylim = c(-.5, .5), axes = F, cex.lab = 1)
axis(1, at = 1:4, cex.axis = 1, labels = c("поле", "залежь 1", "залежь 2", "луг"))
axis(2, cex.axis = 1)
box()

text(x = 1:4, y = 0.45, labels = abc, cex = 1)

#---#

fit <- lm(mpd.a.z ~ site, dfa)
anova(fit)

res <- pairwise.t.test(dfa$mpd.a.z, dfa$site)
abc <- cld(res)

pp <- tapply(dfa$mpd.a.z, dfa$site, function(x) t.test(x)$p.value)

abc[pp < 0.05] <- paste0(abc[pp < 0.05], "(*)")


boxplot(mpd.a.z ~ site, dfa, col = col4, ylab = expression("SES "*MPD[a]), ylim = c(-.5, .5), axes = F, cex.lab = 1)
axis(1, at = 1:4, cex.axis = 1, labels = c("поле", "залежь 1", "залежь 2", "луг"))
axis(2, cex.axis = 1)
box()

text(x = 1:4, y = 0.45, labels = abc, cex = 1)




#________________________________________________________________________________________

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


png("figures/ses-phylo-angiosperm-only.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()



#________________________________________________________________________________________

load("clean.data/meadows-phylo-ses.rda")


df <- data.frame(mpd.z = meta$mpd.z, mpd.a.z = meta$mpd.a.z, mntd.z = meta$mntd.z, mntd.a.z = meta$mntd.a.z,
                 year = meta$year, 
                 site = sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]),
                 id = paste0(sapply(strsplit(as.character(meta$site), "-"), function(x) x[[1]]), "-", c(rep(1:50, 2), rep(1:65, 15))))
df$numsite <- as.numeric(df$site)


dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))

dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)



png("figures/ses-phylo-all-species.png", width = 1000, height = 1000)

op <- par(mfcol = c(2,2), mar = c(4, 4.1, 1, 1), cex = 1.5)

boxplotm(dfa$mpd.z, dfa$site, 1)
boxplotm(dfa$mpd.a.z, dfa$site, 2)
boxplotm(dfa$mntd.z, dfa$site, 3)
boxplotm(dfa$mntd.a.z, dfa$site, 4)

par(op)

dev.off()



