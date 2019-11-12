#load("clean.data/meadows-2018.rda")
#load("clean.data/traits-2018-61.rda")
load("clean.data/meadows-df.rda")

source("R/plot.effect.r")

#---#

dfa <- data.frame(mpd.z = tapply(df$mpd.z, df$id, mean), 
                  mpd.a.z = unname(tapply(df$mpd.a.z, df$id, mean)),
                  mntd.z = unname(tapply(df$mntd.z, df$id, mean)),
                  mntd.a.z = unname(tapply(df$mntd.a.z, df$id, mean)))


dfa <- sapply(df[, 5:50], function(x) tapply(x, df$id, mean, na.rm = T))
dfa <- data.frame(dfa)
dfa$site <- sapply(strsplit(as.character(rownames(dfa)), "-"), function(x) x[[1]])
dfa$site <- factor(dfa$site)

boxplot.effect("cwm.hgt")

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("cwm.hgt")
boxplot.effect("cwm.hgt")
plot.effect("cwm.seed")
boxplot.effect("cwm.seed")
par(op)

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("cwm.sla")
boxplot.effect("cwm.sla")
plot.effect("cwm.ldmc")
boxplot.effect("cwm.ldmc")
par(op)

#---trait CWMs---#

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("cwm.hgt")
plot.effect("cwm.sla")
plot.effect("cwm.ldmc")
plot.effect("cwm.seed")
par(op)

#---trait FRic, FEve and MPDa---#

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

#---quantitative (5) and qualitative (6) traits FRic, FEve, FDiv and MPDa---#

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

#---single traits SES MPD and SES MPDa---#

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

#---quantitative (5) and qualitative (6) traits SES MPD and SES MPDa---#

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("mpd.z.5")
plot.effect("mpd.a.z.5")
plot.effect("mpd.z.6")
plot.effect("mpd.a.z.6")
par(op)

#---phylogenetic SES MPD and SES MPDa---#

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("mpd.z.phylo")
plot.effect("mpd.a.z.phylo")
plot.effect("mpd.z.phylo.ao")
plot.effect("mpd.a.z.phylo.ao")
par(op)