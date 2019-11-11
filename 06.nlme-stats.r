#load("clean.data/meadows-2018.rda")
#load("clean.data/traits-2018-61.rda")
load("clean.data/meadows-df.rda")

source("R/plot.effect.r")

op <- par(mfcol = c(2, 2), mar = c(0.5, 4, 0.5, 0.5))
plot.effect("cwm.hgt")
plot.effect("cwm.sla")
plot.effect("cwm.ldmc")
plot.effect("cwm.seed")
par(op)

#---#

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

