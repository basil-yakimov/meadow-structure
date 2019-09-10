
library(V.PhyloMaker)
library(dplyr)
library(readxl)

load("clean.data/meadows-2018.rda")
rm(meta)

gf <- tips.info[, c(3,4)] %>% group_by(family, genus) %>% slice(1) %>% ungroup()

md.sp.list <- colnames(dat)
md.g.list <- sapply(strsplit(md.sp.list, "_"), function(x) paste(x[1]))

id2 <- which(md.sp.list %in% GBOTB.extended$tip.label)
md.sp.list[-id2]
md.sp.list[112] <- "Silene_flos-cuculi"

md.f.list <- c()
for (i in 1:length(md.g.list)) {
  if (md.g.list[i] %in% gf$genus){
    md.f.list <- c(md.f.list, gf$family[which(gf$genus %in% md.g.list[i])])
  } else {
    md.f.list <- c(md.f.list, NA)
  }
}

md.f.list[105:107] <- "Polygonaceae"

md.df <- data.frame(species = md.sp.list, genus = md.g.list, family = md.f.list)
md.df[, c("species.relative", "genus.relative")] <- NA

rm(id2, md.sp.list, md.g.list, md.f.list)

md.tree <- phylo.maker(sp.list = md.df, scenarios="S1")
md.tree <- md.tree$scenario.1

save(md.tree, file = "clean.data/meadows-phylo.rda")
