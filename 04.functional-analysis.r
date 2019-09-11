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

