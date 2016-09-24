library(RPMM)
gset <- gset.combat[,colnames(gset.combat) %in% outmatrix$sample] # outmatrix has col "sample" with array IDs of interest
bval <- as.matrix(getBeta(gset[rownames(gset) %in% mvps,])) # mvp is list or vector of Most Variant Probes, to cluster on
rpmm <- blcTree(t(bval), verbose = 0,
                splitCriterion = blcSplitCriterionLevelWtdBIC)
rpmmClass <- blcTreeLeafClasses(rpmm)
a <- lapply(levels(rpmmClass), function(k) {colnames(bval)[which(rpmmClass == k)]}) # assign sample IDs from bval colnames
names(a) <- levels(rpmmClass)
m <- matrix(as.character(unlist(a)))

# clu matrix dimensions contingent on cluster count (ie. 4x clusters shown below)
clu <- matrix(c(rep(names(a)[[1]],
                    length(a[[1]])),rep(names(a)[[2]],length(a[[2]])),
                rep(names(a)[[3]], length(a[[3]])),
                rep(names(a)[[4]], length(a[[4]]))))

# construct dataframe of samples and clusters
asg <- data.frame(cbind(m, clu))
colnames(asg) <- c("Sample_name", "RPMM.cluster")
