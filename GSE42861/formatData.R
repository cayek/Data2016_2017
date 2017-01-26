# load data send by OF
load("exp861.RData")
ls()

# save G and X
G <- t(expmat861)
## G
rm(expmat861)
dim(G)
saveRDS(G, "betanormalized_metylationlvl.rds")

# we scale and center data
X <- data.frame(disease.state = as.numeric(disease.state),
                age = as.numeric(age),
                gender = as.numeric(gender),
                smocking.status = as.numeric(smocking.status))
X <- scale(X)
X <- as.matrix(X)
rownames(X) <- rownames(G)
saveRDS(X, "X.rds")

# downsample for test
sample.row <- sample.int(nrow(G), size = 100)
sample.col <- sample.int(ncol(G), size = 2000)
saveRDS(G[sample.row, sample.col], "betanormalized_metylationlvl.sample.rds")
saveRDS(X[sample.row,], "X.sample.rds")
